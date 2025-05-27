import json
import redis
from workflows.celery_worker import celery


@celery.task(bind=True)
def interactive_echo_task(self, payload):
    """Example interactive task that echoes user messages.

    Conversation history is stored in Redis so the task can be stopped and
    later restarted from the same state. The task listens for messages on a
    Redis channel named ``task:{task_id}:input``. Send a JSON object with a
    ``content`` field to deliver a user message or ``{"command": "stop"}`` to
    stop the task.
    """
    r = redis.Redis()
    task_id = payload.get("task_id")
    if not task_id:
        raise ValueError("task_id is required")

    conv_key = f"task:{task_id}:conversation"
    status_key = f"task:{task_id}:status"
    input_channel = f"task:{task_id}:input"

    # Replay existing conversation
    stored = r.lrange(conv_key, 0, -1)
    for raw in stored:
        try:
            message = json.loads(raw)
        except Exception:
            continue
        r.publish(
            "celery_updates",
            json.dumps({"type": "task_message", "task_id": task_id, "data": message}),
        )

    r.set(status_key, "running")

    pubsub = r.pubsub()
    pubsub.subscribe(input_channel)
    try:
        for incoming in pubsub.listen():
            if incoming["type"] != "message":
                continue
            try:
                data = incoming["data"]
                if isinstance(data, bytes):
                    data = data.decode()
                data = json.loads(data)
            except Exception:
                continue

            if data.get("command") == "stop":
                r.set(status_key, "stopped")
                break

            content = data.get("content")
            if not content:
                continue

            user_msg = {"role": "user", "content": content}
            r.rpush(conv_key, json.dumps(user_msg))
            r.publish(
                "celery_updates",
                json.dumps({"type": "task_message", "task_id": task_id, "data": user_msg}),
            )

            assistant_msg = {"role": "assistant", "content": f"ECHO: {content}"}
            r.rpush(conv_key, json.dumps(assistant_msg))
            r.publish(
                "celery_updates",
                json.dumps({"type": "task_message", "task_id": task_id, "data": assistant_msg}),
            )
    finally:
        pubsub.close()

    return {"status": r.get(status_key).decode()}

