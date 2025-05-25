import redis, json, time
from workflows.celery_worker import celery
from webserver.model.message import MessageSchema

@celery.task(bind=True)
def probra_task(self, payload):
    r = redis.Redis()
    task_id = payload['task_id']
    for i in range(5):
        message = MessageSchema(role="assistant", content=f"Step {i}")
        event = { "type": "task_message", "data": message.model_dump(), "task_id": task_id }
        r.publish("celery_updates", json.dumps(event))
        time.sleep(1)
    
    result = {"task_id": task_id, "type": "task_file", "data": "{}"}
    r.publish("celery_updates", json.dumps(result))
    return {"done": True}
