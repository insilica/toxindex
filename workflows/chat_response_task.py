import redis, json
from workflows.celery_worker import celery
from webserver.model.message import MessageSchema

def build_prompt_from_conversation(conversation):
    prompt_parts = ["You are a helpful assistant. Answer the user's question.\n\n"]
    
    for raw_message in conversation:
        message = json.loads(raw_message)
        role = message.get("role", "user")
        content = message.get("content", "")
        prompt_parts.append(f"{role}: {content}\n")
        
    prompt_parts.append("\nassistant: ")
    return "".join(prompt_parts)

@celery.task(bind=True)
def chat_response_task(self, payload):
    r = redis.Redis()
    sid = payload.get("sid")
    task_id = payload.get("task_id")
    user_id = payload.get("user_id")
    conversation = payload.get("conversation", [])

    prompt = build_prompt_from_conversation(conversation)

    assistant_content = "Hello! How can I help you today?"

    response_message = MessageSchema(role="assistant", content=assistant_content)
    event = {
        "event": "task_message",
        "data": response_message.model_dump(),
        "sid": sid,
        "task_id": task_id
    }
    r.publish("celery_updates", json.dumps(event))

    return {"response": assistant_content}
