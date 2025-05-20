import redis
import json
import openai
from workflows.celery_worker import celery

@celery.task(bind=True)
def ber_task(self, payload):
    """Generate a biological evaluation report for acute inhalation toxicity.

    The payload should contain a "payload" field with the user's message and
    optionally a "sid" to stream progress updates via redis pub/sub.
    """
    r = redis.Redis()
    sid = payload.get("sid")
    message = payload.get("payload", "")

    # notify start
    if sid:
        r.publish(
            "celery_updates",
            json.dumps({"sid": sid, "event": "progress", "data": {"status": "generating"}}),
        )

    try:
        response = openai.ChatCompletion.create(
            model="gpt-3.5-turbo",
            messages=[
                {
                    "role": "system",
                    "content": (
                        "You are a toxicology assistant. Identify the chemical of interest "
                        "from the user's request and generate a concise biological evaluation "
                        "report in Markdown focusing on acute inhalation toxicity."
                    ),
                },
                {"role": "user", "content": message},
            ],
            temperature=0.3,
        )
        report = response.choices[0].message["content"]

        if sid:
            r.publish(
                "celery_updates",
                json.dumps({"sid": sid, "event": "complete", "data": {"result": report}}),
            )
        return {"done": True, "report": report}
    except Exception as e:
        if sid:
            r.publish(
                "celery_updates",
                json.dumps({"sid": sid, "event": "error", "data": {"error": str(e)}}),
            )
        raise
