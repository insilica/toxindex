import redis, json, time
from workflows.celery_worker import celery

@celery.task(bind=True)
def probra_task(self, payload):
    r = redis.Redis()
    sid = payload.get("sid")
    for i in range(5):
        r.publish("celery_updates", json.dumps({
            "sid": sid,
            "event": "progress",
            "data": {"step": i}
        }))
        time.sleep(1)
    r.publish("celery_updates", json.dumps({
        "sid": sid,
        "event": "complete",
        "data": {"result": f"Processed: {payload['payload']}"}
    }))
    return {"done": True}
