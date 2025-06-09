import uuid
from webserver.app import redis_listener

if __name__ == "__main__":
    thread_name = f"RedisListenerThread-{uuid.uuid4().hex[:8]}"
    redis_listener(thread_name) 