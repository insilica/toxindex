#!/usr/bin/env python3
"""
Standalone Redis listener for processing Celery events.
This runs in a separate GKE deployment to avoid duplicate listeners.
"""

import os
import sys
import json
import logging
import uuid
import redis

# Add the project root to Python path
sys.path.append(os.path.dirname(os.path.abspath(__file__)))

# Import after setting up path
from webserver.model import Task, Message, File
import webserver.datastore as ds
from webserver.logging_utils import setup_logging, log_service_startup, get_logger

# Setup logging with shared utility
setup_logging("redis-listener", log_level=logging.INFO)
logger = get_logger("redis-listener")

# Log startup information
log_service_startup("redis-listener")

def create_minimal_task_dict(task):
    """Create a minimal task dictionary to reduce payload size"""
    if not task:
        return None
    
    task_dict = task.to_dict()
    # Only include essential fields to reduce payload size
    minimal_task = {
        "task_id": task_dict.get("task_id"),
        "status": task_dict.get("status"),
        "title": task_dict.get("title"),
        "created_at": task_dict.get("created_at"),
        "user_id": task_dict.get("user_id")
    }
    return minimal_task

def create_minimal_message_dict(msg):
    """Create a minimal message dictionary to reduce payload size"""
    if not msg:
        return None
    
    msg_dict = msg.to_dict()
    # Only include essential fields to reduce payload size
    minimal_msg = {
        "message_id": msg_dict.get("message_id"),
        "task_id": msg_dict.get("task_id"),
        "role": msg_dict.get("role"),
        "content": msg_dict.get("content"),
        "created_at": msg_dict.get("created_at")
    }
    return minimal_msg

def redis_listener_standalone():
    """Standalone Redis listener that processes Celery events for database updates only"""
    listener_id = uuid.uuid4().hex[:8]
    
    # Connect to Redis
    r = redis.Redis(
        host=os.environ.get("REDIS_HOST", "localhost"),
        port=int(os.environ.get("REDIS_PORT", "6379"))
    )
    
    pubsub = r.pubsub()
    pubsub.subscribe("celery_updates")
    logger.info(f"[redis_listener] ({listener_id}) Started standalone Redis listener")
    
    try:
        for redis_message in pubsub.listen():
            logger.info(f"[redis_listener] ({listener_id}) Received message: {redis_message}")
            logger.debug(f"[redis_listener] ({listener_id}) Raw redis_message: {redis_message}")
            
            if redis_message["type"] != "message":
                continue
                
            try:
                raw_data = redis_message["data"]
                if isinstance(raw_data, bytes):
                    raw_data = raw_data.decode()
                event = json.loads(raw_data)
                logger.info(f"[redis_listener] ({listener_id}) Redis event: {event}")
                
                event_type = event.get("type")
                event_task_id = event.get("task_id")
                event_data = event.get("data")
                
                if event_task_id is None or event_data is None:
                    logger.warning(f"[redis_listener] ({listener_id}) Redis event missing required fields: event_type={event_type}, task_id={event_task_id}, data={event_data}")
                    continue
                    
                if event_type not in ("task_message", "task_file", "task_status_update"):
                    logger.debug(f"[redis_listener] ({listener_id}) Ignoring event_type={event_type}")
                    continue
                    
                task = Task.get_task(event_task_id)
                if not task:
                    logger.warning(f"[redis_listener] ({listener_id}) No task found for task_id={event_task_id}")
                    continue
                    
                if event_type == "task_message":
                    logger.info(f"[redis_listener] ({listener_id}) Processing task_message for task_id={event_task_id}")
                    msg = Message.process_event(task, event_data)
                    logger.info(f"[redis_listener] ({listener_id}) Processed message for task_id={event_task_id}")
                    
                elif event_type == "task_file":
                    logger.info(f"[redis_listener] ({listener_id}) Processing task_file for task_id={event_task_id}")
                    File.process_event(task, event_data)
                    logger.info(f"[redis_listener] ({listener_id}) Processed file event for task_id={event_task_id}")
                    
            except json.JSONDecodeError:
                logger.error(f"[redis_listener] ({listener_id}) Failed to decode Redis message as JSON.")
            except Exception as e:
                logger.error(f"[redis_listener] ({listener_id}) Redis listener error: {e}", exc_info=True)
                
    except Exception as e:
        logger.error(f"[redis_listener] ({listener_id}) Redis listener failed: {e}")
        raise

if __name__ == "__main__":
    logger.info("Starting standalone Redis listener...")
    redis_listener_standalone() 