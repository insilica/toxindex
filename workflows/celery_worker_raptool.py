import os
import logging
from webserver.logging_utils import setup_logging, log_service_startup, get_logger

# Remove print statements and use env vars for Redis URLs
broker_url = os.environ.get("CELERY_BROKER_URL", "redis://localhost:6379/0")
result_backend = os.environ.get("CELERY_RESULT_BACKEND", "redis://localhost:6379/0")

from workflows.celery_app import celery

# Import tasks so they are registered
import workflows.raptool_task # noqa: F401

def setup_celery_worker():
    """Setup logging and startup for celery worker - only call this when actually starting a worker"""
    worker_name = "celery-worker-raptool_task"

    # Setup logging with shared utility
    setup_logging(worker_name, log_level=logging.INFO)
    logger = get_logger(worker_name)

    # Log startup information
    log_service_startup(worker_name)
    
    # Log registered tasks
    logger.info(f"Registered tasks: {list(celery.tasks.keys())}")

# Only setup logging if this module is run directly (i.e., as a celery worker)
if __name__ == '__main__':
    setup_celery_worker()

