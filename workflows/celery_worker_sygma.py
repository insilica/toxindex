import os
import logging
from webserver.logging_utils import setup_logging, log_service_startup, get_logger

# Remove print statements and use env vars for Redis URLs
broker_url = os.environ.get("CELERY_BROKER_URL", "redis://localhost:6379/0")
result_backend = os.environ.get("CELERY_RESULT_BACKEND", "redis://localhost:6379/0")

from celery import Celery
from . import celery_config

celery = Celery(
    'workflows',
    broker=broker_url,
    backend=result_backend,
)
celery.config_from_object(celery_config)

# Configure Celery to use our logging setup
celery.conf.update(
    worker_hijack_root_logger=False,
    worker_log_format='[%(asctime)s: %(levelname)s/%(processName)s] %(message)s',
    worker_task_log_format='[%(asctime)s: %(levelname)s/%(processName)s][%(task_name)s(%(task_id)s)] %(message)s',
    worker_log_color=True,
    task_track_started=True,
    task_send_sent_event=True,
    worker_send_task_events=True,
)

# Import ONLY the metabolite-sygma task
import workflows.sygma_docker.metabolite_sygma_task # noqa: F401

def setup_celery_worker():
    """Setup logging and startup for celery worker - only call this when actually starting a worker"""
    # Setup logging with shared utility
    setup_logging("celery-worker-sygma", log_level=logging.INFO)
    logger = get_logger("celery-worker-sygma")

    # Log startup information
    log_service_startup("celery-worker-sygma")

    # Log registered tasks
    logger.info(f"Registered tasks: {list(celery.tasks.keys())}")

# Only setup logging if this module is run directly (i.e., as a celery worker)
if __name__ == '__main__':
    setup_celery_worker()