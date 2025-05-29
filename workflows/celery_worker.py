from celery import Celery
from . import celery_config

celery = Celery(
    'workflows',
    broker='redis://localhost:6379/0',
    backend='redis://localhost:6379/0',
)
celery.config_from_object(celery_config)

# Import tasks so they are registered
import workflows.probra
import workflows.chat_response_task
import workflows.interactive_echo_task
