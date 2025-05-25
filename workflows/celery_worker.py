from celery import Celery


def make_celery(app_name=__name__):
    return Celery(
        app_name,
        broker='redis://localhost:6379/0',
        backend='redis://localhost:6379/0',
    )

celery = make_celery()
# celery -A workflows.celery_worker.celery worker --loglevel=info

# task modules
import workflows.probra as probra_task
import workflows.chat_response_task as chat_response_task