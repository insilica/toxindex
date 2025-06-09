import logging
import os
from datetime import datetime

# BD ENC
print("PGUSER:", os.getenv("PGUSER"))
print("PGHOST:", os.getenv("PGHOST"))
print("PGPORT:", os.getenv("PGPORT"))
print("PGDATABASE:", os.getenv("PGDATABASE"))
print("PG_SOCKET_DIR:", os.getenv("PG_SOCKET_DIR"))

# Ensure logs directory exists
os.makedirs(os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'logs')), exist_ok=True)
log_filename = os.path.abspath(os.path.join(
    os.path.dirname(__file__), '..', 'logs', f'app_{datetime.now().strftime("%Y-%m-%d_%H")}.log'))

logging.basicConfig(
    level=logging.DEBUG,
    format='%(asctime)s - %(name)s - %(levelname)s - %(filename)s:%(lineno)d - %(message)s',
    handlers=[
        logging.FileHandler(log_filename),
        logging.StreamHandler()
    ]
)

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
import workflows.plain_openai_tasks
