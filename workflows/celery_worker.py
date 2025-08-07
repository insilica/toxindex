import logging
import os
from datetime import datetime
from webserver.data_paths import LOGS_ROOT

# BD ENC
print("PGUSER:", os.getenv("PGUSER"))
print("PGHOST:", os.getenv("PGHOST"))
print("PGPORT:", os.getenv("PGPORT"))
print("PGDATABASE:", os.getenv("PGDATABASE"))
print("PG_SOCKET_DIR:", os.getenv("PG_SOCKET_DIR"))

# Ensure logs directory exists
LOGS_ROOT().mkdir(parents=True, exist_ok=True)
log_filename = LOGS_ROOT() / f'app_{datetime.now().strftime("%Y-%m-%d_%H")}.log'

# Remove print statements and use env vars for Redis URLs
broker_url = os.environ.get("CELERY_BROKER_URL", "redis://localhost:6379/0")
result_backend = os.environ.get("CELERY_RESULT_BACKEND", "redis://localhost:6379/0")

logging.basicConfig(
    level=logging.INFO,  # Set to INFO for production
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
    broker=broker_url,
    backend=result_backend,
)
celery.config_from_object(celery_config)

# Import tasks so they are registered
import workflows.probra # noqa: F401
import workflows.plain_openai_tasks # noqa: F401
import workflows.raptool_task # noqa: F401
import workflows.pathway_analysis_task # noqa: F401
import workflows.celery_template_simple # noqa: F401
