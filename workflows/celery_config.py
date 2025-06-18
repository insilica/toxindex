import os
broker_url = os.environ.get('CELERY_BROKER_URL', 'redis://localhost:6379/0')
result_backend = os.environ.get('CELERY_RESULT_BACKEND', 'redis://localhost:6379/0')
broker_connection_retry_on_startup = True
worker_max_tasks_per_child = 1000
worker_hijack_root_logger = False
task_serializer = 'json'
result_serializer = 'json'
accept_content = ['json']
enable_utc = True