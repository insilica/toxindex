broker_url = 'redis://localhost:6379/0'
result_backend = 'redis://localhost:6379/0'
broker_connection_retry_on_startup = True
worker_max_tasks_per_child = 1000
task_serializer = 'json'
result_serializer = 'json'
accept_content = ['json']
enable_utc = True