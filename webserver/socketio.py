import os
from flask_socketio import SocketIO

redis_url = f"redis://{os.environ.get('REDIS_HOST', 'localhost')}:{os.environ.get('REDIS_PORT', '6379')}/0"
socketio = SocketIO(
    message_queue=redis_url,
    cors_allowed_origins=["https://www.toxindex.com"],
    logger=False,
    engineio_logger=False,
    manage_session=False,
    # Optimize polling settings
    ping_timeout=60,
    ping_interval=25,
    max_http_buffer_size=1e6,
    async_mode='gevent'
) 