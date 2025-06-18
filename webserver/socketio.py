import os
from flask_socketio import SocketIO
message_queue_url = os.environ.get("SOCKETIO_MESSAGE_QUEUE", "redis://localhost:6379/0")
socketio = SocketIO(
    cors_allowed_origins="*",
    message_queue=message_queue_url,
    manage_session=False,
    logger=True,
    engineio_logger=True
) 