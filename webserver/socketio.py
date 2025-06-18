from flask_socketio import SocketIO
socketio = SocketIO(
    cors_allowed_origins="*",
    message_queue="redis://localhost:6379/0",
    manage_session=False,
    logger=True,
    engineio_logger=True
) 