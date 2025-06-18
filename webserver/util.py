def is_valid_uuid(val):
    import uuid
    try:
        uuid.UUID(str(val))
        return True
    except Exception:
        return False 