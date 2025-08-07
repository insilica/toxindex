# Cloud Logging Setup

This document describes the Cloud Logging integration for toxindex components.

## Overview

The logging system provides consistent logging across all components with automatic Cloud Logging integration when running in GKE.

## Components

### 1. Shared Logging Utility (`webserver/logging_utils.py`)

Provides centralized logging configuration for all services:

- **Local Development**: File + stdout logging
- **GKE Production**: File + stdout + Cloud Logging
- **Automatic Detection**: Detects environment via `KUBERNETES_SERVICE_HOST`

### 2. Service-Specific Logging

Each component uses the shared utility:

- **Web Server**: `webserver/app.py`
- **Celery Worker**: `workflows/celery_worker.py`
- **Redis Listener**: `redis_listener_standalone.py`

## Usage

### Basic Setup

```python
from webserver.logging_utils import setup_logging, get_logger

# Setup logging for your service
setup_logging("my-service", log_level=logging.INFO)
logger = get_logger("my-service")

# Use the logger
logger.info("Service started")
logger.error("An error occurred")
```

### Service Startup/Shutdown

```python
from webserver.logging_utils import log_service_startup, log_service_shutdown

# Log startup with environment info
log_service_startup("my-service", custom_param="value")

# Log shutdown
log_service_shutdown("my-service")
```

## Log Destinations

### Local Development
- **File**: `data/logs/{service_name}_{YYYY-MM-DD_HH}.log`
- **Stdout**: Console output for container logs

### GKE Production
- **File**: `data/logs/{service_name}_{YYYY-MM-DD_HH}.log`
- **Stdout**: Container logs (captured by Kubernetes)
- **Cloud Logging**: Google Cloud Logging with service-specific log names

## Log Format

All logs use the same format:
```
%(asctime)s - %(name)s - %(levelname)s - %(filename)s:%(lineno)d - %(message)s
```

Example:
```
2025-08-05 17:30:45,123 - celery-worker - INFO - probra.py:45 - Task started
```

## Cloud Logging Configuration

### Automatic Setup
- Detects GKE environment via `KUBERNETES_SERVICE_HOST`
- Uses default Google Cloud credentials
- Logs are sent to Cloud Logging with service-specific log names

### Manual Configuration
If you need custom Cloud Logging setup:

```python
from google.cloud import logging as cloud_logging
from google.cloud.logging.handlers import CloudLoggingHandler

client = cloud_logging.Client()
cloud_handler = CloudLoggingHandler(client, name="custom-log-name")
```

## Testing

Run the test script to verify logging setup:

```bash
python test_cloud_logging.py
```

This will test:
- Local vs GKE environment detection
- All log levels (DEBUG, INFO, WARNING, ERROR)
- Structured logging
- Multiple service names

## Dependencies

Add to your `pyproject.toml`:
```toml
dependencies = [
    "google-cloud-logging>=3.0.0",
    # ... other dependencies
]
```

## Troubleshooting

### Cloud Logging Not Working
1. Check if `google-cloud-logging` is installed
2. Verify Google Cloud credentials are configured
3. Ensure running in GKE environment (`KUBERNETES_SERVICE_HOST` set)

### Local Logs Not Appearing
1. Check `data/logs/` directory exists
2. Verify write permissions
3. Check log level configuration

### Performance Issues
- Cloud Logging adds latency for each log message
- Consider using structured logging for better performance
- Use appropriate log levels (INFO for production, DEBUG for development)

## Best Practices

1. **Use Service-Specific Loggers**: Each component should have its own logger name
2. **Structured Logging**: Include relevant context in log messages
3. **Appropriate Log Levels**: Use DEBUG for development, INFO for production
4. **Error Handling**: Always handle logging setup failures gracefully
5. **Environment Detection**: Let the system automatically detect the environment

## Monitoring

### Local Development
- Check log files in `data/logs/`
- Monitor console output

### GKE Production
- Use `kubectl logs <pod-name>` for container logs
- Check Google Cloud Console > Logging for Cloud Logging
- Set up log-based alerts in Cloud Monitoring 