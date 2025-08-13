# Cloud Logging Setup

This document describes the Cloud Logging integration for toxindex components.

## Overview

The logging system provides consistent logging across all components with automatic Cloud Logging integration when running in GKE. The system includes a resilient Cloud Logging handler that gracefully handles timeouts and connection issues.

## Components

### 1. Shared Logging Utility (`webserver/logging_utils.py`)

Provides centralized logging configuration for all services:

- **Local Development**: File + stdout logging
- **GKE Production**: File + stdout + Cloud Logging (with timeout resilience)
- **Automatic Detection**: Detects environment via `KUBERNETES_SERVICE_HOST`
- **Graceful Fallback**: Falls back to local logging if Cloud Logging fails

### 2. Resilient Cloud Logging Handler

The `ResilientCloudLoggingHandler` class provides:
- **Timeout Protection**: Configures Cloud Logging with shorter timeouts
- **Error Handling**: Gracefully handles Cloud Logging failures
- **Fallback Logging**: Automatically falls back to local logging
- **Batch Optimization**: Uses smaller batch sizes to prevent timeouts

### 3. Service-Specific Logging

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
- **Fallback**: Local logging if Cloud Logging fails

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
- Includes timeout and error handling

### Manual Configuration
If you need custom Cloud Logging setup:

```python
from google.cloud import logging as cloud_logging
from google.cloud.logging.handlers import CloudLoggingHandler

client = cloud_logging.Client()
cloud_handler = CloudLoggingHandler(client, name="custom-log-name")
```

### Environment Variables

#### `DISABLE_CLOUD_LOGGING`
Set to `"true"` to disable Cloud Logging entirely:
```bash
export DISABLE_CLOUD_LOGGING=true
```

This is useful for:
- Debugging Cloud Logging issues
- Reducing latency in development
- Emergency situations where Cloud Logging is causing problems

## Testing

Run the test script to verify logging setup:

```bash
python test_logging_fix.py
```

This will test:
- Local vs GKE environment detection
- All log levels (DEBUG, INFO, WARNING, ERROR)
- Structured logging
- Multiple service names
- Cloud Logging timeout resilience
- Fallback logging functionality

## Dependencies

Add to your `pyproject.toml`:
```toml
dependencies = [
    "google-cloud-logging>=3.0.0",
    # ... other dependencies
]
```

## Troubleshooting

### Cloud Logging Timeout Issues
**Problem**: `504 Deadline Exceeded` errors in Cloud Logging
**Solution**: The resilient handler automatically handles this by:
1. Using smaller batch sizes (5 instead of default)
2. Shorter max latency (15 seconds instead of default)
3. Falling back to local logging on failures
4. Using HTTP instead of gRPC to avoid gRPC timeout issues

### Cloud Logging Not Working
1. Check if `google-cloud-logging` is installed
2. Verify Google Cloud credentials are configured
3. Ensure running in GKE environment (`KUBERNETES_SERVICE_HOST` set)
4. Check if `DISABLE_CLOUD_LOGGING` is set to `"true"`

### Local Logs Not Appearing
1. Check `data/logs/` directory exists
2. Verify write permissions
3. Check log level configuration

### Performance Issues
- Cloud Logging adds latency for each log message
- Consider using structured logging for better performance
- Use appropriate log levels (INFO for production, DEBUG for development)
- Set `DISABLE_CLOUD_LOGGING=true` if Cloud Logging is causing performance issues

## Best Practices

1. **Use Service-Specific Loggers**: Each component should have its own logger name
2. **Structured Logging**: Include relevant context in log messages
3. **Appropriate Log Levels**: Use DEBUG for development, INFO for production
4. **Error Handling**: Always handle logging setup failures gracefully
5. **Environment Detection**: Let the system automatically detect the environment
6. **Timeout Resilience**: The system automatically handles Cloud Logging timeouts
7. **Fallback Strategy**: Always have local logging as a fallback

## Monitoring

### Local Development
- Check log files in `data/logs/`
- Monitor console output

### GKE Production
- Use `kubectl logs <pod-name>` for container logs
- Check Google Cloud Console > Logging for Cloud Logging
- Monitor for timeout errors in Cloud Logging
- Check local log files if Cloud Logging is failing

## Emergency Procedures

### If Cloud Logging is Causing Issues
1. **Quick Fix**: Set `DISABLE_CLOUD_LOGGING=true` in the deployment
2. **Redeploy**: Update the Kubernetes deployment
3. **Monitor**: Check that local logging is working
4. **Investigate**: Debug Cloud Logging issues separately

### Deployment Configuration
```yaml
env:
  - name: DISABLE_CLOUD_LOGGING
    value: "true"  # Temporarily disable Cloud Logging
```

## Recent Improvements

### v2.0 - Resilient Cloud Logging
- Added `ResilientCloudLoggingHandler` class
- Implemented timeout protection
- Added graceful fallback to local logging
- Reduced batch sizes to prevent timeouts
- Added environment variable control
- Improved error handling and recovery 