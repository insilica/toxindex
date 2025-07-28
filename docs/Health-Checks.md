# Health Check System Documentation

## 🏥 Overview

The ToxIndex application now includes a comprehensive health check system that monitors all critical components and services. This ensures high availability and early detection of issues in production.

## 📋 Health Check Endpoints

### 1. Basic Health Check (`/api/healthz`)
- **Purpose**: Simple connectivity test
- **Response**: `"ok"` with 200 status
- **Use Case**: Load balancer health checks, basic monitoring
- **Performance**: Very fast (< 1ms)

### 2. Liveness Probe (`/api/health/live`)
- **Purpose**: Kubernetes liveness probe
- **Response**: JSON with application status
- **Use Case**: Determines if pod should be restarted
- **Performance**: Very fast (< 1ms)

### 3. Readiness Probe (`/api/health/ready`)
- **Purpose**: Kubernetes readiness probe
- **Response**: JSON with critical service status
- **Use Case**: Determines if pod should receive traffic
- **Performance**: Fast (< 100ms)
- **Checks**: Database, Redis

### 4. Comprehensive Health Check (`/api/health`)
- **Purpose**: Full system health assessment
- **Response**: Detailed JSON with all component statuses
- **Use Case**: Monitoring dashboards, debugging
- **Performance**: Moderate (1-5 seconds)
- **Checks**: All components

## 🔍 Health Check Components

### Database Health Check
- **Tests**: Connection, table access, basic CRUD operations
- **Tables Checked**: `tasks`, `workflows`, `messages`, `files`, `chat_sessions`
- **Timeout**: 5 seconds
- **Status**: `healthy`, `error`

### Redis Health Check
- **Tests**: Connection, set/get operations, pub/sub
- **Timeout**: 5 seconds
- **Status**: `healthy`, `error`

### GCS Health Check
- **Tests**: Bucket access, storage connectivity
- **Timeout**: 10 seconds
- **Status**: `healthy`, `error`

### Celery Health Check
- **Tests**: Worker connectivity, task queue status
- **Timeout**: 5 seconds
- **Status**: `healthy`, `warning`, `error`

### External APIs Health Check
- **Tests**: OpenAI API connectivity (if configured)
- **Timeout**: 10 seconds
- **Status**: `healthy`, `warning`, `error`

### File Operations Health Check
- **Tests**: Log directory write permissions
- **Timeout**: 2 seconds
- **Status**: `healthy`, `error`

### Memory Usage Health Check
- **Tests**: Application memory consumption
- **Thresholds**: 
  - Warning: 80% of available memory
  - Error: 95% of available memory
- **Status**: `healthy`, `warning`, `error`

### Disk Space Health Check
- **Tests**: Available disk space
- **Thresholds**:
  - Warning: 80% used
  - Error: 95% used
- **Status**: `healthy`, `warning`, `error`

## 🚀 Kubernetes Integration

### Liveness Probe Configuration
```yaml
livenessProbe:
  httpGet:
    path: /api/health/live
    port: 6513
  periodSeconds: 30
  timeoutSeconds: 5
  failureThreshold: 3
```

### Readiness Probe Configuration
```yaml
readinessProbe:
  httpGet:
    path: /api/health/ready
    port: 6513
  periodSeconds: 10
  timeoutSeconds: 5
  failureThreshold: 3
  initialDelaySeconds: 30
```

## 📊 Response Format

### Basic Response
```json
{
  "status": "healthy",
  "message": "All critical services are ready",
  "timestamp": "2024-01-15T10:30:00.123456"
}
```

### Comprehensive Response
```json
{
  "status": "healthy",
  "http_status": 200,
  "checks": {
    "database": {
      "status": "healthy",
      "message": "Database is healthy",
      "response_time": 0.023,
      "tables_accessible": 5,
      "timestamp": "2024-01-15T10:30:00.123456"
    },
    "redis": {
      "status": "healthy",
      "message": "Redis is healthy",
      "response_time": 0.015,
      "timestamp": "2024-01-15T10:30:00.123456"
    },
    "gcs": {
      "status": "healthy",
      "message": "GCS is accessible",
      "response_time": 0.045,
      "bucket_name": "toxindex-uploads",
      "timestamp": "2024-01-15T10:30:00.123456"
    },
    "celery": {
      "status": "healthy",
      "message": "Celery workers detected: 2",
      "response_time": 0.012,
      "worker_count": 2,
      "timestamp": "2024-01-15T10:30:00.123456"
    },
    "external_apis": {
      "status": "healthy",
      "message": "External APIs: 1 checked",
      "response_time": 0.001,
      "checks": [
        {
          "service": "openai",
          "status": "healthy",
          "message": "OpenAI API key configured"
        }
      ],
      "timestamp": "2024-01-15T10:30:00.123456"
    },
    "file_operations": {
      "status": "healthy",
      "message": "File operations: healthy",
      "response_time": 0.003,
      "logs_directory": "/app/logs",
      "timestamp": "2024-01-15T10:30:00.123456"
    },
    "memory_usage": {
      "status": "healthy",
      "message": "Memory usage normal: 45.2%",
      "memory_mb": 256.5,
      "memory_percent": 45.2,
      "timestamp": "2024-01-15T10:30:00.123456"
    },
    "disk_space": {
      "status": "healthy",
      "message": "Disk space normal: 65.3% used",
      "disk_percent": 65.3,
      "disk_free_gb": 45.2,
      "disk_total_gb": 130.5,
      "timestamp": "2024-01-15T10:30:00.123456"
    }
  },
  "_metadata": {
    "total_checks": 8,
    "total_time": 0.104,
    "timestamp": "2024-01-15T10:30:00.123456"
  },
  "timestamp": "2024-01-15T10:30:00.123456"
}
```

## 🛡️ Safety Features

### Timeout Protection
- Each health check has a maximum timeout
- Prevents health checks from hanging
- Graceful degradation on timeouts

### Error Isolation
- Individual check failures don't affect others
- Detailed error messages for debugging
- Graceful handling of missing dependencies

### Performance Optimization
- Caching headers prevent unnecessary requests
- Fast fail on critical service failures
- Minimal resource usage

### Kubernetes Integration
- Proper HTTP status codes for K8s
- Separate liveness and readiness probes
- Configurable thresholds and timeouts

## 🔧 Configuration

### Environment Variables
```bash
# Health check timeouts (in seconds)
HEALTH_CHECK_TIMEOUT=30
HEALTH_CHECK_DB_TIMEOUT=5
HEALTH_CHECK_REDIS_TIMEOUT=5
HEALTH_CHECK_GCS_TIMEOUT=10

# Memory thresholds (percentages)
HEALTH_CHECK_MEMORY_WARNING=80
HEALTH_CHECK_MEMORY_ERROR=95

# Disk space thresholds (percentages)
HEALTH_CHECK_DISK_WARNING=80
HEALTH_CHECK_DISK_ERROR=95
```

### Logging
Health checks log their activities:
```python
# Log levels
logging.info("Running health check: database")
logging.info("Health check database: healthy")
logging.error("Health check database failed with exception: ...")
```

## 🧪 Testing

### Manual Testing
```bash
# Test basic health check
curl http://localhost:6513/api/healthz

# Test comprehensive health check
curl http://localhost:6513/api/health

# Test readiness probe
curl http://localhost:6513/api/health/ready

# Test liveness probe
curl http://localhost:6513/api/health/live
```

### Automated Testing
```bash
# Run the test script
python test_health_checks.py
```

### Kubernetes Testing
```bash
# Test pod health
kubectl get pods -n toxindex-app

# Check pod events
kubectl describe pod <pod-name> -n toxindex-app

# Check pod logs
kubectl logs <pod-name> -n toxindex-app
```

## 📈 Monitoring Integration

### Prometheus Metrics
Health checks can be extended to expose Prometheus metrics:
```python
# Example metric exposure
from prometheus_client import Counter, Histogram

health_check_total = Counter('health_check_total', 'Total health checks', ['check', 'status'])
health_check_duration = Histogram('health_check_duration_seconds', 'Health check duration', ['check'])
```

### Alerting Rules
```yaml
# Example Prometheus alerting rules
groups:
  - name: health_checks
    rules:
      - alert: HealthCheckFailed
        expr: health_check_total{status="error"} > 0
        for: 1m
        labels:
          severity: critical
        annotations:
          summary: "Health check failed"
          description: "Health check {{ $labels.check }} is failing"
```

## 🚨 Troubleshooting

### Common Issues

1. **Database Connection Failed**
   - Check database server status
   - Verify connection credentials
   - Check network connectivity

2. **Redis Connection Failed**
   - Check Redis server status
   - Verify Redis host/port configuration
   - Check network connectivity

3. **GCS Access Failed**
   - Verify GCS credentials
   - Check bucket permissions
   - Verify network connectivity

4. **Memory Usage High**
   - Check for memory leaks
   - Review application memory usage
   - Consider increasing pod memory limits

5. **Disk Space Low**
   - Clean up old logs
   - Check for large temporary files
   - Consider increasing disk size

### Debug Commands
```bash
# Check pod health
kubectl get pods -n toxindex-app -o wide

# View health check logs
kubectl logs <pod-name> -n toxindex-app | grep health

# Test health check manually
kubectl exec -it <pod-name> -n toxindex-app -- curl localhost:6513/api/health

# Check resource usage
kubectl top pods -n toxindex-app
```

## 🎯 Best Practices

1. **Use Appropriate Endpoints**
   - Use `/api/healthz` for load balancers
   - Use `/api/health/live` for liveness probes
   - Use `/api/health/ready` for readiness probes
   - Use `/api/health` for monitoring dashboards

2. **Set Proper Timeouts**
   - Liveness probe: 30s period, 5s timeout
   - Readiness probe: 10s period, 5s timeout
   - Comprehensive check: 60s period, 30s timeout

3. **Monitor Health Check Performance**
   - Track response times
   - Set up alerts for slow health checks
   - Monitor health check failure rates

4. **Regular Testing**
   - Test health checks after deployments
   - Verify health checks work in all environments
   - Document any environment-specific issues

## 🔄 Continuous Improvement

### Metrics to Track
- Health check response times
- Health check failure rates
- Component-specific health status
- Resource usage trends

### Regular Reviews
- Review health check thresholds monthly
- Update health checks when adding new components
- Optimize health check performance
- Add new health checks as needed

This comprehensive health check system ensures your application remains reliable and provides early warning of potential issues in production. 