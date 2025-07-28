# Data Paths Migration Guide

## 🚀 Overview

The ToxIndex application has migrated from local file storage to cloud-native services. This document outlines the changes and what remains local.

## 📊 Migration Summary

| **Data Type** | **Before** | **After** | **Reason** |
|---------------|------------|-----------|------------|
| **File Uploads** | Local filesystem | GCS bucket | Multi-pod access |
| **Chat Sessions** | Local files | Cloud SQL | Persistent, queryable |
| **Caching** | Local files | Redis Memorystore | Fast, shared |
| **Outputs** | Local files | GCS bucket | Multi-pod access |
| **Logs** | Local files | Local files | Debugging, monitoring |
| **Temp Files** | Local files | Local files | Processing only |

## 🔄 What Changed

### **1. File Storage Migration**
```python
# Before: Local file storage
uploads_dir = "/app/data/uploads"
outputs_dir = "/app/data/outputs"

# After: GCS storage
from webserver.storage import GCSFileStorage
gcs_storage = GCSFileStorage()
gcs_path = f"environments/{env_id}/files/{filename}"
```

### **2. Database Migration**
```python
# Before: Local SQLite/PostgreSQL
DATABASE_URL = "sqlite:///local.db"

# After: Cloud SQL
DATABASE_URL = "postgresql://user:pass@cloud-sql-host:5432/db"
```

### **3. Caching Migration**
```python
# Before: Local file cache
cache_dir = "/app/data/cache"
cache_file = cache_dir / "function_cache.json"

# After: Redis cache
import redis
r = redis.Redis(host=REDIS_HOST, port=REDIS_PORT)
r.set("cache_key", "cache_value")
```

### **4. Session Storage Migration**
```python
# Before: Local session files
session_dir = "/app/data/sessions"
session_file = session_dir / f"session_{user_id}.json"

# After: Cloud SQL database
session = ChatSession.get_session(session_id)
session.messages = messages
session.save()
```

## 📁 Current Local Paths

### **What Remains Local**

#### **1. Logs (`/data/logs`)**
```python
# Still local for debugging and monitoring
LOGS_ROOT().mkdir(exist_ok=True)
log_filename = LOGS_ROOT() / f"app_{datetime.now().strftime('%Y-%m-%d_%H')}.log"
```

**Why Local:**
- Easy access for debugging
- No network dependency for logging
- Real-time monitoring
- Kubernetes log aggregation

#### **2. Temporary Files (`/data/tmp`)**
```python
# Used for processing temporary files
import tempfile
with tempfile.NamedTemporaryFile(dir=TMP_ROOT()) as temp_file:
    # Process file
    pass
```

**Why Local:**
- Fast I/O during processing
- Automatic cleanup
- No network overhead for temporary data

#### **3. Local Cache (`/data/cache`)**
```python
# Only for file-based caches that can't use Redis
cache_path = get_cache_path("chemprop_predictions")
```

**Why Local:**
- Large file caches that don't fit in Redis
- Binary data that's expensive to serialize
- Legacy cache compatibility

## ☁️ Cloud Service Configuration

### **1. Redis Memorystore**
```python
# Configuration
REDIS_HOST = os.environ.get("REDIS_HOST", "localhost")
REDIS_PORT = int(os.environ.get("REDIS_PORT", "6379"))

# Usage
import redis
r = redis.Redis(host=REDIS_HOST, port=REDIS_PORT)
```

**Purpose:**
- Application caching
- Session storage
- Task queue (Celery)
- Real-time messaging

### **2. Google Cloud Storage**
```python
# Configuration
GCS_BUCKET_NAME = os.environ.get("GCS_BUCKET_NAME", "toxindex-uploads")

# Usage
from webserver.storage import GCSFileStorage
gcs_storage = GCSFileStorage()
```

**Purpose:**
- File uploads
- Analysis outputs
- Generated reports
- Multi-pod file sharing

### **3. Cloud SQL**
```python
# Configuration
DATABASE_URL = f"postgresql://{user}:{pass}@{host}:{port}/{db}"

# Usage
from webserver.model import Task, Workflow, Message
task = Task.get_task(task_id)
```

**Purpose:**
- User data
- Chat sessions
- Task history
- Application state

## 🔧 Updated Configuration

### **Environment Variables**
```bash
# Cloud SQL
PGHOST=your-cloud-sql-host
PGPORT=5432
PGDATABASE=toxindex
PGUSER=toxindex_user
PGPASSWORD=your-password

# Redis Memorystore
REDIS_HOST=your-redis-host
REDIS_PORT=6379

# GCS
GCS_BUCKET_NAME=toxindex-uploads
```

### **Kubernetes Secrets**
```yaml
# k8s/backend-secrets.yaml
apiVersion: v1
kind: Secret
metadata:
  name: backend-secrets
  namespace: toxindex-app
type: Opaque
data:
  PGHOST: <base64-encoded-host>
  PGPASSWORD: <base64-encoded-password>
  GCS_BUCKET_NAME: <base64-encoded-bucket-name>
```

## 📈 Benefits of Migration

### **1. Scalability**
- **Multi-pod access**: Files accessible from any pod
- **Horizontal scaling**: Add more pods without data conflicts
- **Load distribution**: Share load across multiple instances

### **2. Reliability**
- **Persistent storage**: Data survives pod restarts
- **Backup and recovery**: Cloud-managed backups
- **High availability**: Cloud service SLAs

### **3. Performance**
- **Fast caching**: Redis for high-speed data access
- **CDN integration**: GCS with Cloud CDN
- **Optimized queries**: Cloud SQL performance tuning

### **4. Cost Efficiency**
- **Pay-per-use**: Only pay for what you use
- **No storage management**: Cloud handles infrastructure
- **Automatic scaling**: Scale based on demand

## 🛠️ Code Changes Required

### **1. File Operations**
```python
# Before
def upload_file(file, filename):
    file_path = os.path.join(uploads_dir, filename)
    file.save(file_path)
    return file_path

# After
def upload_file(file, filename):
    gcs_storage = GCSFileStorage()
    gcs_path = f"uploads/{uuid.uuid4()}_{filename}"
    gcs_storage.upload_file(file, gcs_path)
    return gcs_path
```

### **2. Caching**
```python
# Before
def get_cached_result(key):
    cache_file = cache_dir / f"{key}.json"
    if cache_file.exists():
        return json.loads(cache_file.read_text())

# After
def get_cached_result(key):
    cached = r.get(f"cache:{key}")
    if cached:
        return json.loads(cached)
```

### **3. Session Management**
```python
# Before
def save_session(session_id, data):
    session_file = sessions_dir / f"{session_id}.json"
    session_file.write_text(json.dumps(data))

# After
def save_session(session_id, data):
    session = ChatSession.get_or_create(session_id)
    session.data = data
    session.save()
```

## 🔍 Monitoring and Debugging

### **1. Cloud Service Health**
```python
# Health checks for cloud services
def check_redis_health():
    r = redis.Redis(host=REDIS_HOST, port=REDIS_PORT)
    return r.ping()

def check_gcs_health():
    gcs_storage = GCSFileStorage()
    return gcs_storage.bucket_exists()

def check_sql_health():
    # Test database connection
    pass
```

### **2. Local Path Monitoring**
```python
# Monitor local disk usage
def check_local_storage():
    logs_usage = LOGS_ROOT().stat().st_size
    tmp_usage = TMP_ROOT().stat().st_size
    return {"logs_mb": logs_usage / 1024 / 1024, "tmp_mb": tmp_usage / 1024 / 1024}
```

## 🚨 Troubleshooting

### **Common Issues**

#### **1. Redis Connection Failed**
```bash
# Check Redis connectivity
kubectl exec -it <pod-name> -n toxindex-app -- redis-cli -h $REDIS_HOST ping

# Check Redis logs
kubectl logs <redis-pod> -n redis
```

#### **2. GCS Access Denied**
```bash
# Check GCS permissions
gcloud storage ls gs://toxindex-uploads/

# Check service account
kubectl describe secret backend-secrets -n toxindex-app
```

#### **3. Cloud SQL Connection Failed**
```bash
# Check database connectivity
kubectl exec -it <pod-name> -n toxindex-app -- psql $DATABASE_URL

# Check Cloud SQL logs
gcloud sql logs tail --instance=toxindex-sql
```

### **Debug Commands**
```bash
# Check all cloud services
python -c "from webserver.data_paths import get_cloud_services_summary; print(get_cloud_services_summary())"

# Check local paths
python -c "from webserver.data_paths import LOGS_ROOT, TMP_ROOT; print(f'Logs: {LOGS_ROOT()}, Temp: {TMP_ROOT()}')"

# Test Redis
python -c "import redis; r = redis.Redis(); print(r.ping())"
```

## 📋 Migration Checklist

### **✅ Completed**
- [x] File uploads migrated to GCS
- [x] Chat sessions migrated to Cloud SQL
- [x] Caching migrated to Redis
- [x] Outputs migrated to GCS
- [x] Health checks updated
- [x] Environment variables configured
- [x] Kubernetes secrets updated

### **🔄 In Progress**
- [ ] Monitor cloud service performance
- [ ] Optimize cache strategies
- [ ] Set up cloud service alerts
- [ ] Document cloud service costs

### **📝 Future Improvements**
- [ ] Implement GCS lifecycle policies
- [ ] Set up Redis persistence
- [ ] Configure Cloud SQL read replicas
- [ ] Implement data archival strategies

## 🎯 Best Practices

### **1. Local vs Cloud Decision**
- **Local**: Logs, temp files, processing data
- **Cloud**: Persistent data, shared data, user data

### **2. Caching Strategy**
- **Redis**: Small, frequently accessed data
- **Local**: Large files, binary data, legacy caches

### **3. Error Handling**
- **Graceful degradation**: Fall back to local storage if cloud services fail
- **Retry logic**: Implement exponential backoff for cloud operations
- **Monitoring**: Track cloud service performance and costs

### **4. Security**
- **Environment variables**: Never hardcode credentials
- **IAM roles**: Use least privilege access
- **Encryption**: Enable encryption at rest and in transit

This migration provides a robust, scalable foundation for your cloud-native application while maintaining the simplicity of local storage where appropriate. 