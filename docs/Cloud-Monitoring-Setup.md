# Google Cloud Monitoring Setup Guide

## 🎯 Overview

This guide sets up Google Cloud Monitoring for your ToxIndex application with deployment verification to ensure **all functions are working** before marking a deployment as successful.

## 🚀 **Why This Approach?**

### **✅ Deployment Verification Strategy**
```yaml
# Traditional Deployment
1. Deploy code
2. Check if pods are running
3. Mark as successful

# Our Enhanced Deployment
1. Deploy code
2. Wait for pods to be ready
3. Run comprehensive health checks
4. Verify Google Cloud Monitoring metrics
5. Test all critical functions
6. Only then mark as successful
```

### **✅ Benefits**
- **Zero broken deployments**: Only successful if all functions work
- **Real-time monitoring**: Cloud Monitoring integration
- **Comprehensive verification**: Tests all components
- **Automatic rollback**: Failed verification triggers rollback

## 📋 **Setup Steps**

### **1. Enable Google Cloud Monitoring**

```bash
# Enable required APIs
gcloud services enable monitoring.googleapis.com
gcloud services enable cloudresourcemanager.googleapis.com

# Verify monitoring is enabled
gcloud monitoring metrics list --limit=5
```

### **2. Set Up Service Account Permissions**

```bash
# Get your project number
PROJECT_NUMBER=$(gcloud projects describe $(gcloud config get-value project) --format='value(projectNumber)')

# Create service account for monitoring
gcloud iam service-accounts create toxindex-monitoring \
  --display-name="ToxIndex Monitoring Service Account"

# Grant monitoring permissions
gcloud projects add-iam-policy-binding $(gcloud config get-value project) \
  --member="serviceAccount:toxindex-monitoring@$(gcloud config get-value project).iam.gserviceaccount.com" \
  --role="roles/monitoring.metricWriter"

gcloud projects add-iam-policy-binding $(gcloud config get-value project) \
  --member="serviceAccount:toxindex-monitoring@$(gcloud config get-value project).iam.gserviceaccount.com" \
  --role="roles/monitoring.viewer"
```

### **3. Update Kubernetes Secrets**

```yaml
# Add to k8s/backend-secrets.yaml
apiVersion: v1
kind: Secret
metadata:
  name: backend-secrets
  namespace: toxindex-app
type: Opaque
data:
  # ... existing secrets ...
  GOOGLE_CLOUD_PROJECT: <base64-encoded-project-id>
```

### **4. Create Monitoring Dashboard**

```bash
# Create the dashboard
gcloud monitoring dashboards create --config-from-file=monitoring/dashboard.yaml

# Verify dashboard creation
gcloud monitoring dashboards list
```

### **5. Set Up Cloud Build Trigger**

```bash
# Create trigger for monitoring-enabled deployment
gcloud builds triggers create github \
  --repo-name=toxindex \
  --repo-owner=your-github-username \
  --branch-pattern="^gke$" \
  --build-config=cloudbuild-with-monitoring.yaml \
  --name="deploy-with-monitoring"
```

## 🔧 **Configuration Details**

### **1. Cloud Build with Monitoring (`cloudbuild-with-monitoring.yaml`)**

#### **Key Features:**
```yaml
# Deployment phases:
1. Backend deployment (Docker build → GKE)
2. Frontend deployment (Build → GCS)
3. Verification phase (Pod readiness)
4. Functional verification (Health checks)
5. Cloud Monitoring verification (Metrics collection)
6. Critical function testing (Database, Redis, GCS, Celery)
```

#### **Verification Steps:**
```bash
# 1. Basic health check
curl /api/healthz

# 2. Comprehensive health check
curl /api/health

# 3. Frontend accessibility
curl https://your-frontend-domain.com

# 4. Cloud Monitoring metrics
gcloud monitoring metrics list

# 5. Function verification
- Database connectivity
- Redis connectivity  
- GCS connectivity
- Celery worker status
```

### **2. Custom Metrics Exporter (`webserver/metrics_exporter.py`)**

#### **Dependencies:**
```toml
# Added to pyproject.toml
psutil>=5.9.0
google-cloud-monitoring>=2.0.0
```

#### **Metrics Exported:**
```python
# Health check metrics
- health_check_status (overall and per component)
- database_connections
- redis_operations
- celery_worker_count
- memory_usage_percent
- disk_usage_percent

# HTTP metrics
- http_requests_total
- http_request_duration_seconds
- http_errors_total

# Custom business metrics
- business_* (any custom metrics)
```

### **3. Monitoring Dashboard (`monitoring/dashboard.yaml`)**

#### **Dashboard Widgets:**
```yaml
# Infrastructure metrics
- Pod health status
- CPU usage by pod
- Memory usage by pod
- Network traffic

# Application metrics
- Application health checks
- HTTP error rate
- Response time
- Celery worker status

# Custom metrics
- Database connections
- Redis memory usage
- GCS operations
- Application metrics summary
```

## 📊 **Monitoring Dashboard**

### **1. Access the Dashboard**

```bash
# Open in browser
gcloud monitoring dashboards list --format="table(name,displayName)"

# Or access via Google Cloud Console
# https://console.cloud.google.com/monitoring/dashboards
```

### **2. Dashboard Sections**

#### **Infrastructure Health:**
- Pod status and resource usage
- Container metrics
- Network traffic patterns

#### **Application Health:**
- Health check status for all components
- Response times and error rates
- Custom application metrics

#### **Business Metrics:**
- User activity
- Feature usage
- Performance indicators

## 🚨 **Alerting Configuration**

### **1. Create Alerting Policies**

```bash
# High error rate alert
gcloud alpha monitoring policies create \
  --policy-from-file=monitoring/alerts/high-error-rate.yaml

# Health check failure alert
gcloud alpha monitoring policies create \
  --policy-from-file=monitoring/alerts/health-check-failure.yaml

# Resource usage alert
gcloud alpha monitoring policies create \
  --policy-from-file=monitoring/alerts/high-resource-usage.yaml
```

### **2. Alert Configuration Files**

#### **High Error Rate Alert:**
```yaml
# monitoring/alerts/high-error-rate.yaml
displayName: "High Error Rate Alert"
conditions:
  - displayName: "Error rate > 5%"
    conditionThreshold:
      filter: 'metric.type="custom.googleapis.com/http_errors_total"'
      comparison: COMPARISON_GREATER_THAN
      thresholdValue: 0.05
      duration: 300s
notificationChannels:
  - projects/your-project/notificationChannels/your-channel
```

#### **Health Check Failure Alert:**
```yaml
# monitoring/alerts/health-check-failure.yaml
displayName: "Health Check Failure Alert"
conditions:
  - displayName: "Health check status = 0"
    conditionThreshold:
      filter: 'metric.type="custom.googleapis.com/health_check_status"'
      comparison: COMPARISON_LESS_THAN
      thresholdValue: 1.0
      duration: 60s
notificationChannels:
  - projects/your-project/notificationChannels/your-channel
```

## 🔍 **Deployment Verification Process**

### **1. Pre-Deployment Checks**
```bash
# Verify monitoring is accessible
gcloud monitoring metrics list --limit=1

# Check service account permissions
gcloud auth list
```

### **2. Deployment Process**
```bash
# Trigger deployment with monitoring
git push origin gke  # Triggers Cloud Build

# Monitor deployment progress
gcloud builds list --limit=5
```

### **3. Post-Deployment Verification**
```bash
# Check deployment status
kubectl get pods -n toxindex-app

# Verify health checks
curl https://your-app.com/api/health

# Check monitoring metrics
gcloud monitoring metrics list --filter="metric.type:custom.googleapis.com"
```

## 📈 **Performance Monitoring**

### **1. Key Metrics to Track**

#### **Infrastructure Metrics:**
```yaml
- CPU usage per pod
- Memory usage per pod
- Network traffic
- Pod restart count
```

#### **Application Metrics:**
```yaml
- Health check status
- Response times
- Error rates
- Custom business metrics
```

#### **Business Metrics:**
```yaml
- Active users
- Feature usage
- Performance indicators
- User satisfaction scores
```

### **2. SLO/SLI Definition**

```yaml
# Service Level Objectives
availability: 99.9%
response_time: < 1 second (95th percentile)
error_rate: < 1%

# Service Level Indicators
- Health check success rate
- HTTP response time
- Error rate by endpoint
- Database connection success rate
```

## 🛠️ **Troubleshooting**

### **1. Common Issues**

#### **Monitoring API Access:**
```bash
# Check permissions
gcloud auth list
gcloud projects get-iam-policy $(gcloud config get-value project)

# Test monitoring access
gcloud monitoring metrics list --limit=1
```

#### **Custom Metrics Not Appearing:**
```bash
# Check metric descriptors
gcloud monitoring metrics list --filter="metric.type:custom.googleapis.com"

# Verify metric creation
gcloud monitoring metrics list --filter="metric.type:custom.googleapis.com/health_check_status"
```

#### **Dashboard Not Loading:**
```bash
# Check dashboard configuration
gcloud monitoring dashboards list

# Verify dashboard permissions
gcloud monitoring dashboards describe <dashboard-name>
```

### **2. Debug Commands**

```bash
# Check monitoring setup
gcloud services list --enabled --filter="name:monitoring"

# Test metric writing
python -c "
from webserver.metrics_exporter import export_custom_metric
export_custom_metric('test_metric', 1.0)
"

# Verify dependencies are installed
pip list | grep -E "(psutil|google-cloud-monitoring)"

# Verify metrics in monitoring
gcloud monitoring metrics list --filter="metric.type:custom.googleapis.com/test_metric"
```

## 🎯 **Best Practices**

### **1. Metric Design**
- Use descriptive metric names
- Include appropriate labels
- Set reasonable thresholds
- Monitor metric cardinality

### **2. Alerting Strategy**
- Set up alerts for critical issues
- Use different severity levels
- Include context in alert messages
- Test alerting regularly

### **3. Dashboard Design**
- Group related metrics
- Use appropriate visualizations
- Include thresholds and targets
- Keep dashboards focused

### **4. Deployment Verification**
- Test all critical functions
- Verify monitoring integration
- Check alerting setup
- Document any issues

## 🚀 **Next Steps**

### **1. Immediate Actions**
```bash
# 1. Enable monitoring APIs
gcloud services enable monitoring.googleapis.com

# 2. Create dashboard
gcloud monitoring dashboards create --config-from-file=monitoring/dashboard.yaml

# 3. Set up alerts
gcloud alpha monitoring policies create --policy-from-file=monitoring/alerts/health-check-failure.yaml

# 4. Test deployment
git push origin gke
```

### **2. Ongoing Monitoring**
- Review dashboard regularly
- Adjust alert thresholds
- Add new metrics as needed
- Optimize performance

### **3. Advanced Features**
- Set up custom dashboards
- Implement advanced alerting
- Add business metrics
- Create performance reports

This setup ensures your deployments are **only marked successful when all functions are working** and provides comprehensive monitoring for ongoing operations. 