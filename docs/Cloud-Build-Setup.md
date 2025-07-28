# Cloud Build CI/CD Setup Guide

This guide explains how to set up automated deployment pipelines using Google Cloud Build for the ToxIndex application.

## 🚀 Overview

We have three main Cloud Build configurations:

1. **Backend Deploy** (`cloudbuild.yaml`) - Deploys backend to GKE
2. **Frontend Deploy** (`cloudbuild-frontend.yaml`) - Deploys frontend to GCS
3. **Full Stack Deploy** (`cloudbuild-fullstack.yaml`) - Coordinates both deployments

## 🔧 Setup Instructions

### 1. Enable Required APIs

```bash
# Enable Cloud Build API
gcloud services enable cloudbuild.googleapis.com

# Enable Container Registry API
gcloud services enable containerregistry.googleapis.com

# Enable Artifact Registry API
gcloud services enable artifactregistry.googleapis.com

# Enable Kubernetes Engine API
gcloud services enable container.googleapis.com
```

### 2. Set Up Cloud Build Service Account

Cloud Build automatically uses a service account. Grant it the necessary permissions:

```bash
# Get the Cloud Build service account
PROJECT_NUMBER=$(gcloud projects describe toxindex --format='value(projectNumber)')
CLOUDBUILD_SA="${PROJECT_NUMBER}@cloudbuild.gserviceaccount.com"

# Grant necessary roles
gcloud projects add-iam-policy-binding toxindex \
  --member="serviceAccount:${CLOUDBUILD_SA}" \
  --role="roles/container.developer"

gcloud projects add-iam-policy-binding toxindex \
  --member="serviceAccount:${CLOUDBUILD_SA}" \
  --role="roles/storage.admin"

gcloud projects add-iam-policy-binding toxindex \
  --member="serviceAccount:${CLOUDBUILD_SA}" \
  --role="roles/artifactregistry.writer"

gcloud projects add-iam-policy-binding toxindex \
  --member="serviceAccount:${CLOUDBUILD_SA}" \
  --role="roles/compute.admin"
```

### 3. Create Cloud Build Triggers

```bash
# Import triggers from configuration
gcloud builds triggers import --source=cloudbuild-triggers.yaml
```

Or create them manually:

```bash
# Backend trigger
gcloud builds triggers create github \
  --repo-name=toxindex \
  --repo-owner=your-github-username \
  --branch-pattern="^main$" \
  --build-config=cloudbuild.yaml \
  --included-files="webserver/**,workflows/**,requirements.txt,pyproject.toml,Dockerfile" \
  --name="backend-deploy"

# Frontend trigger
gcloud builds triggers create github \
  --repo-name=toxindex \
  --repo-owner=your-github-username \
  --branch-pattern="^main$" \
  --build-config=cloudbuild-frontend.yaml \
  --included-files="frontend/**,package.json,package-lock.json" \
  --name="frontend-deploy"

# Full stack trigger
gcloud builds triggers create github \
  --repo-name=toxindex \
  --repo-owner=your-github-username \
  --branch-pattern="^main$" \
  --build-config=cloudbuild-fullstack.yaml \
  --name="full-stack-deploy"
```

## 📋 Configuration Details

### Backend Deployment (`cloudbuild.yaml`)

**Process:**
1. Build Docker image
2. Push to Artifact Registry
3. Get GKE credentials
4. Update Kubernetes deployment
5. Restart backend and Celery worker
6. Verify deployment

**Manual Commands Replaced:**
```bash
# Before (manual)
docker system prune -af
docker volume prune -f
docker build -t us-docker.pkg.dev/toxindex/toxindex-backend/backend:latest .
docker push us-docker.pkg.dev/toxindex/toxindex-backend/backend:latest
kubectl rollout restart deployment backend -n toxindex-app
kubectl rollout restart deployment celery-worker -n toxindex-app

# After (automated)
git push origin main  # Triggers Cloud Build automatically
```

### Frontend Deployment (`cloudbuild-frontend.yaml`)

**Process:**
1. Install Node.js
2. Install dependencies
3. Build React app
4. Sync to GCS bucket
5. Set cache headers
6. Invalidate CDN cache
7. Verify deployment

**Manual Commands Replaced:**
```bash
# Before (manual)
npm run build
gsutil -m rsync -r dist/ gs://toxindex-react
gsutil -m setmeta -h "Cache-Control:no-store" gs://toxindex-react/**
gcloud compute url-maps invalidate-cdn-cache frontend-url-map --path "/*"

# After (automated)
git push origin main  # Triggers Cloud Build automatically
```

### Full Stack Deployment (`cloudbuild-fullstack.yaml`)

**Process:**
1. Deploy backend (Docker build → GKE)
2. Deploy frontend (Build → GCS)
3. Verify both deployments
4. Show deployment summary

## 🎯 Benefits

### ✅ Native Google Cloud Integration
- No external CI/CD service needed
- Native integration with GKE, GCS, Artifact Registry
- Automatic authentication and permissions

### ✅ Cost Effective
- Pay only for build time
- No GitHub Actions minutes
- Efficient resource usage

### ✅ Scalable
- Automatic scaling based on demand
- Parallel builds supported
- High-performance build machines

### ✅ Secure
- Service account-based authentication
- No external secrets needed
- Built-in security features

## 🔍 Monitoring

### Cloud Build Console
- View all builds and their status
- Check build logs in real-time
- Monitor build performance

### Cloud Build CLI
```bash
# List recent builds
gcloud builds list

# View build details
gcloud builds describe BUILD_ID

# View build logs
gcloud builds log BUILD_ID

# List triggers
gcloud builds triggers list
```

### Kubernetes Monitoring
```bash
# Check deployment status
kubectl get pods -n toxindex-app

# View deployment logs
kubectl logs deployment/backend -n toxindex-app

# Check rollout status
kubectl rollout status deployment/backend -n toxindex-app
```

## 🚨 Troubleshooting

### Common Issues

1. **Authentication Errors**
   ```bash
   # Check service account permissions
   gcloud projects get-iam-policy toxindex
   
   # Verify Cloud Build service account
   gcloud projects describe toxindex --format='value(projectNumber)'
   ```

2. **Build Failures**
   ```bash
   # Check build logs
   gcloud builds log BUILD_ID
   
   # Test locally
   docker build -t test-image .
   ```

3. **Deployment Failures**
   ```bash
   # Check GKE cluster status
   gcloud container clusters describe toxindex-cluster --zone=us-central1-a
   
   # Check pod status
   kubectl get pods -n toxindex-app
   kubectl describe pod POD_NAME -n toxindex-app
   ```

4. **GCS Sync Issues**
   ```bash
   # Check bucket permissions
   gsutil iam get gs://toxindex-react
   
   # Test GCS access
   gsutil ls gs://toxindex-react/
   ```

### Debug Commands

```bash
# Check trigger status
gcloud builds triggers list

# Manually trigger a build
gcloud builds triggers run TRIGGER_NAME --branch=main

# Check build history
gcloud builds list --limit=10

# View build configuration
gcloud builds describe BUILD_ID --format='value(steps)'
```

## 🔄 Rollback Process

If a deployment fails:

1. **Automatic**: Cloud Build includes health checks
2. **Manual**: Use previous Docker image tag
3. **Emergency**: Revert to last known good commit

```bash
# Manual rollback
kubectl set image deployment/backend backend=us-docker.pkg.dev/toxindex/toxindex-backend/backend:previous-tag -n toxindex-app
kubectl rollout restart deployment/backend -n toxindex-app
```

## 📊 Performance

### Build Times
- **Backend**: 3-5 minutes (Docker build + K8s deploy)
- **Frontend**: 2-3 minutes (Node.js build + GCS sync)
- **Full Stack**: 5-8 minutes (both deployments)

### Resource Usage
- **Machine Type**: E2_HIGHCPU_8 (8 vCPUs, 8GB RAM)
- **Disk Size**: 100GB for full stack, 50GB for frontend
- **Concurrent Builds**: Up to 10 (configurable)

## 🎉 Getting Started

1. **Set up the project**:
   ```bash
   gcloud config set project toxindex
   ```

2. **Enable APIs**:
   ```bash
   gcloud services enable cloudbuild.googleapis.com containerregistry.googleapis.com artifactregistry.googleapis.com container.googleapis.com
   ```

3. **Configure permissions**:
   ```bash
   # Run the service account setup commands above
   ```

4. **Create triggers**:
   ```bash
   gcloud builds triggers import --source=cloudbuild-triggers.yaml
   ```

5. **Test deployment**:
   ```bash
   git push origin main
   ```

Your Cloud Build CI/CD pipeline is now ready! 🚀 