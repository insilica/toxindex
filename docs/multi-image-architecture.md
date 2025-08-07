# Multi-Image Docker Architecture

## Overview

This document describes the multi-image Docker architecture that provides **complete environment isolation** for different task types, solving Python dependency conflicts in production.

## Architecture Benefits

### ✅ **Complete Dependency Isolation**
- Each worker type has its own Docker image with specific dependencies
- No conflicts between external repositories
- Independent version management per task type

### ✅ **Production Best Practices**
- Follows microservices principles
- Enables independent scaling and deployment
- Improves security through isolation

### ✅ **Resource Optimization**
- Smaller, focused images per task type
- Faster builds and deployments
- Better resource utilization

## Image Structure

```
┌─────────────────────────────────────────────────────────────┐
│                    Docker Images                            │
├─────────────────────────────────────────────────────────────┤
│  ┌─────────────┐  ┌─────────────┐  ┌─────────────┐       │
│  │   Base      │  │   General   │  │   RAPtool   │       │
│  │   Image     │  │   Image     │  │   Image     │       │
│  │ (Common)    │  │ (Standard)  │  │ (RAPtool)   │       │
│  └─────────────┘  └─────────────┘  └─────────────┘       │
│           │              │              │                 │
│           └──────────────┼──────────────┘                 │
│                          │                                │
│                    ┌─────────────┐                        │
│                    │   Pathway   │                        │
│                    │   Image     │                        │
│                    │ (Pathway)   │                        │
│                    └─────────────┘                        │
└─────────────────────────────────────────────────────────────┘
```

## Image Details

### Base Image (`Dockerfile.base`)
- **Purpose**: Common dependencies and setup
- **Contains**: System dependencies, Python setup, project files
- **Used by**: All other images

### General Image (`Dockerfile.general`)
- **Purpose**: Standard tasks (probra, openai, default)
- **Dependencies**: Core application dependencies
- **Entrypoint**: Gunicorn for web server

### RAPtool Image (`Dockerfile.raptool`)
- **Purpose**: RAPtool-specific tasks
- **Dependencies**: RAPtool + external repositories
- **Entrypoint**: Celery worker for raptool queue

### Pathway Image (`Dockerfile.pathway`)
- **Purpose**: Pathway analysis tasks
- **Dependencies**: Pathway analysis + external repositories
- **Entrypoint**: Celery worker for pathway queue

## Requirements Files

### `requirements-general.txt`
- Core application dependencies
- Standard ML/MLOps packages
- No external repository conflicts

### `requirements-raptool.txt`
- RAPtool-specific dependencies
- **Add your external repository here**:
```txt
# Example:
git+https://github.com/example/external-repo.git@v1.0.0#egg=external-repo
```

### `requirements-pathway.txt`
- Pathway analysis dependencies
- **Add your external repository here**:
```txt
# Example:
git+https://github.com/example/pathway-repo.git@v1.0.0#egg=pathway-repo
```

## Building Images

### Quick Build
```bash
# Build all images
./scripts/build-images.sh
```

### Manual Build
```bash
# Build base image first
docker build -f Dockerfile.base -t toxindex-base:latest .

# Build specific images
docker build -f Dockerfile.general -t toxindex-general:latest .
docker build -f Dockerfile.raptool -t toxindex-raptool:latest .
docker build -f Dockerfile.pathway -t toxindex-pathway:latest .
```

### Push to Google Container Registry
```bash
# Tag images
docker tag toxindex-general:latest us-docker.pkg.dev/toxindex/toxindex-backend/general:latest
docker tag toxindex-raptool:latest us-docker.pkg.dev/toxindex/toxindex-backend/raptool:latest
docker tag toxindex-pathway:latest us-docker.pkg.dev/toxindex/toxindex-backend/pathway:latest

# Push images
docker push us-docker.pkg.dev/toxindex/toxindex-backend/general:latest
docker push us-docker.pkg.dev/toxindex/toxindex-backend/raptool:latest
docker push us-docker.pkg.dev/toxindex/toxindex-backend/pathway:latest
```

## Adding External Repositories

### Step 1: Identify the Task Type
- **RAPtool tasks**: Add to `requirements-raptool.txt`
- **Pathway tasks**: Add to `requirements-pathway.txt`
- **General tasks**: Add to `requirements-general.txt`

### Step 2: Add with Version Constraints
```txt
# In requirements-raptool.txt
git+https://github.com/example/external-repo.git@v1.0.0#egg=external-repo
```

### Step 3: Rebuild the Specific Image
```bash
# Rebuild only the affected image
docker build -f Dockerfile.raptool -t toxindex-raptool:latest .
docker push us-docker.pkg.dev/toxindex/toxindex-backend/raptool:latest
```

### Step 4: Update Kubernetes Deployment
```bash
# Restart the specific deployment
kubectl rollout restart deployment celery-worker-raptool
```

## Kubernetes Integration

### Deployment Files
- `k8s/celery-worker-general-deployment.yaml` → `general:latest`
- `k8s/celery-worker-raptool-deployment.yaml` → `raptool:latest`
- `k8s/celery-worker-pathway-deployment.yaml` → `pathway:latest`

### Deployment Commands
```bash
# Deploy all workers
kubectl apply -f k8s/celery-worker-general-deployment.yaml
kubectl apply -f k8s/celery-worker-raptool-deployment.yaml
kubectl apply -f k8s/celery-worker-pathway-deployment.yaml

# Scale specific workers
kubectl scale deployment celery-worker-raptool --replicas=2
kubectl scale deployment celery-worker-pathway --replicas=1
```

## Monitoring and Debugging

### Check Image Sizes
```bash
docker images | grep toxindex
```

### Inspect Image Contents
```bash
# Check what's installed in each image
docker run --rm toxindex-raptool:latest pip list
docker run --rm toxindex-pathway:latest pip list
```

### Debug Build Issues
```bash
# Build with verbose output
docker build -f Dockerfile.raptool -t toxindex-raptool:latest . --progress=plain --no-cache
```

## Migration Strategy

### Phase 1: Build and Test
1. Build all images locally
2. Test each image individually
3. Verify dependency isolation

### Phase 2: Deploy Gradually
1. Deploy new workers alongside existing ones
2. Route a small percentage of tasks to new workers
3. Monitor performance and stability

### Phase 3: Full Migration
1. Route all tasks to new workers
2. Scale down old workers
3. Remove old deployment

## Troubleshooting

### Common Issues

#### Build Failures
```bash
# Check if base image exists
docker images | grep toxindex-base

# Rebuild base image first
docker build -f Dockerfile.base -t toxindex-base:latest .
```

#### Dependency Conflicts
```bash
# Check installed packages in image
docker run --rm toxindex-raptool:latest pip check

# Test specific package installation
docker run --rm toxindex-raptool:latest python -c "import your_package"
```

#### Kubernetes Issues
```bash
# Check pod logs
kubectl logs -f deployment/celery-worker-raptool

# Check image pull issues
kubectl describe pod -l app=celery-worker-raptool
```

### Performance Optimization

#### Image Size Optimization
- Use multi-stage builds
- Remove unnecessary dependencies
- Use `.dockerignore` to exclude files

#### Build Time Optimization
- Use Docker layer caching
- Build base image once
- Use build cache effectively

## Best Practices

1. **Always specify versions** for external repositories
2. **Test images locally** before pushing to registry
3. **Monitor resource usage** per worker type
4. **Use health checks** in Kubernetes deployments
5. **Implement proper logging** for each worker type
6. **Set up monitoring** for queue lengths and processing times
