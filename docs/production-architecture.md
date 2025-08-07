# Production Architecture Guide

## Overview

This document outlines the recommended production architecture for handling Python dependency conflicts in GKE.

## Problem Statement

When adding external repositories to `requirements.txt`, you may encounter dependency conflicts that cannot be resolved by Nix or uv. This is because:

1. **Nix** only provides the Python environment and system dependencies
2. **uv** handles Python package resolution, but conflicts arise when packages have incompatible requirements
3. **Single pod architecture** forces all tasks to share the same dependency environment

## Solution: Multi-Pod Architecture

### Architecture Overview

```
┌─────────────────────────────────────────────────────────────┐
│                    GKE Cluster                             │
├─────────────────────────────────────────────────────────────┤
│  ┌─────────────┐  ┌─────────────┐  ┌─────────────┐       │
│  │   General   │  │   RAPtool   │  │  Pathway    │       │
│  │   Worker    │  │   Worker    │  │   Worker    │       │
│  │   Pods      │  │   Pod       │  │   Pod       │       │
│  │ (2 replicas)│  │ (1 replica) │  │ (1 replica) │       │
│  └─────────────┘  └─────────────┘  └─────────────┘       │
│           │              │              │                 │
│           └──────────────┼──────────────┘                 │
│                          │                                │
│                    ┌─────────────┐                        │
│                    │   Redis     │                        │
│                    │  (Broker)   │                        │
│                    └─────────────┘                        │
└─────────────────────────────────────────────────────────────┘
```

### Queue Assignment

| Task Type | Queue | Worker Pod | Dependencies |
|-----------|-------|------------|--------------|
| `probra_task` | `probra` | General | Standard ML stack |
| `plain_openai_task` | `openai` | General | OpenAI, basic deps |
| `raptool_task` | `raptool` | RAPtool | RAPtool-specific deps |
| `pathway_analysis_task` | `pathway` | Pathway | Pathway-specific deps |

## Implementation Steps

### 1. Update Task Definitions

All tasks now specify their queue:

```python
@celery.task(bind=True, queue='raptool')
def raptool_task(self, payload):
    # RAPtool specific code
```

### 2. Deploy Separate Worker Pods

```bash
# Deploy the new worker pods
kubectl apply -f k8s/celery-worker-raptool-deployment.yaml
kubectl apply -f k8s/celery-worker-pathway-deployment.yaml
kubectl apply -f k8s/celery-worker-general-deployment.yaml

# Scale down the old single worker
kubectl scale deployment celery-worker --replicas=0
```

### 3. Handle External Repository Dependencies

For each external repository with conflicts:

#### Option A: Separate Requirements Files

Create separate requirements files for each worker type:

```bash
# requirements-raptool.txt
-e ./RAPtool
# Add RAPtool-specific external deps here

# requirements-pathway.txt  
-e ./pathway_analysis_tool
# Add pathway-specific external deps here

# requirements-general.txt
# Standard dependencies only
```

#### Option B: Docker Multi-Stage Builds

Create separate Docker images for each worker type:

```dockerfile
# Dockerfile.raptool
FROM python:3.12-slim
COPY requirements-raptool.txt .
RUN pip install -r requirements-raptool.txt
# ... rest of setup
```

### 4. Update Task Routing

Ensure tasks are routed to the correct queues:

```python
# In your task submission code
celery.send_task('workflows.raptool_task', args=[payload], queue='raptool')
celery.send_task('workflows.pathway_analysis_task', args=[payload], queue='pathway')
```

## Benefits

1. **Dependency Isolation**: Each worker type has its own dependency environment
2. **Resource Optimization**: Scale workers based on task type demand
3. **Fault Isolation**: Issues in one task type don't affect others
4. **Independent Scaling**: Scale different task types independently
5. **Easier Maintenance**: Update dependencies for specific task types without affecting others

## Monitoring and Scaling

### Resource Monitoring

Monitor each worker type separately:

```bash
# Check queue lengths
kubectl exec -it <pod-name> -- celery -A workflows.celery_worker inspect active_queues

# Monitor worker status
kubectl exec -it <pod-name> -- celery -A workflows.celery_worker inspect stats
```

### Auto-scaling

Set up Horizontal Pod Autoscalers for each worker type:

```yaml
apiVersion: autoscaling/v2
kind: HorizontalPodAutoscaler
metadata:
  name: celery-worker-general-hpa
spec:
  scaleTargetRef:
    apiVersion: apps/v1
    kind: Deployment
    name: celery-worker-general
  minReplicas: 1
  maxReplicas: 10
  metrics:
  - type: Resource
    resource:
      name: cpu
      target:
        type: Utilization
        averageUtilization: 70
```

## Migration Strategy

1. **Phase 1**: Deploy new worker pods alongside existing ones
2. **Phase 2**: Gradually migrate tasks to new queues
3. **Phase 3**: Monitor and scale down old workers
4. **Phase 4**: Remove old worker deployment

## Troubleshooting

### Common Issues

1. **Tasks stuck in queue**: Check if workers are consuming from correct queues
2. **Memory issues**: Adjust resource limits per worker type
3. **Dependency conflicts**: Ensure each worker has correct requirements file

### Debugging Commands

```bash
# Check worker queues
kubectl exec -it <pod-name> -- celery -A workflows.celery_worker inspect active

# Check task routing
kubectl exec -it <pod-name> -- celery -A workflows.celery_worker inspect registered

# Monitor Redis queues
kubectl exec -it redis-pod -- redis-cli LLEN celery
```
