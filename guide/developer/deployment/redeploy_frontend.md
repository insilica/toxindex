---
title: redeploy frontend
sidebar_position: 2
---

### Step 1. add a field to resources>default_workflows.json

```
{
  "workflow_id": 6,
  "frontend_id": "toxindex-sygma",
  "title": "toxindex-sygma",
  "label": "Sygma Analysis",
  "description": "Sygma analysis is blah blah.",
  "initial_prompt": "Enter pathway ID (e.g., WP3657) or upload data file",
  "celery_task": "metabolite-sygma"
}
```

### Step 2. update database

```bash
python scripts/seed_workflows.py
```

### Step 2. rebuild static files

```bash
cd frontend && npm run build
```

### Step 3. update the google cloud storage

```bash
gsutil -m rsync -r dist/ gs://toxindex-react
```

### Step 4. remove cache to make the change effective

```bash
gsutil -m setmeta -h "Cache-Control:no-store" gs://toxindex-react/**
gcloud compute url-maps invalidate-cdn-cache frontend-url-map --path "/*"
```