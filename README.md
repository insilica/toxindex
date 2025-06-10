# Toxindex Portal
A collection of toxicology workflows

## Services
1. postgres
2. flyway - runs at the start of docker compose up and performs all necessary migrations on postgres service.
3. webserver - a flask web user interface, it has users, projects, and projects have views for each services

## Development
deactivate
rm -rf .venv
rm -rf ~/.cache/uv
unset LD_LIBRARY_PATH
nix flake update
nix develop

1. nix develop

2-1 (production) gunicorn webserver.app:app --bind 0.0.0.0:8000 --worker-class eventlet
2-2 (production) .venv/bin/python redis_listener_service.py
2. python -m webserver.app
3. celery -A workflows.celery_worker worker --loglevel=info

lsof -i :8000
kill PID

git add . && git commit -m "frontend update" && git push

## Roadmap
The basic idea is that: 

1. Users create environments
2. Users execute workflows in their environments

We should drive towards a ProbRA workflow.see tasks in pfoa_tasks.md. We (WP5) were asked to:

1. predict hazard with toxtransformer

To do this we can:

1. create a workflow run for probra
2. user puts in a prompt "Is PFOA hepatotoxic?"
3. create a task, for now there is no back and forth. 
    1. Adds a task to a table of tasks under the prompt
    2. Send task to workflow
    3. Create a webhook for activities and results

To get this working we can:

1. create a docker container that 





Frontend
React, Tainwind CSS with Vite, TypeScript

Backend
Flask, Eventlet, Celery

Data
Langchain

Gunicorn on EC2