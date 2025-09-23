#!/bin/bash
# make executable with: chmod +x dev_startup.sh

# # Enter nix environment
# nix develop || exit 1

# Start Flask app
python -m webserver.app &

# Start redis listener
python redis_listener_standalone.py &

# Start Celery workers
celery -A workflows.celery_worker worker --loglevel=info -Q celery,probra,openai,raptool,pathway,sygma &
# celery -A workflows.celery_worker_probra worker --loglevel=info -Q probra &
# celery -A workflows.celery_worker_raptool worker --loglevel=info -Q raptool &
# celery -A workflows.celery_worker_ranking worker --loglevel=info -Q ranking &
celery -A workflows.celery_worker_sygma worker --loglevel=info -Q sygma &

# Start frontend (development mode)
(cd frontend && npm run dev) &