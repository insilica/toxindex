#!/usr/bin/env bash

OLD_OPTS=$(set +o)
set -euo pipefail

source scripts/flake/setup_postgres.sh \
  "postgres" \
  "devpassword" \
  "5433" \
  "toxindex" \
  "localhost" \
  "socket" \
  "/nix/store/3pzlrs5nddszkpgasnrcpf4ifrzm76lb-postgresql-15.13/bin"

# Drop and recreate the public schema (dev only!)
# psql -U postgres -d toxindex -c "DROP SCHEMA public CASCADE; CREATE SCHEMA public;"

source scripts/flake/run_flyway.sh

# Seed default workflows from JSON
echo "Seeding default workflows..."
python3 scripts/seed_workflows.py

# Set up default admin user
echo "Setting up default admin user..."
python3 scripts/setup_default_admin.py

source scripts/flake/start_redis.sh

source scripts/load_gcp_environment.sh 

eval "$OLD_OPTS"

export FLASK_APP=webserver.app
# export FLASK_ENV=development
# export FLASK_DEBUG=1
export FLASK_ENV=production
export FLASK_DEBUG=0
export PREFERRED_URL_SCHEME=http

# export SERVER_NAME=localhost:6513
export SERVER_NAME=https://toxindex.com
source .env

export FRONTEND_URL=https://toxindex.com

# Redis (Redis is running on the same EC2 instance as the webserver)
export CELERY_BROKER_URL=redis://localhost:6379/0
export CELERY_RESULT_BACKEND=redis://localhost:6379/0
export SOCKETIO_MESSAGE_QUEUE=redis://localhost:6379/0

# Python
export PYTHONPATH="${PYTHONPATH:+${PYTHONPATH}:}${PWD}"