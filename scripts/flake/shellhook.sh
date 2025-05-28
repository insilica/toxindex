#!/usr/bin/env bash

source scripts/flake/setup_postgres.sh \
  "postgres" \
  "devpassword" \
  "5433" \
  "toxindex" \
  "localhost" \
  "socket" \
  "/nix/store/..." # use actual path or pass as env

echo -e "\nLogging in to AWS..."
source scripts/flake/aws_login.sh

echo -e "\nRunning flyway migrations..."
source scripts/flake/run_flyway.sh

echo -e "\nStarting redis..."
source scripts/flake/start_redis.sh

echo -e "\nLoading application environment script..."
source scripts/load_environment.sh "$AWS_PROFILE" "insilica/toxindex+development"

echo -e "\nDevelopment environment ready."
export FLASK_APP=webserver.app
export FLASK_ENV=development
export DEBUG=1
export PREFERRED_URL_SCHEME=http
export SERVER_NAME=localhost:6513
source .env
export PYTHONPATH="$PWD"
