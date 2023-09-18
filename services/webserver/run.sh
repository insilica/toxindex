#!/usr/bin/env bash
set -euo pipefail

cd "$DEVENV_ROOT/services/webserver"

poetry lock
poetry install

poetry run \
	flask --app webserver.app --debug \
	run --host=0.0.0.0 --port="$WEBSERVER_PORT"
