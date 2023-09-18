#!/usr/bin/env bash
set -euo pipefail

cd "$DEVENV_ROOT/services/report"

# rdkit python library needs this for libstdc++.so.6
export LD_LIBRARY_PATH="$__STDCXX_PATH:$LD_LIBRARY_PATH"

poetry lock
poetry install

poetry run \
	flask --app report.app --debug \
	run --host=0.0.0.0 --port="$REPORTS_PORT"
