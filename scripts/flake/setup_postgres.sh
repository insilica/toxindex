#!/usr/bin/env bash
set -euo pipefail

# Arguments:
#   1 - PostgreSQL user
#   2 - PostgreSQL password
#   3 - PostgreSQL port
#   4 - Database name
#   5 - Host
#   6 - Socket directory subpath
#   7 - PostgreSQL bin path

PG_SETTINGS_USER="$1"
PG_SETTINGS_PASSWORD="$2"
PG_SETTINGS_PORT="$3"
PG_SETTINGS_DB_NAME="$4"
PG_SETTINGS_HOST="$5"
PG_SETTINGS_SOCKET_DIR_SUBPATH="$6"
POSTGRESQL_BIN_PATH="$7"

export PG_SETTINGS_USER PG_SETTINGS_PASSWORD PG_SETTINGS_PORT \
       PG_SETTINGS_DB_NAME PG_SETTINGS_HOST PG_SETTINGS_SOCKET_DIR_SUBPATH \
       POSTGRESQL_BIN_PATH

# Call the existing PostgreSQL helper
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$SCRIPT_DIR/.."

bash "$ROOT_DIR/shellhook_postgres.sh"
