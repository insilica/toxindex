#!/usr/bin/env bash
set -euo pipefail

echo "Setting up PostgreSQL development environment..."

export PGDATA="$PWD/.postgres/data"
export PG_SOCKET_DIR="$PWD/.postgres/socket"
export PGUSER="postgres"
export PGPASSWORD="devpassword"
export PGHOST="localhost"
export PGPORT="5433"
export PGDATABASE="toxindex"

PGLOGFILE="$PWD/.postgres/logfile"
mkdir -p "$PWD/.postgres"
mkdir -p "$PG_SOCKET_DIR"

if [ ! -f "$PGDATA/postgresql.conf" ]; then
  echo "Initializing new PostgreSQL instance in $PGDATA..."
  mkdir -p "$PGDATA"
  initdb --username="$PGUSER" --pwfile=<(echo "$PGPASSWORD") -D "$PGDATA"

  echo "listen_addresses = 'localhost'" >> "$PGDATA/postgresql.conf"
  echo "port = 5433" >> "$PGDATA/postgresql.conf"
  echo "unix_socket_directories = '$PG_SOCKET_DIR'" >> "$PGDATA/postgresql.conf"
fi

if ! pg_ctl -D "$PGDATA" status > /dev/null 2>&1; then
  echo "Starting PostgreSQL server (Log: $PGLOGFILE)..."
  if ! pg_ctl -D "$PGDATA" -l "$PGLOGFILE" -w -t 30 start; then
    echo "--------------------------------------------------------------------"
    echo "PostgreSQL FAILED to start. Contents of $PGLOGFILE:"
    cat "$PGLOGFILE"
    echo "--------------------------------------------------------------------"
    exit 1
  fi
  echo "PostgreSQL server started successfully."
else
  echo "PostgreSQL already running."
fi

if ! psql -lqt | cut -d \| -f 1 | grep -qw "$PGDATABASE"; then
  echo "Creating database: $PGDATABASE..."
  createdb "$PGDATABASE"
else
  echo "Database $PGDATABASE already exists."
fi

echo "Running AWS login helper..."
bash scripts/shellhook_aws_login.sh

echo "Running Flyway migrations..."
flyway \
  -url="jdbc:postgresql://$PGHOST:$PGPORT/$PGDATABASE" \
  -user="$PGUSER" \
  -password="$PGPASSWORD" \
  -locations="filesystem:./flyway/sql" \
  migrate

echo "Loading application environment script..."
bash scripts/load_environment.sh "$AWS_PROFILE" "insilica/toxindex+development"

echo "Development environment ready."
echo "PostgreSQL is accessible on: $PGHOST:$PGPORT (TCP/IP)"
echo "PostgreSQL server socket directory: $PG_SOCKET_DIR"
echo "Database: $PGDATABASE, User: $PGUSER"
