#!/usr/bin/env bash
set -euo pipefail

if ! pgrep redis-server > /dev/null; then
  echo "Starting Redis server..."
  redis-server --daemonize yes
else
  echo "Redis is already running."
fi
