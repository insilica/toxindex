#!/usr/bin/env bash

if ! pgrep redis-server > /dev/null; then
  echo "Starting Redis server..."
  redis-server --save "" --stop-writes-on-bgsave-error no --daemonize yes
  sleep 1
  if redis-cli ping > /dev/null 2>&1; then
    echo "Redis started successfully."
  else
    echo "Redis failed to start!" >&2
    exit 1
  fi
else
  echo "Redis is already running."
fi