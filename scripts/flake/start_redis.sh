#!/usr/bin/env bash

if ! pgrep redis-server > /dev/null; then
  echo "Starting Redis server..."
  redis-server --save "" --stop-writes-on-bgsave-error no --daemonize yes
  sleep 1
  redis-cli ping > /dev/null 2>&1 && echo "Redis started successfully." || echo "Redis failed to start!"
else
  echo "Redis is already running."
fi