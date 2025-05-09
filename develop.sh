#!/bin/bash

# --- Configuration ---
# Docker Compose file location (relative to this script)
COMPOSE_FILE="docker-compose.yml" # Assumes docker-compose.yml is in the same directory

# Service names defined in docker-compose.yml
DB_SERVICE_NAME="db"
FLYWAY_SERVICE_NAME="flyway"
WEB_SERVICE_NAME="webserver"

# --- Cleanup Function ---
# This function will be called when the script is interrupted (e.g., Ctrl+C)
cleanup() {
    echo
    echo "Interrupt received, stopping and removing Docker Compose services..."
    # Stops and removes containers, networks, and volumes (if not named externally)
    docker-compose -f "$COMPOSE_FILE" down -v --remove-orphans
    echo "Cleanup complete."
    exit 0
}

# --- Trap Interrupts ---
# Call the cleanup function on SIGINT (Ctrl+C) and SIGTERM
trap cleanup SIGINT SIGTERM

# --- Start PostgreSQL using Docker Compose ---
echo "Starting PostgreSQL service ($DB_SERVICE_NAME) via Docker Compose..."
# Start in detached mode
docker-compose -f "$COMPOSE_FILE" up -d "$DB_SERVICE_NAME"

# Check if the DB service started successfully
# A more robust check would be to wait for the healthcheck if defined
echo "Waiting for $DB_SERVICE_NAME to be healthy..."
# This loop attempts to wait for the service to be healthy.
# Adjust timeout as needed.
timeout_seconds=60
start_time=$(date +%s)
while true; do
    current_time=$(date +%s)
    elapsed_time=$((current_time - start_time))

    if [ "$elapsed_time" -ge "$timeout_seconds" ]; then
        echo "Timeout waiting for $DB_SERVICE_NAME to become healthy. Please check logs: docker-compose -f \"$COMPOSE_FILE\" logs $DB_SERVICE_NAME"
        cleanup
        exit 1
    fi

    # Check health status (requires healthcheck in docker-compose.yml)
    # If no healthcheck, you might check for logs indicating readiness or try to connect.
    health_status=$(docker-compose -f "$COMPOSE_FILE" ps -q "$DB_SERVICE_NAME" | xargs -I {} docker inspect -f '{{if .State.Health}}{{.State.Health.Status}}{{else}}{{.State.Status}}{{end}}' {})
    
    if [ "$health_status" == "healthy" ]; then
        echo "$DB_SERVICE_NAME is healthy."
        break
    elif [ "$health_status" == "running" ] && ! docker inspect --format='{{.State.Health}}' "$(docker-compose -f "$COMPOSE_FILE" ps -q "$DB_SERVICE_NAME")" | grep -q "Status"; then
        # If no healthcheck defined but service is running, assume ready after a short delay (less reliable)
        echo "$DB_SERVICE_NAME is running (no healthcheck, proceeding with caution after short delay)."
        sleep 5 # Adjust as needed
        break
    elif [[ "$health_status" == "exited" || "$health_status" == "dead" ]]; then
        echo "$DB_SERVICE_NAME failed to start. Status: $health_status. Please check logs: docker-compose -f \"$COMPOSE_FILE\" logs $DB_SERVICE_NAME"
        cleanup
        exit 1
    fi
    
    echo "Still waiting for $DB_SERVICE_NAME... (Status: $health_status)"
    sleep 5
done


# --- Run Flyway Migrations using Docker Compose ---
echo "Running Flyway migrations ($FLYWAY_SERVICE_NAME) via Docker Compose..."
# 'run --rm' creates a temporary container for the command and removes it afterward
docker-compose -f "$COMPOSE_FILE" run --rm "$FLYWAY_SERVICE_NAME" migrate

# Check if Flyway command was successful
if [ $? -ne 0 ]; then
    echo "Flyway migrations failed. Exiting."
    cleanup # Or handle error differently
    exit 1
fi
echo "Flyway migrations completed."

# --- Start Flask Webserver using Docker Compose ---
echo "Starting Flask webserver service ($WEB_SERVICE_NAME) via Docker Compose..."
echo "Flask app will be available on http://localhost:5000 (Press CTRL+C to stop all services)"
# Start the webserver in the foreground to see logs directly.
# The 'trap' will handle cleanup when you press Ctrl+C.
docker-compose -f "$COMPOSE_FILE" up "$WEB_SERVICE_NAME"

# If 'docker-compose up' for the webserver exits (e.g., Flask crashes),
# the script will reach here.
echo "Webserver process ended."
cleanup # Perform cleanup
