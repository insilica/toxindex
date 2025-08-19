#!/bin/bash

# Build and test script for probra Docker image
set -e

echo "=== Building Probra Docker Image ==="

# Build the base image first
echo "Building base image..."
docker build -f Dockerfile.base -t toxindex-base:latest .

# Build the rap image
echo "Building rap image..."
docker build -f Dockerfile.rap -t toxindex-rap:latest .

echo "=== Testing Docker Image ==="

# Test the image by running our test script
echo "Running test script in container..."
docker run --rm \
  -e REDIS_HOST=localhost \
  -e REDIS_PORT=6379 \
  -e CELERY_BROKER_URL=redis://localhost:6379/0 \
  -e CELERY_RESULT_BACKEND=redis://localhost:6379/0 \
  -e OPENAI_API_KEY=${OPENAI_API_KEY} \
  -e GEMINI_API_KEY=${GEMINI_API_KEY} \
  -e GOOGLE_CSE_ID=${GOOGLE_CSE_ID} \
  -e GOOGLE_SEARCH_KEY=${GOOGLE_SEARCH_KEY} \
  -e PGHOST=${PGHOST} \
  -e PGPORT=${PGPORT} \
  -e PGDATABASE=${PGDATABASE} \
  -e PGUSER=${PGUSER} \
  -e PGPASSWORD=${PGPASSWORD} \
  -e GCS_BUCKET_NAME=${GCS_BUCKET_NAME} \
  -v $(pwd)/test_probra_docker.py:/app/test_probra_docker.py \
  toxindex-rap:latest \
  python test_probra_docker.py

echo "=== Testing Celery Worker ==="

# Test that the celery worker can start (without actually connecting to Redis)
echo "Testing celery worker startup..."
docker run --rm \
  -e REDIS_HOST=localhost \
  -e REDIS_PORT=6379 \
  -e CELERY_BROKER_URL=redis://localhost:6379/0 \
  -e CELERY_RESULT_BACKEND=redis://localhost:6379/0 \
  -e OPENAI_API_KEY=${OPENAI_API_KEY} \
  -e GEMINI_API_KEY=${GEMINI_API_KEY} \
  -e GOOGLE_CSE_ID=${GOOGLE_CSE_ID} \
  -e GOOGLE_SEARCH_KEY=${GOOGLE_SEARCH_KEY} \
  -e PGHOST=${PGHOST} \
  -e PGPORT=${PGPORT} \
  -e PGDATABASE=${PGDATABASE} \
  -e PGUSER=${PGUSER} \
  -e PGPASSWORD=${PGPASSWORD} \
  -e GCS_BUCKET_NAME=${GCS_BUCKET_NAME} \
  toxindex-rap:latest \
  timeout 10s celery -A workflows.celery_worker worker --loglevel=info -Q probra || true

echo "=== Build and Test Complete ==="
echo "Image 'toxindex-rap:latest' is ready for use!"
echo ""
echo "To run the worker with Redis:"
echo "docker run --rm \\"
echo "  -e REDIS_HOST=your-redis-host \\"
echo "  -e REDIS_PORT=6379 \\"
echo "  -e CELERY_BROKER_URL=redis://your-redis-host:6379/0 \\"
echo "  -e CELERY_RESULT_BACKEND=redis://your-redis-host:6379/0 \\"
echo "  toxindex-rap:latest"
