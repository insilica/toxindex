#!/bin/bash

# Build script for multi-image Docker setup
set -e

echo "🚀 Building ToxIndex Docker images..."

# Build base image first
echo "📦 Building base image..."
docker build -f Dockerfile.base -t toxindex-base:latest .

# Build general worker image
echo "🔧 Building general worker image..."
docker build -f Dockerfile.general -t toxindex-general:latest .

# Build RAPtool worker image
echo "🧪 Building RAPtool worker image..."
docker build -f Dockerfile.raptool -t toxindex-raptool:latest .

# Build pathway worker image
echo "🔄 Building pathway worker image..."
docker build -f Dockerfile.pathway -t toxindex-pathway:latest .

echo "✅ All images built successfully!"
echo ""
echo "📋 Available images:"
echo "  - toxindex-base:latest (base image)"
echo "  - toxindex-general:latest (general tasks)"
echo "  - toxindex-raptool:latest (RAPtool tasks)"
echo "  - toxindex-pathway:latest (pathway tasks)"
echo ""
echo "🚀 To push to Google Container Registry:"
echo "  docker tag toxindex-general:latest us-docker.pkg.dev/toxindex/toxindex-backend/general:latest"
echo "  docker tag toxindex-raptool:latest us-docker.pkg.dev/toxindex/toxindex-backend/raptool:latest"
echo "  docker tag toxindex-pathway:latest us-docker.pkg.dev/toxindex/toxindex-backend/pathway:latest"
echo ""
echo "  docker push us-docker.pkg.dev/toxindex/toxindex-backend/general:latest"
echo "  docker push us-docker.pkg.dev/toxindex/toxindex-backend/raptool:latest"
echo "  docker push us-docker.pkg.dev/toxindex/toxindex-backend/pathway:latest"
