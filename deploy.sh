#!/usr/bin/env bash

set -e  # Exit on error

echo "Pulling latest code..."
git pull

echo "Building frontend..."
cd frontend
npm install
npm run build
cd ..

echo "Restarting backend and celery services..."
sudo systemctl restart toxindex-backend
sudo systemctl restart toxindex-celery

echo "Deployment complete!"