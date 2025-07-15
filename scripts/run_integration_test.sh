#!/usr/bin/env bash

# Integration test runner for toxindex
# This script sets up the environment and runs the full integration test

set -euo pipefail

echo "🧪 Starting toxindex integration test..."
echo "========================================"

# Check if we're in the right directory
if [ ! -f "pyproject.toml" ]; then
    echo "❌ Error: Please run this script from the toxindex project root"
    exit 1
fi

# Check if virtual environment is activated
if [ -z "${VIRTUAL_ENV:-}" ]; then
    echo "❌ Error: Please activate the virtual environment first"
    echo "   Run: source .venv/bin/activate"
    exit 1
fi

# Check if required services are running
echo "🔍 Checking required services..."

# Check Redis
if ! redis-cli ping > /dev/null 2>&1; then
    echo "❌ Error: Redis is not running"
    echo "   Please start Redis first: redis-server"
    exit 1
fi
echo "✅ Redis is running"

# Check PostgreSQL
if ! pg_isready -h localhost -p 5433 > /dev/null 2>&1; then
    echo "❌ Error: PostgreSQL is not running on port 5433"
    echo "   Please start PostgreSQL first"
    exit 1
fi
echo "✅ PostgreSQL is running"

# Check if Blazegraph is running
if ! curl -s http://localhost:8889/bigdata/namespace/kb/sparql > /dev/null 2>&1; then
    echo "⚠️  Warning: Blazegraph is not running on port 8889"
    echo "   This may affect pathway analysis features"
fi

# Run the integration test
echo "🚀 Running integration test..."
python tests/test_integration.py

# Check the exit code
if [ $? -eq 0 ]; then
    echo ""
    echo "🎉 Integration test completed successfully!"
    echo "✅ All services are working correctly"
    echo "✅ Login with test@test.com works"
    echo "✅ Workflow creation with 'is acetaminophen toxic?' works"
else
    echo ""
    echo "❌ Integration test failed!"
    echo "Please check the logs above for details"
    exit 1
fi 