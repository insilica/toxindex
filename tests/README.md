# Toxindex Integration Tests

This directory contains integration tests for the toxindex fullstack application.

## Test Overview

The integration test verifies that all components of the toxindex stack work together correctly:

1. **Celery Worker** - Background task processing
2. **Flask Webapp** - Backend API server
3. **React Frontend** - User interface
4. **Authentication** - Login with test@test.com
5. **Workflow Creation** - Creating tasks with the toxindex-rap workflow
6. **Task Execution** - Running "is acetaminophen toxic?" query

## Prerequisites

Before running the integration test, ensure you have:

1. **Virtual environment activated**:
   ```bash
   source .venv/bin/activate
   ```

2. **Required services running**:
   - Redis: `redis-server`
   - PostgreSQL on port 5433
   - Blazegraph on port 8889 (optional, for pathway analysis)

3. **Dependencies installed**:
   ```bash
   pip install requests
   ```

## Running the Test

### Option 1: Using the test runner script (Recommended)

```bash
./run_integration_test.sh
```

This script will:
- Check if all required services are running
- Install test dependencies
- Run the full integration test
- Provide clear success/failure feedback

### Option 2: Running the test directly

```bash
python tests/test_integration.py
```

## Test Details

### What the test does:

1. **Starts Celery Worker**: Launches a Celery worker process for background task processing
2. **Starts Flask Webapp**: Launches the Flask backend server on port 6513
3. **Starts React Frontend**: Launches the React development server on port 5173
4. **Logs in**: Authenticates with test@test.com / test credentials
5. **Creates Workflow**: Creates a new task using the toxindex-rap workflow
6. **Submits Query**: Sends "is acetaminophen toxic?" to the workflow
7. **Monitors Progress**: Tracks the task status until completion
8. **Cleans up**: Stops all started processes

### Expected Output:

```
🧪 Starting toxindex integration test...
========================================
🔍 Checking required services...
✅ Redis is running
✅ PostgreSQL is running
⚠️  Warning: Blazegraph is not running on port 8889
📦 Installing test dependencies...
🚀 Running integration test...
=== Step 1: Starting Celery worker ===
✅ Celery worker started
=== Step 2: Starting Flask webapp ===
✅ Flask webapp started
=== Step 3: Starting React frontend ===
✅ React frontend started
=== Step 4: Logging in with test@test.com ===
✅ Login successful
=== Step 5: Getting available workflows ===
✅ Found workflow: toxindex-rap (ID: 1)
=== Step 6: Creating task with 'is acetaminophen toxic?' ===
✅ Task created with ID: abc123...
=== Step 7: Monitoring task status ===
Task status: starting
Task status: running
Task status: done
✅ Task completed successfully
🎉 All integration tests passed!

🎉 Integration test completed successfully!
✅ All services are working correctly
✅ Login with test@test.com works
✅ Workflow creation with 'is acetaminophen toxic?' works
```

## Troubleshooting

### Common Issues:

1. **Redis not running**:
   ```bash
   redis-server
   ```

2. **PostgreSQL not running**:
   ```bash
   # Start PostgreSQL (check your setup)
   pg_ctl -D /path/to/postgres/data start
   ```

3. **Port conflicts**:
   - Check if ports 6513, 5173, 6379, 5433 are available
   - Kill any conflicting processes

4. **Virtual environment not activated**:
   ```bash
   source .venv/bin/activate
   ```

5. **Test dependencies missing**:
   ```bash
   pip install requests
   ```

### Debug Mode:

To see more detailed output, you can modify the logging level in `test_integration.py`:

```python
logging.basicConfig(
    level=logging.DEBUG,  # Change from INFO to DEBUG
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
```

## Test Architecture

The test uses two main classes:

- **ServiceManager**: Handles starting/stopping services (Celery, Flask, React)
- **ToxindexClient**: Provides API client for interacting with the webapp

The test is designed to be:
- **Isolated**: Starts its own service instances
- **Self-cleaning**: Stops all processes when done
- **Comprehensive**: Tests the entire stack end-to-end
- **Informative**: Provides clear feedback on what's working/failing 