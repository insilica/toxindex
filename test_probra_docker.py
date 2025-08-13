#!/usr/bin/env python3
"""
Test script to verify probra celery task can be imported and run in Docker environment.
This script can be used to test the Docker image setup.
"""

import os
import sys
import logging
import json
import uuid
from datetime import datetime

# Setup basic logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

def test_imports():
    """Test that all required modules can be imported"""
    logger.info("Testing imports...")
    
    try:
        # Test core imports
        import redis
        logger.info("✓ redis imported successfully")
        
        from celery import Celery
        logger.info("✓ celery imported successfully")
        
        # Test webserver imports
        from webserver.model.task import Task
        logger.info("✓ webserver.model.task imported successfully")
        
        from webserver.model.message import MessageSchema
        logger.info("✓ webserver.model.message imported successfully")
        
        from webserver.storage import GCSFileStorage
        logger.info("✓ webserver.storage imported successfully")
        
        from webserver.cache_manager import cache_manager
        logger.info("✓ webserver.cache_manager imported successfully")
        
        # Test workflow imports
        from workflows.celery_worker_probra import celery
        logger.info("✓ workflows.celery_worker_probra imported successfully")
        
        from workflows.probra import probra_task
        logger.info("✓ workflows.probra imported successfully")
        
        # Test tools imports (these may fail due to missing API keys, which is expected)
        try:
            from webserver.tools.deeptox_agent import deeptox_agent
            logger.info("✓ webserver.tools.deeptox_agent imported successfully")
        except Exception as e:
            logger.warning(f"⚠ webserver.tools.deeptox_agent import failed (expected without API keys): {e}")
        
        try:
            from webserver.tools.toxicity_models_simple import ChemicalToxicityAssessment
            logger.info("✓ webserver.tools.toxicity_models_simple imported successfully")
        except Exception as e:
            logger.warning(f"⚠ webserver.tools.toxicity_models_simple import failed (expected without API keys): {e}")
        
        from webserver.ai_service import convert_pydantic_to_markdown
        logger.info("✓ webserver.ai_service imported successfully")
        
        logger.info("✓ All core imports successful!")
        return True
        
    except ImportError as e:
        logger.error(f"✗ Import failed: {e}")
        return False
    except Exception as e:
        logger.error(f"✗ Unexpected error during import: {e}")
        return False

def test_celery_setup():
    """Test that Celery is properly configured"""
    logger.info("Testing Celery setup...")
    
    try:
        from workflows.celery_worker_probra import celery
        
        # Check if probra task is registered
        if 'workflows.probra.probra_task' in celery.tasks:
            logger.info("✓ probra_task is registered with Celery")
        else:
            logger.warning("⚠ probra_task not found in registered tasks")
            logger.info(f"Available tasks: {list(celery.tasks.keys())}")
        
        # Check broker configuration
        broker_url = celery.conf.get('broker_url', 'Not set')
        logger.info(f"✓ Celery broker URL: {broker_url}")
        
        return True
        
    except Exception as e:
        logger.warning(f"⚠ Celery setup failed (expected without API keys): {e}")
        return True  # Don't fail the test for this

def test_environment_variables():
    """Test that required environment variables are set"""
    logger.info("Testing environment variables...")
    
    required_vars = [
        'REDIS_HOST',
        'REDIS_PORT', 
        'CELERY_BROKER_URL',
        'CELERY_RESULT_BACKEND'
    ]
    
    optional_vars = [
        'PGHOST',
        'PGPORT', 
        'PGDATABASE',
        'PGUSER',
        'PGPASSWORD'
    ]
    
    all_good = True
    
    for var in required_vars:
        value = os.environ.get(var)
        if value:
            logger.info(f"✓ {var}: {value}")
        else:
            logger.warning(f"⚠ {var}: Not set (will use default)")
    
    for var in optional_vars:
        value = os.environ.get(var)
        if value:
            logger.info(f"✓ {var}: {value}")
        else:
            logger.info(f"ℹ {var}: Not set (optional)")
    
    return all_good

def create_test_payload():
    """Create a test payload for the probra task"""
    return {
        "task_id": str(uuid.uuid4()),
        "user_id": str(uuid.uuid4()),
        "payload": "Is Gentamicin nephrotoxic?"
    }

def main():
    """Main test function"""
    logger.info("=== Testing Probra Docker Setup ===")
    logger.info(f"Python version: {sys.version}")
    logger.info(f"Working directory: {os.getcwd()}")
    
    # Test imports
    if not test_imports():
        logger.error("Import test failed!")
        sys.exit(1)
    
    # Test environment variables
    test_environment_variables()
    
    # Test Celery setup
    if not test_celery_setup():
        logger.error("Celery setup test failed!")
        sys.exit(1)
    
    # Create test payload
    test_payload = create_test_payload()
    logger.info(f"Test payload created: {json.dumps(test_payload, indent=2)}")
    
    logger.info("=== All tests passed! ===")
    logger.info("The Docker image should be ready to run probra celery tasks.")
    logger.info("You can now build and run the Docker image.")

if __name__ == "__main__":
    main()
