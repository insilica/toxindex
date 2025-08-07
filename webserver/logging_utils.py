"""
Shared logging utilities for toxindex components.
Provides consistent Cloud Logging integration across all services.
"""

import os
import logging
from datetime import datetime
from pathlib import Path
from typing import Optional


def setup_logging(
    service_name: str,
    log_level: int = logging.INFO,
    logs_root: Optional[Path] = None
) -> None:
    """
    Setup logging with local file, stdout, and Cloud Logging.
    
    Args:
        service_name: Name of the service (e.g., 'webserver', 'celery-worker', 'redis-listener')
        log_level: Logging level (default: INFO)
        logs_root: Optional custom logs directory path
    """
    # Ensure logs directory exists
    if logs_root is None:
        from webserver.data_paths import LOGS_ROOT
        logs_root = LOGS_ROOT()
    
    logs_root.mkdir(parents=True, exist_ok=True)
    log_filename = logs_root / f"{service_name}_{datetime.now().strftime('%Y-%m-%d_%H')}.log"
    
    handlers = [
        logging.FileHandler(log_filename),  # Local file for debugging
        logging.StreamHandler()             # Stdout for container logs
    ]
    
    # Add Cloud Logging handler if running in GKE
    if os.environ.get("KUBERNETES_SERVICE_HOST"):
        try:
            from google.cloud import logging as cloud_logging
            from google.cloud.logging.handlers import CloudLoggingHandler
            
            # Initialize Cloud Logging client
            client = cloud_logging.Client()
            cloud_handler = CloudLoggingHandler(client, name=service_name)
            
            # Set Cloud Logging level to INFO
            cloud_handler.setLevel(logging.INFO)
            
            handlers.append(cloud_handler)
            print(f"✅ Cloud Logging enabled for {service_name} in GKE environment")
        except ImportError:
            print(f"⚠️  google-cloud-logging not available, skipping Cloud Logging for {service_name}")
        except Exception as e:
            print(f"⚠️  Failed to setup Cloud Logging for {service_name}: {e}")
    else:
        print(f"ℹ️  Running {service_name} locally, Cloud Logging disabled")

    logging.basicConfig(
        level=log_level,
        format='%(asctime)s - %(name)s - %(levelname)s - %(filename)s:%(lineno)d - %(message)s',
        handlers=handlers
    )


def get_logger(name: str) -> logging.Logger:
    """
    Get a logger with the specified name.
    
    Args:
        name: Logger name
        
    Returns:
        Configured logger instance
    """
    return logging.getLogger(name)


def log_service_startup(service_name: str, **kwargs):
    """
    Log service startup information.
    
    Args:
        service_name: Name of the service
        **kwargs: Additional key-value pairs to log
    """
    logger = get_logger(service_name)
    logger.info(f"🚀 {service_name} starting up")
    
    # Log environment information
    env_info = {
        "KUBERNETES_SERVICE_HOST": os.environ.get("KUBERNETES_SERVICE_HOST", "Not set"),
        "REDIS_HOST": os.environ.get("REDIS_HOST", "localhost"),
        "REDIS_PORT": os.environ.get("REDIS_PORT", "6379"),
        "PGHOST": os.environ.get("PGHOST", "Not set"),
        "PGPORT": os.environ.get("PGPORT", "5432"),
        "PGDATABASE": os.environ.get("PGDATABASE", "Not set"),
        "PGUSER": os.environ.get("PGUSER", "Not set"),
    }
    
    # Add any additional kwargs
    env_info.update(kwargs)
    
    for key, value in env_info.items():
        logger.info(f"📋 {key}: {value}")


def log_service_shutdown(service_name: str):
    """
    Log service shutdown information.
    
    Args:
        service_name: Name of the service
    """
    logger = get_logger(service_name)
    logger.info(f"🛑 {service_name} shutting down") 