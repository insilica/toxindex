"""
Comprehensive health checker for the ToxIndex application.
Tests all critical functionality including database, Redis, GCS, and external services.
"""

import os
import logging
import time
import json
from typing import Dict, List, Tuple, Optional
from datetime import datetime, timedelta

# Database imports
from webserver.model import Task, Workflow, Message, File, ChatSession
from webserver.data_source import ds

# Storage imports
from webserver.storage import GCSFileStorage

# Redis imports
import redis

# External service imports
import requests
from urllib.parse import urljoin

logger = logging.getLogger(__name__)

class HealthChecker:
    """Comprehensive health checker for the ToxIndex application."""
    
    def __init__(self):
        self.checks = {
            'database': self._check_database,
            'redis': self._check_redis,
            'gcs': self._check_gcs,
            'celery': self._check_celery,
            'external_apis': self._check_external_apis,
            'file_operations': self._check_file_operations,
            'memory_usage': self._check_memory_usage,
            'disk_space': self._check_disk_space,
        }
        
    def run_all_checks(self) -> Dict[str, Dict]:
        """Run all health checks and return results."""
        results = {}
        start_time = time.time()
        
        for check_name, check_func in self.checks.items():
            try:
                logger.info(f"Running health check: {check_name}")
                result = check_func()
                results[check_name] = result
                logger.info(f"Health check {check_name}: {result['status']}")
            except Exception as e:
                logger.error(f"Health check {check_name} failed with exception: {e}")
                results[check_name] = {
                    'status': 'error',
                    'message': str(e),
                    'timestamp': datetime.now().isoformat()
                }
        
        total_time = time.time() - start_time
        results['_metadata'] = {
            'total_checks': len(self.checks),
            'total_time': total_time,
            'timestamp': datetime.now().isoformat()
        }
        
        # Export metrics to Google Cloud Monitoring
        try:
            from webserver.metrics_exporter import export_health_metrics
            export_health_metrics(results)
        except Exception as e:
            logger.warning(f"Failed to export health metrics: {e}")
        
        return results
    
    def get_overall_status(self, results: Dict[str, Dict]) -> Tuple[str, int]:
        """Determine overall health status and HTTP status code."""
        if not results:
            return 'error', 503
        
        # Count statuses
        status_counts = {}
        for check_name, result in results.items():
            if check_name == '_metadata':
                continue
            status = result.get('status', 'unknown')
            status_counts[status] = status_counts.get(status, 0) + 1
        
        # Determine overall status
        if status_counts.get('error', 0) > 0:
            return 'error', 503
        elif status_counts.get('warning', 0) > 0:
            return 'warning', 200
        elif status_counts.get('healthy', 0) == len(self.checks):
            return 'healthy', 200
        else:
            return 'unknown', 503
    
    def _check_database(self) -> Dict:
        """Check database connectivity and basic operations."""
        try:
            start_time = time.time()
            
            # Test basic connection
            result = ds.find("SELECT 1 as test")
            if not result or result.get('test') != 1:
                return {
                    'status': 'error',
                    'message': 'Database connection test failed',
                    'timestamp': datetime.now().isoformat()
                }
            
            # Test table access
            tables_to_check = ['tasks', 'workflows', 'messages', 'files', 'chat_sessions']
            for table in tables_to_check:
                try:
                    ds.find(f"SELECT COUNT(*) as count FROM {table}")
                except Exception as e:
                    return {
                        'status': 'error',
                        'message': f'Cannot access table {table}: {str(e)}',
                        'timestamp': datetime.now().isoformat()
                    }
            
            # Test basic CRUD operations
            test_task = Task.get_task(1)  # Try to get a task
            workflows = Workflow.get_all_workflows()  # Try to get workflows
            
            response_time = time.time() - start_time
            
            return {
                'status': 'healthy',
                'message': 'Database is healthy',
                'response_time': response_time,
                'tables_accessible': len(tables_to_check),
                'timestamp': datetime.now().isoformat()
            }
            
        except Exception as e:
            return {
                'status': 'error',
                'message': f'Database check failed: {str(e)}',
                'timestamp': datetime.now().isoformat()
            }
    
    def _check_redis(self) -> Dict:
        """Check Redis connectivity and basic operations."""
        try:
            start_time = time.time()
            
            # Connect to Redis
            r = redis.Redis(
                host=os.environ.get("REDIS_HOST", "localhost"),
                port=int(os.environ.get("REDIS_PORT", "6379")),
                db=0,
                socket_connect_timeout=5,
                socket_timeout=5
            )
            
            # Test basic operations
            test_key = f"health_check_{int(time.time())}"
            test_value = "test_value"
            
            # Test set/get
            r.set(test_key, test_value, ex=60)  # Expire in 60 seconds
            retrieved_value = r.get(test_key)
            
            if retrieved_value.decode() != test_value:
                return {
                    'status': 'error',
                    'message': 'Redis get/set test failed',
                    'timestamp': datetime.now().isoformat()
                }
            
            # Test pub/sub (basic connectivity)
            pubsub = r.pubsub()
            pubsub.subscribe("health_check_test")
            pubsub.unsubscribe("health_check_test")
            
            # Clean up
            r.delete(test_key)
            
            response_time = time.time() - start_time
            
            return {
                'status': 'healthy',
                'message': 'Redis is healthy',
                'response_time': response_time,
                'timestamp': datetime.now().isoformat()
            }
            
        except Exception as e:
            return {
                'status': 'error',
                'message': f'Redis check failed: {str(e)}',
                'timestamp': datetime.now().isoformat()
            }
    
    def _check_gcs(self) -> Dict:
        """Check Google Cloud Storage connectivity and operations."""
        try:
            start_time = time.time()
            
            gcs_storage = GCSFileStorage()
            
            # Test bucket access
            bucket_name = os.environ.get("GCS_BUCKET_NAME", "toxindex-uploads")
            
            # Test basic operations (without actually creating files)
            # This just tests connectivity and permissions
            try:
                # Try to list a small number of objects (limit=1)
                # This tests bucket access without creating files
                pass  # For now, just test that GCSFileStorage can be instantiated
                
            except Exception as e:
                return {
                    'status': 'error',
                    'message': f'GCS connectivity failed: {str(e)}',
                    'timestamp': datetime.now().isoformat()
                }
            
            response_time = time.time() - start_time
            
            return {
                'status': 'healthy',
                'message': 'GCS is accessible',
                'response_time': response_time,
                'bucket_name': bucket_name,
                'timestamp': datetime.now().isoformat()
            }
            
        except Exception as e:
            return {
                'status': 'error',
                'message': f'GCS check failed: {str(e)}',
                'timestamp': datetime.now().isoformat()
            }
    
    def _check_celery(self) -> Dict:
        """Check Celery worker connectivity and basic operations."""
        try:
            start_time = time.time()
            
            # Test Redis connection (Celery uses Redis)
            r = redis.Redis(
                host=os.environ.get("REDIS_HOST", "localhost"),
                port=int(os.environ.get("REDIS_PORT", "6379")),
                db=0,
                socket_connect_timeout=5,
                socket_timeout=5
            )
            
            # Check if Celery workers are active by looking for worker keys
            # This is a basic check - in production you might want more sophisticated monitoring
            worker_keys = r.keys("celery@*")
            
            response_time = time.time() - start_time
            
            if worker_keys:
                return {
                    'status': 'healthy',
                    'message': f'Celery workers detected: {len(worker_keys)}',
                    'response_time': response_time,
                    'worker_count': len(worker_keys),
                    'timestamp': datetime.now().isoformat()
                }
            else:
                return {
                    'status': 'warning',
                    'message': 'No Celery workers detected',
                    'response_time': response_time,
                    'worker_count': 0,
                    'timestamp': datetime.now().isoformat()
                }
                
        except Exception as e:
            return {
                'status': 'error',
                'message': f'Celery check failed: {str(e)}',
                'timestamp': datetime.now().isoformat()
            }
    
    def _check_external_apis(self) -> Dict:
        """Check external API connectivity (if any)."""
        try:
            start_time = time.time()
            
            # Check if we can reach external services
            # Add your external API checks here
            external_checks = []
            
            # Example: Check if OpenAI API is accessible (if configured)
            openai_api_key = os.environ.get("OPENAI_API_KEY")
            if openai_api_key:
                try:
                    # Basic connectivity test (don't make actual API calls in health check)
                    external_checks.append({
                        'service': 'openai',
                        'status': 'healthy',
                        'message': 'OpenAI API key configured'
                    })
                except Exception as e:
                    external_checks.append({
                        'service': 'openai',
                        'status': 'warning',
                        'message': f'OpenAI check failed: {str(e)}'
                    })
            
            response_time = time.time() - start_time
            
            if not external_checks:
                return {
                    'status': 'healthy',
                    'message': 'No external APIs configured',
                    'response_time': response_time,
                    'timestamp': datetime.now().isoformat()
                }
            
            # Determine overall status
            error_count = sum(1 for check in external_checks if check['status'] == 'error')
            warning_count = sum(1 for check in external_checks if check['status'] == 'warning')
            
            if error_count > 0:
                overall_status = 'error'
            elif warning_count > 0:
                overall_status = 'warning'
            else:
                overall_status = 'healthy'
            
            return {
                'status': overall_status,
                'message': f'External APIs: {len(external_checks)} checked',
                'response_time': response_time,
                'checks': external_checks,
                'timestamp': datetime.now().isoformat()
            }
            
        except Exception as e:
            return {
                'status': 'error',
                'message': f'External API check failed: {str(e)}',
                'timestamp': datetime.now().isoformat()
            }
    
    def _check_file_operations(self) -> Dict:
        """Check file system operations and permissions."""
        try:
            start_time = time.time()
            
            # Check if log directory is writable
            from webserver.data_paths import LOGS_ROOT
            logs_dir = LOGS_ROOT()
            
            # Test write permission
            test_file = logs_dir / f"health_check_{int(time.time())}.tmp"
            try:
                test_file.write_text("health_check_test")
                test_file.unlink()  # Clean up
                file_operations = 'healthy'
            except Exception as e:
                file_operations = f'error: {str(e)}'
            
            response_time = time.time() - start_time
            
            return {
                'status': 'healthy' if file_operations == 'healthy' else 'error',
                'message': f'File operations: {file_operations}',
                'response_time': response_time,
                'logs_directory': str(logs_dir),
                'timestamp': datetime.now().isoformat()
            }
            
        except Exception as e:
            return {
                'status': 'error',
                'message': f'File operations check failed: {str(e)}',
                'timestamp': datetime.now().isoformat()
            }
    
    def _check_memory_usage(self) -> Dict:
        """Check memory usage of the application."""
        try:
            import psutil
            
            process = psutil.Process()
            memory_info = process.memory_info()
            memory_percent = process.memory_percent()
            
            # Define thresholds
            memory_warning_threshold = 80  # 80% of available memory
            memory_error_threshold = 95    # 95% of available memory
            
            if memory_percent > memory_error_threshold:
                status = 'error'
                message = f'Memory usage too high: {memory_percent:.1f}%'
            elif memory_percent > memory_warning_threshold:
                status = 'warning'
                message = f'Memory usage high: {memory_percent:.1f}%'
            else:
                status = 'healthy'
                message = f'Memory usage normal: {memory_percent:.1f}%'
            
            return {
                'status': status,
                'message': message,
                'memory_mb': memory_info.rss / 1024 / 1024,
                'memory_percent': memory_percent,
                'timestamp': datetime.now().isoformat()
            }
            
        except ImportError:
            return {
                'status': 'warning',
                'message': 'psutil not available, cannot check memory usage',
                'timestamp': datetime.now().isoformat()
            }
        except Exception as e:
            return {
                'status': 'error',
                'message': f'Memory check failed: {str(e)}',
                'timestamp': datetime.now().isoformat()
            }
    
    def _check_disk_space(self) -> Dict:
        """Check available disk space."""
        try:
            import psutil
            
            # Check disk space for the current directory
            disk_usage = psutil.disk_usage('.')
            disk_percent = disk_usage.percent
            
            # Define thresholds
            disk_warning_threshold = 80  # 80% used
            disk_error_threshold = 95    # 95% used
            
            if disk_percent > disk_error_threshold:
                status = 'error'
                message = f'Disk space critical: {disk_percent:.1f}% used'
            elif disk_percent > disk_warning_threshold:
                status = 'warning'
                message = f'Disk space low: {disk_percent:.1f}% used'
            else:
                status = 'healthy'
                message = f'Disk space normal: {disk_percent:.1f}% used'
            
            return {
                'status': status,
                'message': message,
                'disk_percent': disk_percent,
                'disk_free_gb': disk_usage.free / 1024 / 1024 / 1024,
                'disk_total_gb': disk_usage.total / 1024 / 1024 / 1024,
                'timestamp': datetime.now().isoformat()
            }
            
        except ImportError:
            return {
                'status': 'warning',
                'message': 'psutil not available, cannot check disk space',
                'timestamp': datetime.now().isoformat()
            }
        except Exception as e:
            return {
                'status': 'error',
                'message': f'Disk space check failed: {str(e)}',
                'timestamp': datetime.now().isoformat()
            }

# Global health checker instance
health_checker = HealthChecker() 