"""
Custom metrics exporter for Google Cloud Monitoring.
Sends application-specific metrics to Cloud Monitoring for dashboard visualization.
"""

import os
import time
import logging
from datetime import datetime
from typing import Dict, Any, Optional

# Google Cloud Monitoring imports
from google.cloud import monitoring_v3
from google.api_core import exceptions

logger = logging.getLogger(__name__)

class MetricsExporter:
    """Exports custom metrics to Google Cloud Monitoring."""
    
    def __init__(self, project_id: Optional[str] = None):
        """
        Initialize the metrics exporter.
        
        Args:
            project_id: Google Cloud project ID. If None, will be read from environment.
        """
        self.project_id = project_id or os.environ.get('GOOGLE_CLOUD_PROJECT')
        if not self.project_id:
            raise ValueError("GOOGLE_CLOUD_PROJECT environment variable must be set")
        
        self.client = monitoring_v3.MetricServiceClient()
        self.project_name = f"projects/{self.project_id}"
        
        # Metric type prefixes
        self.metric_prefix = "custom.googleapis.com"
        
        # Cache for metric descriptors
        self._metric_descriptors = {}
    
    def create_metric_descriptor(self, metric_type: str, description: str, unit: str = "") -> None:
        """
        Create a metric descriptor in Cloud Monitoring.
        
        Args:
            metric_type: The metric type (e.g., 'health_check_status')
            description: Description of the metric
            unit: Unit of measurement
        """
        full_metric_type = f"{self.metric_prefix}/{metric_type}"
        
        if full_metric_type in self._metric_descriptors:
            return  # Already created
        
        descriptor = monitoring_v3.MetricDescriptor(
            type=full_metric_type,
            description=description,
            unit=unit,
            metric_kind=monitoring_v3.MetricDescriptor.MetricKind.GAUGE,
            value_type=monitoring_v3.MetricDescriptor.ValueType.DOUBLE,
        )
        
        try:
            self.client.create_metric_descriptor(
                name=self.project_name,
                metric_descriptor=descriptor
            )
            self._metric_descriptors[full_metric_type] = descriptor
            logger.info(f"Created metric descriptor: {full_metric_type}")
        except exceptions.AlreadyExists:
            logger.debug(f"Metric descriptor already exists: {full_metric_type}")
        except Exception as e:
            logger.error(f"Failed to create metric descriptor {full_metric_type}: {e}")
    
    def write_metric(self, metric_type: str, value: float, labels: Dict[str, str] = None) -> None:
        """
        Write a metric to Cloud Monitoring.
        
        Args:
            metric_type: The metric type (e.g., 'health_check_status')
            value: The metric value
            labels: Optional labels for the metric
        """
        try:
            full_metric_type = f"{self.metric_prefix}/{metric_type}"
            
            # Create metric descriptor if it doesn't exist
            if full_metric_type not in self._metric_descriptors:
                self.create_metric_descriptor(metric_type, f"Custom metric: {metric_type}")
            
            # Create the time series
            series = monitoring_v3.TimeSeries()
            series.metric.type = full_metric_type
            
            # Add labels
            if labels:
                series.metric.labels.update(labels)
            
            # Set the resource
            series.resource.type = "k8s_pod"
            series.resource.labels["project_id"] = self.project_id
            series.resource.labels["location"] = os.environ.get("CLOUDSDK_COMPUTE_REGION", "us-central1")
            series.resource.labels["cluster_name"] = os.environ.get("CLOUDSDK_CONTAINER_CLUSTER", "toxindex-cluster")
            series.resource.labels["namespace_name"] = "toxindex-app"
            series.resource.labels["pod_name"] = os.environ.get("HOSTNAME", "unknown")
            
            # Set the point
            point = monitoring_v3.Point()
            point.value.double_value = value
            point.interval.end_time.seconds = int(time.time())
            series.points = [point]
            
            # Write the time series
            self.client.create_time_series(
                name=self.project_name,
                time_series=[series]
            )
            
            logger.debug(f"Wrote metric {metric_type}: {value}")
            
        except Exception as e:
            logger.error(f"Failed to write metric {metric_type}: {e}")
    
    def export_health_check_metrics(self, health_results: Dict[str, Any]) -> None:
        """
        Export health check results as metrics.
        
        Args:
            health_results: Results from health checker
        """
        try:
            # Overall health status (1.0 = healthy, 0.0 = unhealthy)
            overall_status = 1.0 if health_results.get('status') == 'healthy' else 0.0
            self.write_metric('health_check_status', overall_status, {'component': 'overall'})
            
            # Individual component health
            checks = health_results.get('checks', {})
            for check_name, check_result in checks.items():
                if check_name == '_metadata':
                    continue
                
                status_value = 1.0 if check_result.get('status') == 'healthy' else 0.0
                self.write_metric('health_check_status', status_value, {'component': check_name})
                
                # Export response time if available
                if 'response_time' in check_result:
                    self.write_metric(f'{check_name}_response_time', check_result['response_time'])
                
                # Export specific metrics for each component
                if check_name == 'database':
                    self.write_metric('database_connections', check_result.get('tables_accessible', 0))
                elif check_name == 'redis':
                    self.write_metric('redis_operations', 1.0)  # Simple connectivity indicator
                elif check_name == 'celery':
                    self.write_metric('celery_worker_count', check_result.get('worker_count', 0))
                elif check_name == 'memory_usage':
                    self.write_metric('memory_usage_percent', check_result.get('memory_percent', 0))
                elif check_name == 'disk_space':
                    self.write_metric('disk_usage_percent', check_result.get('disk_percent', 0))
                    
        except Exception as e:
            logger.error(f"Failed to export health check metrics: {e}")
    
    def export_http_metrics(self, method: str, endpoint: str, status_code: int, duration: float) -> None:
        """
        Export HTTP request metrics.
        
        Args:
            method: HTTP method (GET, POST, etc.)
            endpoint: API endpoint
            status_code: HTTP status code
            duration: Request duration in seconds
        """
        try:
            # Request count
            self.write_metric('http_requests_total', 1.0, {
                'method': method,
                'endpoint': endpoint,
                'status_code': str(status_code)
            })
            
            # Request duration
            self.write_metric('http_request_duration_seconds', duration, {
                'method': method,
                'endpoint': endpoint
            })
            
            # Error rate indicator
            is_error = 1.0 if status_code >= 400 else 0.0
            self.write_metric('http_errors_total', is_error, {
                'method': method,
                'endpoint': endpoint,
                'status_code': str(status_code)
            })
            
        except Exception as e:
            logger.error(f"Failed to export HTTP metrics: {e}")
    
    def export_business_metrics(self, metrics: Dict[str, float]) -> None:
        """
        Export business-specific metrics.
        
        Args:
            metrics: Dictionary of business metrics
        """
        try:
            for metric_name, value in metrics.items():
                self.write_metric(f'business_{metric_name}', value)
                
        except Exception as e:
            logger.error(f"Failed to export business metrics: {e}")
    
    def export_custom_metric(self, metric_name: str, value: float, labels: Dict[str, str] = None) -> None:
        """
        Export a custom metric.
        
        Args:
            metric_name: Name of the metric
            value: Metric value
            labels: Optional labels
        """
        try:
            self.write_metric(metric_name, value, labels)
        except Exception as e:
            logger.error(f"Failed to export custom metric {metric_name}: {e}")


# Global metrics exporter instance
_metrics_exporter: Optional[MetricsExporter] = None

def get_metrics_exporter() -> MetricsExporter:
    """Get the global metrics exporter instance."""
    global _metrics_exporter
    if _metrics_exporter is None:
        _metrics_exporter = MetricsExporter()
    return _metrics_exporter

def export_health_metrics(health_results: Dict[str, Any]) -> None:
    """Export health check results as metrics."""
    try:
        exporter = get_metrics_exporter()
        exporter.export_health_check_metrics(health_results)
    except Exception as e:
        logger.error(f"Failed to export health metrics: {e}")

def export_http_metrics(method: str, endpoint: str, status_code: int, duration: float) -> None:
    """Export HTTP request metrics."""
    try:
        exporter = get_metrics_exporter()
        exporter.export_http_metrics(method, endpoint, status_code, duration)
    except Exception as e:
        logger.error(f"Failed to export HTTP metrics: {e}")

def export_custom_metric(metric_name: str, value: float, labels: Dict[str, str] = None) -> None:
    """Export a custom metric."""
    try:
        exporter = get_metrics_exporter()
        exporter.export_custom_metric(metric_name, value, labels)
    except Exception as e:
        logger.error(f"Failed to export custom metric: {e}") 