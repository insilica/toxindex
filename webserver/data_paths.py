"""
Data path configuration utility for toxindex.

This module provides centralized access to data directory paths
defined in config/data_paths.yaml to avoid hardcoding paths throughout the codebase.
"""

import os
import yaml
from pathlib import Path
from typing import Dict, Any, Optional


class DataPaths:
    """Centralized data path configuration manager."""
    
    def __init__(self, config_file: str = "config/data_paths.yaml"):
        """
        Initialize DataPaths with configuration file.
        
        Args:
            config_file: Path to the YAML configuration file
        """
        self.config_file = config_file
        self.project_root = Path(__file__).resolve().parent.parent
        self.config = self._load_config()
        self._resolve_paths()
    
    def _load_config(self) -> Dict[str, Any]:
        """Load configuration from YAML file."""
        config_path = self.project_root / self.config_file
        if not config_path.exists():
            raise FileNotFoundError(f"Configuration file not found: {config_path}")
        
        with open(config_path, 'r') as f:
            return yaml.safe_load(f)
    
    def _resolve_paths(self):
        """Resolve variable references in path definitions."""
        # Get base values
        data_root = self.config['data_root']
        directories = self.config['directories']
        
        # Resolve paths with variable substitution
        paths = self.config['paths']
        resolved_paths = {}
        
        for key, path_template in paths.items():
            # Simple variable substitution
            resolved_path = path_template
            resolved_path = resolved_path.replace('${data_root}', data_root)
            
            # Handle nested directory references
            for dir_key, dir_value in directories.items():
                resolved_path = resolved_path.replace(f'${{directories.{dir_key}}}', dir_value)
            
            resolved_paths[key] = resolved_path
        
        self.config['resolved_paths'] = resolved_paths
    
    def get_path(self, path_key: str) -> Path:
        """
        Get a resolved path for the given key.
        
        Args:
            path_key: Key from the paths section of the config
            
        Returns:
            Path object relative to project root
        """
        if 'resolved_paths' not in self.config:
            self._resolve_paths()
        
        if path_key not in self.config['resolved_paths']:
            raise KeyError(f"Path key '{path_key}' not found in configuration")
        
        return self.project_root / self.config['resolved_paths'][path_key]
    
    def get_absolute_path(self, path_key: str) -> Path:
        """
        Get an absolute path for the given key.
        
        Args:
            path_key: Key from the paths section of the config
            
        Returns:
            Absolute Path object
        """
        return self.get_path(path_key).resolve()
    
    def ensure_directory(self, path_key: str) -> Path:
        """
        Ensure a directory exists and return its path.
        
        Args:
            path_key: Key from the paths section of the config
            
        Returns:
            Path object for the created/verified directory
        """
        path = self.get_path(path_key)
        path.mkdir(parents=True, exist_ok=True)
        return path
    
    def get_cache_path(self, cache_type: str) -> Path:
        """
        Get a cache subdirectory path.
        
        Args:
            cache_type: Key from cache_subdirs section
            
        Returns:
            Path object for the cache directory
        """
        if cache_type not in self.config['cache_subdirs']:
            raise KeyError(f"Cache type '{cache_type}' not found in configuration")
        
        cache_subdir = self.config['cache_subdirs'][cache_type]
        return self.get_path('cache') / cache_subdir
    
    def get_session_path(self, session_type: str) -> Path:
        """
        Get a session subdirectory path.
        
        Args:
            session_type: Key from session_subdirs section
            
        Returns:
            Path object for the session directory
        """
        if session_type not in self.config['session_subdirs']:
            raise KeyError(f"Session type '{session_type}' not found in configuration")
        
        session_subdir = self.config['session_subdirs'][session_type]
        return self.get_path('sessions') / session_subdir
    
    def list_available_paths(self) -> Dict[str, str]:
        """List all available path keys and their resolved values."""
        if 'resolved_paths' not in self.config:
            self._resolve_paths()
        
        return {
            key: str(self.project_root / path)
            for key, path in self.config['resolved_paths'].items()
        }


# Global instance for easy access
_data_paths: Optional[DataPaths] = None


def get_data_paths() -> DataPaths:
    """Get the global DataPaths instance, creating it if necessary."""
    global _data_paths
    if _data_paths is None:
        _data_paths = DataPaths()
    return _data_paths


# Convenience functions for common operations
def get_path(path_key: str) -> Path:
    """Get a path for the given key."""
    return get_data_paths().get_path(path_key)


def get_absolute_path(path_key: str) -> Path:
    """Get an absolute path for the given key."""
    return get_data_paths().get_absolute_path(path_key)


def ensure_directory(path_key: str) -> Path:
    """Ensure a directory exists and return its path."""
    return get_data_paths().ensure_directory(path_key)


def get_cache_path(cache_type: str) -> Path:
    """Get a cache subdirectory path."""
    return get_data_paths().get_cache_path(cache_type)


def get_session_path(session_type: str) -> Path:
    """Get a session subdirectory path."""
    return get_data_paths().get_session_path(session_type)


# Common path constants for easy access
def OUTPUTS_ROOT() -> Path:
    """Get the outputs directory path."""
    return get_path('outputs')


def TMP_ROOT() -> Path:
    """Get the temporary files directory path."""
    return get_path('tmp')


def CHATS_ROOT() -> Path:
    """Get the chats directory path."""
    return get_path('chats')


def CACHE_ROOT() -> Path:
    """Get the cache directory path."""
    return get_path('cache')


def LOGS_ROOT() -> Path:
    """Get the logs directory path."""
    return get_path('logs')


def UPLOADS_ROOT() -> Path:
    """Get the uploads directory path."""
    return get_path('uploads') 