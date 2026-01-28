"""
Configuration file parser module.

Handles parsing and validation of YAML configuration files.
"""

import logging
from pathlib import Path
from typing import Any, Dict, Optional

import yaml

logger = logging.getLogger(__name__)


class ConfigParser:
    """Parse and manage YAML configuration files."""

    def __init__(self, config_path: Path):
        """
        Initialize configuration parser.

        Args:
            config_path: Path to YAML config file
        """
        self.config_path = Path(config_path)
        self.config = self._load_config()
        logger.info(f"Loaded configuration from {config_path}")

    def _load_config(self) -> Dict[str, Any]:
        """
        Load YAML configuration file.

        Returns:
            Configuration dictionary
        """
        try:
            with open(self.config_path, "r") as f:
                config = yaml.safe_load(f)
            return config or {}
        except FileNotFoundError:
            logger.error(f"Configuration file not found: {self.config_path}")
            raise
        except yaml.YAMLError as e:
            logger.error(f"Error parsing YAML config: {e}")
            raise

    def get(self, key: str, default: Any = None) -> Any:
        """
        Get configuration value by key.

        Args:
            key: Configuration key (supports dot notation for nested keys)
            default: Default value if key not found

        Returns:
            Configuration value or default
        """
        keys = key.split(".")
        value = self.config

        for k in keys:
            if isinstance(value, dict):
                value = value.get(k)
            else:
                return default

        return value if value is not None else default

    def get_section(self, section: str) -> Dict[str, Any]:
        """
        Get entire configuration section.

        Args:
            section: Section name

        Returns:
            Section dictionary
        """
        return self.config.get(section, {})

    def validate_required_keys(self, required_keys: list) -> bool:
        """
        Validate that required keys exist in config.

        Args:
            required_keys: List of required key paths

        Returns:
            True if all keys present, raises ValueError otherwise
        """
        missing_keys = [k for k in required_keys if self.get(k) is None]
        if missing_keys:
            raise ValueError(f"Missing required configuration keys: {missing_keys}")
        return True
