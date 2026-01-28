"""
File handling utilities.

Provides common file operations like reading, writing, and validation.
"""

import logging
import shutil
from pathlib import Path
from typing import Optional, Union

logger = logging.getLogger(__name__)


class FileHandler:
    """Handle common file operations."""

    @staticmethod
    def ensure_dir(dirpath: Path) -> Path:
        """
        Ensure directory exists, create if necessary.

        Args:
            dirpath: Path to directory

        Returns:
            Path object of directory
        """
        dirpath = Path(dirpath)
        dirpath.mkdir(parents=True, exist_ok=True)
        logger.debug(f"Ensured directory exists: {dirpath}")
        return dirpath

    @staticmethod
    def safe_write(filepath: Path, content: Union[str, bytes], mode: str = "w") -> bool:
        """
        Safely write content to file with backup.

        Args:
            filepath: Path to file
            content: Content to write
            mode: Write mode ('w' for text, 'wb' for binary)

        Returns:
            True if successful, False otherwise
        """
        try:
            filepath = Path(filepath)
            FileHandler.ensure_dir(filepath.parent)

            # Backup existing file
            if filepath.exists():
                backup_path = filepath.with_suffix(filepath.suffix + ".backup")
                shutil.copy2(filepath, backup_path)
                logger.debug(f"Created backup: {backup_path}")

            with open(filepath, mode) as f:
                f.write(content)

            logger.info(f"Successfully wrote to {filepath}")
            return True
        except Exception as e:
            logger.error(f"Error writing to {filepath}: {e}")
            return False

    @staticmethod
    def safe_read(filepath: Path, mode: str = "r") -> Optional[Union[str, bytes]]:
        """
        Safely read file content.

        Args:
            filepath: Path to file
            mode: Read mode ('r' for text, 'rb' for binary)

        Returns:
            File content or None if read fails
        """
        try:
            with open(filepath, mode) as f:
                content = f.read()
            logger.debug(f"Successfully read from {filepath}")
            return content
        except FileNotFoundError:
            logger.error(f"File not found: {filepath}")
            return None
        except Exception as e:
            logger.error(f"Error reading {filepath}: {e}")
            return None
