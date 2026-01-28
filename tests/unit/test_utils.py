"""Unit tests for utils module."""

import pytest
from pathlib import Path
from mutation_scan.utils import ConfigParser, FileHandler
import tempfile

class TestConfigParser:
    """Test suite for ConfigParser class."""

    def test_get_simple_key(self, tmp_path):
        """Test retrieving simple configuration key."""
        config_file = tmp_path / "config.yaml"
        config_file.write_text("test_key: test_value\n")
        
        parser = ConfigParser(config_file)
        assert parser.get("test_key") == "test_value"

    def test_get_missing_key_returns_default(self, tmp_path):
        """Test that missing key returns default value."""
        config_file = tmp_path / "config.yaml"
        config_file.write_text("existing_key: value\n")
        
        parser = ConfigParser(config_file)
        assert parser.get("missing_key", "default") == "default"

    def test_get_section(self, tmp_path):
        """Test retrieving entire configuration section."""
        config_file = tmp_path / "config.yaml"
        config_file.write_text("section:\n  key1: value1\n  key2: value2\n")
        
        parser = ConfigParser(config_file)
        section = parser.get_section("section")
        assert section["key1"] == "value1"
        assert section["key2"] == "value2"


class TestFileHandler:
    """Test suite for FileHandler class."""

    def test_ensure_dir_creates_directory(self, tmp_path):
        """Test that ensure_dir creates directory."""
        new_dir = tmp_path / "new" / "nested" / "dir"
        result = FileHandler.ensure_dir(new_dir)
        
        assert new_dir.exists()
        assert result == new_dir

    def test_safe_write_creates_file(self, tmp_path):
        """Test that safe_write creates file."""
        file_path = tmp_path / "test.txt"
        content = "test content"
        
        success = FileHandler.safe_write(file_path, content)
        
        assert success
        assert file_path.exists()
        assert file_path.read_text() == content

    def test_safe_write_creates_backup(self, tmp_path):
        """Test that safe_write creates backup of existing file."""
        file_path = tmp_path / "test.txt"
        
        # First write
        FileHandler.safe_write(file_path, "original content")
        
        # Second write (should create backup)
        FileHandler.safe_write(file_path, "new content")
        
        backup_path = file_path.with_suffix(".txt.backup")
        assert backup_path.exists()
        assert backup_path.read_text() == "original content"

    def test_safe_read_reads_file(self, tmp_path):
        """Test that safe_read reads file content."""
        file_path = tmp_path / "test.txt"
        content = "test content"
        file_path.write_text(content)
        
        result = FileHandler.safe_read(file_path)
        
        assert result == content

    def test_safe_read_nonexistent_file(self, tmp_path):
        """Test that safe_read returns None for nonexistent file."""
        file_path = tmp_path / "nonexistent.txt"
        
        result = FileHandler.safe_read(file_path)
        
        assert result is None
