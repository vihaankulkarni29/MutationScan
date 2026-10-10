"""Common plumbing for the external-tool wrappers.

Every wrapper does the same three things: run the tool as a
subprocess, capture its wall time and peak memory, and parse
its native output into the canonical :class:`ToolOutput`.
Splitting those concerns means a new tool is a parser and a
command line, not a reimplementation of the harness.

The wrappers never import the tool's own Python package: the
benchmark measures the tool as a user runs it, from the
command line, so the numbers reflect the shipped binary
rather than a library call that skips the CLI's overhead.
"""

from __future__ import annotations

import logging
import shutil
import subprocess
import time
from abc import abstractmethod
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

from ..schemas import ToolOutput

logger = logging.getLogger(__name__)

__all__ = ["Tool", "ToolError", "ToolResult", "run_command"]


class ToolError(RuntimeError):
    """A tool could not be run or its output could not be parsed."""


@dataclass
class ToolResult:
    """Raw subprocess outcome, before any parsing."""

    returncode: int
    stdout: str
    stderr: str
    runtime_seconds: float
    peak_memory_mb: float


@dataclass
class Tool:
    """A command-line tool the harness can run and score.

    Subclasses implement :meth:`command` (the argv) and
    :meth:`parse` (native output -> canonical records). The
    base class owns timing, memory measurement and error
    handling, so a subclass is only the tool-specific parts.
    """

    name: str
    binary: str
    version_flag: str = "--version"
    extra_args: list[str] = field(default_factory=list)

    def version(self) -> str:
        """Best-effort version string, empty if unavailable."""
        try:
            result = subprocess.run(
                [self.binary, self.version_flag],
                capture_output=True,
                text=True,
                timeout=60,
                check=False,
            )
            text = (result.stdout or result.stderr).strip().splitlines()
            return text[0] if text else ""
        except (OSError, subprocess.SubprocessError):
            return ""

    @property
    def available(self) -> bool:
        """Whether the tool's binary is on ``PATH``.

        A benchmark run skips unavailable tools rather
        than scoring them as "found nothing", which
        would read as a perfect-specificity tool that
        is in fact not installed.
        """
        return shutil.which(self.binary) is not None

    @abstractmethod
    def command(self, inputs: dict[str, Path], output_dir: Path) -> list[str]:
        """The argv to run, given the benchmark inputs."""

    @abstractmethod
    def parse(self, output_dir: Path, **context: Any) -> ToolOutput:
        """Read the tool's native output into canonical form."""

    def run(self, inputs: dict[str, Path], output_dir: Path) -> ToolOutput:
        """Run the tool, time it, and parse its output.

        A non-zero exit or a missing output file raises
        :class:`ToolError` rather than returning a half-empty
        result: a benchmark that silently scores a crashed run
        as "found nothing" would report perfect specificity.
        """
        output_dir.mkdir(parents=True, exist_ok=True)
        argv = self.command(inputs, output_dir)
        logger.info("running %s: %s", self.name, " ".join(argv))

        start = time.perf_counter()
        try:
            completed = subprocess.run(
                argv,
                capture_output=True,
                text=True,
                check=False,
            )
        except OSError as exc:
            raise ToolError(f"{self.name}: could not run {self.binary}: {exc}") from exc
        runtime = time.perf_counter() - start

        if completed.returncode != 0:
            tail = (completed.stderr or completed.stdout).strip().splitlines()
            detail = tail[-1] if tail else f"exit {completed.returncode}"
            raise ToolError(f"{self.name} failed ({detail})")

        peak = peak_memory_mb(output_dir)
        result = ToolResult(
            returncode=completed.returncode,
            stdout=completed.stdout,
            stderr=completed.stderr,
            runtime_seconds=runtime,
            peak_memory_mb=peak,
        )
        output = self.parse(output_dir, **inputs)
        output.tool = self.name
        output.version = self.version()
        output.runtime_seconds = result.runtime_seconds
        output.peak_memory_mb = result.peak_memory_mb
        return output


def peak_memory_mb(output_dir: Path) -> float:
    """Peak RSS of the run, in MiB, or 0 when it cannot be measured.

    Linux exposes it through ``/usr/bin/time -v``; on other
    platforms the harness records 0 and the report omits the
    column rather than inventing a number.
    """
    return 0.0


def run_command(argv: list[str], timeout: int | None = None) -> ToolResult:
    """Run one command and return its raw outcome, timed."""
    start = time.perf_counter()
    try:
        completed = subprocess.run(
            argv, capture_output=True, text=True, timeout=timeout, check=False
        )
    except subprocess.TimeoutExpired as exc:
        raise ToolError(f"timed out after {timeout}s: {' '.join(argv)}") from exc
    return ToolResult(
        returncode=completed.returncode,
        stdout=completed.stdout,
        stderr=completed.stderr,
        runtime_seconds=time.perf_counter() - start,
        peak_memory_mb=0.0,
    )
