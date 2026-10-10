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
import threading
import time
from abc import abstractmethod
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

from ..schemas import ToolOutput

logger = logging.getLogger(__name__)

__all__ = ["Tool", "ToolError", "ToolResult", "run_command"]

#: Polling interval for the memory sampler, in seconds. Fast
#: enough to catch a spike, coarse enough to not tax the host.
_MEM_POLL_INTERVAL = 0.05

#: Bytes per MiB.
_MIB = 1024 * 1024


def _psutil() -> Any:
    """Import psutil lazily; the benchmark is usable without it."""
    try:
        import psutil

        return psutil
    except ImportError:
        return None


class _MemorySampler:
    """Samples the RSS of a process tree on a background thread.

    The sampler must be started *before* the parent waits on
    the child: once ``communicate()`` returns the process is
    gone and there is nothing left to measure. ``stop()``
    returns the peak summed RSS, in MiB, of the process and
    every descendant -- the honest footprint of a tool that
    shells out to ``tblastn``.
    """

    def __init__(self, popen: subprocess.Popen) -> None:
        self._popen = popen
        self._proc: Any = None
        self._peak = 0.0
        self._stop = threading.Event()
        self._thread: threading.Thread | None = None

    def start(self) -> None:
        psutil = _psutil()
        if psutil is None:
            return
        try:
            self._proc = psutil.Process(self._popen.pid)
        except psutil.NoSuchProcess:
            self._proc = None
            return
        self._thread = threading.Thread(target=self._sample, daemon=True)
        self._thread.start()

    def _tree_rss(self) -> float:
        """Summed RSS, in MiB, of the process and its descendants."""
        if self._proc is None:
            return 0.0
        psutil = _psutil()
        total = 0
        try:
            members = [self._proc, *self._proc.children(recursive=True)]
        except psutil.NoSuchProcess:
            return 0.0
        for member in members:
            try:
                total += member.memory_info().rss
            except (psutil.NoSuchProcess, psutil.AccessDenied):
                continue
        return total / _MIB

    def _sample(self) -> None:
        while not self._stop.is_set():
            value = self._tree_rss()
            if value > self._peak:
                self._peak = value
            self._stop.wait(_MEM_POLL_INTERVAL)

    def stop(self) -> float:
        self._stop.set()
        if self._thread is not None:
            self._thread.join(timeout=1.0)
        return self._peak


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
    def command(self, inputs: dict[str, Any], output_dir: Path) -> list[str]:
        """The argv to run, given the benchmark inputs."""

    @abstractmethod
    def parse(self, output_dir: Path, **context: Any) -> ToolOutput:
        """Read the tool's native output into canonical form."""

    def run(self, inputs: dict[str, Any], output_dir: Path) -> ToolOutput:
        """Run the tool, time it, and parse its output.

        A non-zero exit or a missing output file raises
        :class:`ToolError` rather than returning a half-empty
        result: a benchmark that silently scores a crashed run
        as "found nothing" would report perfect specificity.

        Peak RSS is sampled across the whole process tree on a
        background thread, so a tool's memory footprint is
        measured while it runs rather than reported as 0.
        """
        output_dir.mkdir(parents=True, exist_ok=True)
        argv = self.command(inputs, output_dir)
        logger.info("running %s: %s", self.name, " ".join(argv))

        start = time.perf_counter()
        proc: subprocess.Popen[bytes] | None = None
        try:
            proc = subprocess.Popen(argv, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            # Sample memory while the tool runs: once
            # communicate() returns the process is gone,
            # so the sampler has to be running first.
            sampler = _MemorySampler(proc)
            sampler.start()
            try:
                stdout, stderr = proc.communicate()
            finally:
                peak = sampler.stop()
            runtime = time.perf_counter() - start
        except OSError as exc:
            if proc is not None:
                proc.kill()
            raise ToolError(f"{self.name}: could not run {self.binary}: {exc}") from exc

        text = stdout.decode("utf-8", errors="replace") if stdout else ""
        errtext = stderr.decode("utf-8", errors="replace") if stderr else ""

        if proc.returncode != 0:
            tail = (errtext or text).strip().splitlines()
            detail = tail[-1] if tail else f"exit {proc.returncode}"
            raise ToolError(f"{self.name} failed ({detail})")

        result = ToolResult(
            returncode=proc.returncode,
            stdout=text,
            stderr=errtext,
            runtime_seconds=runtime,
            peak_memory_mb=peak,
        )
        output = self.parse(output_dir, **inputs)
        output.tool = self.name
        output.version = self.version()
        output.runtime_seconds = result.runtime_seconds
        output.peak_memory_mb = result.peak_memory_mb
        return output


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
