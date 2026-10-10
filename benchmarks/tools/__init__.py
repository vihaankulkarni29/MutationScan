"""External-tool wrappers for the AMR benchmark.

Each module wraps one command-line tool and normalises
its native output into the canonical :mod:`schemas`.
The harness imports the registry below; a new tool is
added by writing a wrapper and registering it here.
"""

from __future__ import annotations

from typing import Any, cast

from .amrfinderplus import AmrfinderPlus
from .base import Tool, ToolError, ToolResult
from .kleborate import Kleborate
from .mutationscan import MutationScanTool
from .mykrobe import Mykrobe
from .staramr import StarAMR

__all__ = [
    "AmrfinderPlus",
    "Kleborate",
    "Mykrobe",
    "MutationScanTool",
    "StarAMR",
    "Tool",
    "ToolError",
    "ToolResult",
]

#: The tools the harness knows how to run, by name.
REGISTRY: dict[str, type[Tool]] = {
    "amrfinderplus": AmrfinderPlus,
    "mykrobe": Mykrobe,
    "kleborate": Kleborate,
    "staramr": StarAMR,
    "mutationscan": MutationScanTool,
}


def build_tool(name: str, **kwargs: Any) -> Tool:
    """Instantiate a tool by its registry name."""
    try:
        cls = REGISTRY[name]
    except KeyError:
        available = ", ".join(sorted(REGISTRY))
        raise ValueError(f"unknown tool '{name}'; available: {available}") from None
    return cast(Tool, cls(**kwargs))
