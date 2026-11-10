"""Explicit dependency marker for compiled Bride controls only."""
import importlib.util
import pytest

requires_lightweaver = pytest.mark.skipif(
    importlib.util.find_spec("lightweaver") is None,
    reason="Compiled Bride control requires the optional Lightweaver backend",
)
