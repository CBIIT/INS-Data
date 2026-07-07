"""
conftest.py
Shared pytest configuration and fixtures for the INS-Data test suite.
"""

import pytest


def pytest_collection_modifyitems(items):
    """Reorder test collection so live_api tests always run last."""
    live, rest = [], []
    for item in items:
        (live if item.get_closest_marker("live_api") else rest).append(item)
    items[:] = rest + live


_live_banner_shown = False

def pytest_runtest_protocol(item, nextitem):
    """Print a visual separator before the first live_api test."""
    global _live_banner_shown
    if not _live_banner_shown and item.get_closest_marker("live_api"):
        _live_banner_shown = True
        tw = item.config.get_terminal_writer()
        tw.line()
        tw.sep("=", "LIVE API TESTS — dependent on external services", bold=True)
