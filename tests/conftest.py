"""
conftest.py
Shared pytest configuration and fixtures for the INS-Data test suite.
"""

import concurrent.futures
import pytest


# ---------------------------------------------------------------------------
# Guard against flaky tests that depend on ThreadPoolExecutor ordering.
#
# concurrent.futures.as_completed() returns futures in *completion* order,
# which is non-deterministic.  Tests that use list-based mock side_effects
# with concurrent workers will silently pass when threads happen to run in
# submission order and fail otherwise — classic intermittent flake.
#
# Note: production code is safe because each worker returns its identity
# alongside its result as a tuple (e.g. (pmid, sra_ids)), so completion
# order doesn't affect which key maps to which value.  The vulnerability
# is only in tests that use list-based mock side_effects, where call order
# determines which input gets which return value.
#
# Fix: monkeypatch as_completed to reverse the order during every test.
# Any order-dependent test will then fail consistently and immediately
# instead of only on unlucky CI runs.
# ---------------------------------------------------------------------------

_real_as_completed = concurrent.futures.as_completed


def _reversed_as_completed(fs, timeout=None):
    """Yield completed futures in reversed order to surface order-dependent bugs."""
    return reversed(list(_real_as_completed(fs, timeout=timeout)))


@pytest.fixture(autouse=True)
def _shuffle_as_completed(monkeypatch):
    """Automatically reverse as_completed order for every test."""
    monkeypatch.setattr(concurrent.futures, "as_completed", _reversed_as_completed)


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
