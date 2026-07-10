"""
conftest.py
Shared pytest configuration and fixtures for the INS-Data test suite.
"""

import concurrent.futures
import pytest


# ---------------------------------------------------------------------------
# Guard against flaky tests that depend on the order in which results from
# ThreadPoolExecutor futures are consumed.
#
# concurrent.futures.as_completed() yields futures in *completion* order,
# which is non-deterministic.  Downstream code that assembles results by
# consumption order (rather than by an identity carried in each result) can
# silently pass when threads happen to complete in submission order and fail
# otherwise — classic intermittent flake.
#
# Note: production code is safe because each worker returns its identity
# alongside its result as a tuple (e.g. (pmid, sra_ids)), so completion
# order doesn't affect which key maps to which value.  The vulnerability
# is only in tests/code paths that rely on the consumption order of
# as_completed to correlate results with inputs.
#
# Fix: monkeypatch as_completed to force a deterministic reversed
# consumption order during every test.  Note that this does NOT influence
# the order in which worker threads invoke the mocked callable (that is
# controlled by the executor's scheduling of submitted tasks); it only
# changes the order in which their completed futures are handed back to
# the caller.  It also changes semantics slightly: because the reversed
# order requires materializing the iterator, callers will block until all
# futures have completed before receiving any result, instead of receiving
# results incrementally as each future completes.
# ---------------------------------------------------------------------------

_real_as_completed = concurrent.futures.as_completed


def _reversed_as_completed(fs, timeout=None):
    """Yield completed futures in reversed order to surface order-dependent bugs."""
    return reversed(list(_real_as_completed(fs, timeout=timeout)))


@pytest.fixture(autouse=True)
def _reverse_as_completed(monkeypatch):
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
