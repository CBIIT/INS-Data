"""
test_check_url_health.py
2026-10-13 ZD

Pytest test suite for the check_url_health.py module.
"""

import csv
from unittest.mock import MagicMock, patch

import pytest
import requests

from modules.check_url_health import (
    split_url_field,
    read_target_urls,
    collect_urls,
    request_once,
    describe_error,
    check_url,
    normalize_host,
    is_root_path,
    is_benign_redirect,
    rule_connection_error,
    rule_too_many_redirects,
    rule_client_error,
    rule_server_error,
    rule_unexpected_status,
    rule_cross_domain,
    rule_homepage_bounce,
    rule_non_benign_redirect,
    rule_long_redirect_chain,
    rule_slow_response,
    classify_result,
    format_chain,
    build_report_row,
    write_report,
    REPORT_COLUMNS,
)


# ============================================================
# Helpers & Fixtures
# ============================================================

def make_response(status, location=None):
    """Build a fake requests.Response."""
    resp = MagicMock()
    resp.status_code = status
    resp.headers = {'Location': location} if location else {}
    return resp


def make_result(original, final=None, status=200,
                chain=None, response_time=0.1,
                error=None, too_many=False):
    """Build a check_url-style result dict."""
    final = final or original
    if chain is None:
        chain = [(status, original)]
    return {
        'original_url': original,
        'final_url': final,
        'status_code': status,
        'chain': chain,
        'response_time': response_time,
        'error': error,
        'too_many_redirects': too_many,
    }


@pytest.fixture
def mock_config():
    """Patch module config with default settings."""
    with patch("modules.check_url_health.config") as mock:
        mock.URL_HEALTH_TIMEOUT_SECONDS = 5
        mock.URL_HEALTH_MAX_REDIRECT_HOPS = 5
        mock.URL_HEALTH_SLOW_THRESHOLD_FRACTION = 0.5
        mock.URL_HEALTH_YELLOW_REDIRECT_HOPS = 2
        mock.URL_HEALTH_USER_AGENT = "test-agent"
        yield mock


@pytest.fixture
def sample_tsv(tmp_path):
    """TSV with single, semicolon-list, and empty URL values."""
    path = tmp_path / "datasets.tsv"
    path.write_text(
        "id\turl\tlinks\n"
        "1\thttps://a.gov\t https://b.gov ; https://c.gov \n"
        "2\t\t;;\n"
        "3\thttps://a.gov\thttps://d.gov\n",
        encoding="utf-8",
    )
    return str(path)


# ============================================================
# split_url_field
# ============================================================

class TestSplitUrlField:

    def test_single_url(self):
        assert split_url_field("https://a.gov") == ["https://a.gov"]

    def test_semicolon_list_trimmed(self):
        result = split_url_field(" https://a.gov ;https://b.gov ")
        assert result == ["https://a.gov", "https://b.gov"]

    def test_skips_empty_parts(self):
        assert split_url_field("https://a.gov;; ;") == ["https://a.gov"]

    @pytest.mark.parametrize("value", ["", None, "   ", ";"])
    def test_empty_values(self, value):
        assert split_url_field(value) == []


# ============================================================
# read_target_urls / collect_urls
# ============================================================

class TestReadTargetUrls:

    def test_extracts_urls_with_locations(self, sample_tsv):
        pairs = read_target_urls(sample_tsv, ["url", "links"])
        urls = [u for u, _ in pairs]
        assert urls == ["https://a.gov", "https://b.gov",
                        "https://c.gov", "https://a.gov",
                        "https://d.gov"]
        assert pairs[0][1] == "datasets.tsv | url | row 2"
        assert pairs[3][1] == "datasets.tsv | url | row 4"

    def test_missing_file_returns_empty(self, tmp_path, capsys):
        result = read_target_urls(str(tmp_path / "nope.tsv"), ["url"])
        assert result == []
        assert "File not found" in capsys.readouterr().out

    def test_missing_column_warns_and_continues(self, sample_tsv, capsys):
        pairs = read_target_urls(sample_tsv, ["url", "bogus"])
        assert len(pairs) == 2
        assert "bogus" in capsys.readouterr().out


class TestCollectUrls:

    def test_dedupes_and_tracks_occurrences(self, sample_tsv):
        result = collect_urls({sample_tsv: ["url", "links"]})
        assert list(result) == ["https://a.gov", "https://b.gov",
                                "https://c.gov", "https://d.gov"]
        assert len(result["https://a.gov"]) == 2
        assert len(result["https://d.gov"]) == 1


# ============================================================
# request_once / describe_error / check_url
# ============================================================

class TestRequestOnce:

    def test_head_success_skips_get(self):
        session = MagicMock()
        session.head.return_value = make_response(200)
        assert request_once(session, "https://a.gov", 5).status_code == 200
        session.get.assert_not_called()

    @pytest.mark.parametrize("status", [405, 501, 403])
    def test_falls_back_to_get(self, status):
        session = MagicMock()
        session.head.return_value = make_response(status)
        session.get.return_value = make_response(200)
        assert request_once(session, "https://a.gov", 5).status_code == 200
        session.get.assert_called_once()


class TestDescribeError:

    @pytest.mark.parametrize("exc, expected", [
        (requests.exceptions.SSLError("bad cert"), "SSL/TLS error"),
        (requests.exceptions.ConnectTimeout("t"), "timeout"),
        (requests.exceptions.ReadTimeout("t"), "timeout"),
        (requests.exceptions.ConnectionError(
            "NameResolutionError: failed"), "DNS failure"),
        (requests.exceptions.ConnectionError(
            "[Errno 111] Connection refused"), "connection refused"),
        (requests.exceptions.ConnectionError("reset"), "connection error"),
        (requests.exceptions.MissingSchema("no"),
         "invalid or unsupported URL"),
        (requests.exceptions.InvalidSchema("ftp"),
         "invalid or unsupported URL"),
    ])
    def test_maps_exceptions(self, exc, expected):
        assert describe_error(exc) == expected


class TestCheckUrl:

    def test_no_redirect(self):
        session = MagicMock()
        session.head.return_value = make_response(200)
        result = check_url(session, "https://a.gov/x", 5, 5)
        assert result['status_code'] == 200
        assert result['final_url'] == "https://a.gov/x"
        assert result['chain'] == [(200, "https://a.gov/x")]
        assert result['error'] is None

    def test_follows_relative_redirect(self):
        session = MagicMock()
        session.head.side_effect = [
            make_response(301, "/new"),
            make_response(200),
        ]
        result = check_url(session, "https://a.gov/old", 5, 5)
        assert result['final_url'] == "https://a.gov/new"
        assert result['chain'] == [(301, "https://a.gov/old"),
                                   (200, "https://a.gov/new")]

    def test_too_many_redirects(self):
        session = MagicMock()
        session.head.return_value = make_response(302, "/loop")
        result = check_url(session, "https://a.gov/", 5, 2)
        assert result['too_many_redirects'] is True
        assert len(result['chain']) == 3

    def test_captures_connection_error(self):
        session = MagicMock()
        session.head.side_effect = requests.exceptions.ConnectTimeout()
        result = check_url(session, "https://a.gov", 5, 5)
        assert result['error'] == "timeout"
        assert result['status_code'] is None

    def test_unexpected_exception_does_not_raise(self):
        session = MagicMock()
        session.head.side_effect = RuntimeError("boom")
        result = check_url(session, "https://a.gov", 5, 5)
        assert result['error'] == "unexpected error (RuntimeError)"


# ============================================================
# URL comparison helpers
# ============================================================

class TestUrlHelpers:

    def test_normalize_host_strips_www_and_case(self):
        assert normalize_host("https://WWW.Cancer.gov/x") == "cancer.gov"

    @pytest.mark.parametrize("url, expected", [
        ("https://a.gov", True),
        ("https://a.gov/", True),
        ("https://a.gov/?q=1", True),
        ("https://a.gov/page", False),
    ])
    def test_is_root_path(self, url, expected):
        assert is_root_path(url) is expected

    @pytest.mark.parametrize("original, final", [
        ("http://a.gov/x", "https://a.gov/x"),
        ("https://a.gov/x", "https://a.gov/x/"),
        ("https://a.gov/x/", "https://a.gov/x"),
        ("https://a.gov/x", "https://www.a.gov/x"),
        ("https://www.a.gov/x", "https://a.gov/x"),
        ("http://a.gov/x", "https://www.a.gov/x/"),
    ])
    def test_benign(self, original, final):
        assert is_benign_redirect(original, final) is True

    @pytest.mark.parametrize("original, final", [
        ("https://a.gov/x", "http://a.gov/x"),
        ("https://a.gov/x", "https://a.gov/y"),
        ("https://a.gov/x", "https://b.gov/x"),
        ("https://a.gov/x?id=1", "https://a.gov/x?id=2"),
        ("https://a.gov/x", "https://sub.a.gov/x"),
    ])
    def test_not_benign(self, original, final):
        assert is_benign_redirect(original, final) is False


# ============================================================
# Red rules
# ============================================================

class TestRedRules:

    def test_connection_error(self):
        result = make_result("https://a.gov", status=None,
                             chain=[], error="timeout")
        assert rule_connection_error(result) == "timeout"
        assert rule_connection_error(make_result("https://a.gov")) is None

    def test_too_many_redirects(self):
        result = make_result("https://a.gov", status=302, too_many=True)
        assert "too many redirects" in rule_too_many_redirects(result)
        assert rule_too_many_redirects(make_result("https://a.gov")) is None

    @pytest.mark.parametrize("status, text", [
        (404, "404 not found"),
        (403, "403 forbidden"),
        (410, "410 gone"),
    ])
    def test_client_error(self, status, text):
        result = make_result("https://a.gov", status=status)
        assert rule_client_error(result) == text
        assert rule_server_error(result) is None

    @pytest.mark.parametrize("status", [500, 502, 503])
    def test_server_error(self, status):
        result = make_result("https://a.gov", status=status)
        assert rule_server_error(result).startswith(str(status))
        assert rule_client_error(result) is None

    def test_redirect_ending_in_error(self):
        result = make_result(
            "https://a.gov/x", "https://a.gov/y", status=404,
            chain=[(301, "https://a.gov/x"), (404, "https://a.gov/y")])
        assert rule_client_error(result) == "redirect ends in 404 not found"

    def test_none_status_is_not_http_error(self):
        result = make_result("https://a.gov", status=None, chain=[])
        assert rule_client_error(result) is None
        assert rule_server_error(result) is None


# ============================================================
# Yellow rules
# ============================================================

class TestYellowRules:

    def test_unexpected_status(self):
        assert rule_unexpected_status(
            make_result("https://a.gov", status=204)) == \
            "unexpected status 204"
        assert rule_unexpected_status(make_result("https://a.gov")) is None

    def test_cross_domain(self):
        result = make_result("https://a.gov/x", "https://b.org/x")
        assert rule_cross_domain(result) == "cross-domain redirect to b.org"

    def test_cross_domain_ignores_www(self):
        result = make_result("https://a.gov/x", "https://www.a.gov/x")
        assert rule_cross_domain(result) is None

    def test_cross_domain_flags_subdomain_change(self):
        result = make_result("https://dceg.a.gov/x", "https://www.a.gov/x")
        assert rule_cross_domain(result) is not None

    def test_homepage_bounce(self):
        result = make_result("https://a.gov/study/1", "https://a.gov/")
        assert rule_homepage_bounce(result) == "redirect to homepage"

    def test_homepage_bounce_not_flagged_for_root_original(self):
        result = make_result("https://a.gov", "https://www.a.gov/")
        assert rule_homepage_bounce(result) is None

    def test_non_benign_redirect(self):
        chain = [(301, "https://a.gov/x"), (200, "https://a.gov/y")]
        result = make_result("https://a.gov/x", "https://a.gov/y",
                             chain=chain)
        assert rule_non_benign_redirect(result) is not None

    def test_non_benign_ignores_benign(self):
        chain = [(301, "http://a.gov/x"), (200, "https://a.gov/x/")]
        result = make_result("http://a.gov/x", "https://a.gov/x/",
                             chain=chain)
        assert rule_non_benign_redirect(result) is None

    def test_long_redirect_chain(self):
        chain = [(301, "u")] * 3 + [(200, "u")]
        result = make_result("https://a.gov", chain=chain)
        assert rule_long_redirect_chain(result, 2) == \
            "long redirect chain (3 hops)"
        assert rule_long_redirect_chain(result, 3) is None

    def test_slow_response(self):
        assert rule_slow_response(
            make_result("https://a.gov", response_time=3.0), 5, 0.5) == \
            "slow response (3.0s)"
        assert rule_slow_response(
            make_result("https://a.gov", response_time=2.0), 5, 0.5) is None


# ============================================================
# classify_result
# ============================================================

class TestClassifyResult:

    def test_green_no_redirect(self, mock_config):
        assert classify_result(make_result("https://a.gov/x")) == ('G', 'ok')

    def test_green_benign_redirect(self, mock_config):
        chain = [(301, "http://a.gov/x"), (200, "https://www.a.gov/x/")]
        result = make_result("http://a.gov/x", "https://www.a.gov/x/",
                             chain=chain)
        label, _ = classify_result(result)
        assert label == 'G'

    def test_red_beats_yellow(self, mock_config):
        result = make_result("https://a.gov/x", "https://b.org/",
                             status=404, response_time=4.0)
        label, reason = classify_result(result)
        assert label == 'R'
        assert "404" in reason

    def test_yellow_cross_domain_suppresses_generic(self, mock_config):
        chain = [(301, "https://a.gov/x"), (200, "https://b.org/y")]
        result = make_result("https://a.gov/x", "https://b.org/y",
                             chain=chain)
        label, reason = classify_result(result)
        assert label == 'Y'
        assert reason == "cross-domain redirect to b.org"

    def test_yellow_multiple_reasons(self, mock_config):
        chain = [(301, "https://a.gov/x"), (200, "https://a.gov/")]
        result = make_result("https://a.gov/x", "https://a.gov/",
                             chain=chain, response_time=3.0)
        label, reason = classify_result(result)
        assert label == 'Y'
        assert "redirect to homepage" in reason
        assert "slow response" in reason

    def test_uses_config_thresholds(self, mock_config):
        mock_config.URL_HEALTH_SLOW_THRESHOLD_FRACTION = 0.1
        result = make_result("https://a.gov", response_time=1.0)
        assert classify_result(result)[0] == 'Y'


# ============================================================
# Reporting
# ============================================================

class TestReporting:

    def test_format_chain(self):
        chain = [(301, "http://a.gov"), (200, "https://a.gov")]
        assert format_chain(chain) == \
            "301 http://a.gov -> 200 https://a.gov"

    def test_build_report_row(self, mock_config):
        result = make_result("https://a.gov", status=404)
        row = build_report_row(result, ["f | url | row 2",
                                        "f | url | row 5"])
        assert set(row) == set(REPORT_COLUMNS)
        assert row['classification'] == 'R'
        assert row['occurrence_count'] == 2
        assert row['source_locations'] == "f | url | row 2; f | url | row 5"
        assert row['url_changed'] == 'not_changed'

    def test_build_report_row_url_changed(self, mock_config):
        chain = [(301, "http://a.gov"), (200, "https://a.gov")]
        result = make_result("http://a.gov", "https://a.gov", chain=chain)
        row = build_report_row(result, ["loc"])
        assert row['url_changed'] == 'changed'

    def test_build_report_row_blank_status_on_error(self, mock_config):
        result = make_result("https://a.gov", status=None, chain=[],
                             error="DNS failure")
        row = build_report_row(result, ["loc"])
        assert row['final_status_code'] == ''
        assert row['reason'] == "DNS failure"

    def test_write_report_orders_red_first(self, tmp_path):
        rows = [
            {c: '' for c in REPORT_COLUMNS} | {'classification': 'G'},
            {c: '' for c in REPORT_COLUMNS} | {'classification': 'R'},
            {c: '' for c in REPORT_COLUMNS} | {'classification': 'Y'},
        ]
        path = tmp_path / "out" / "report.csv"
        write_report(rows, str(path))

        with open(path, newline='', encoding='utf-8') as f:
            labels = [r['classification'] for r in csv.DictReader(f)]
        assert labels == ['R', 'Y', 'G']
