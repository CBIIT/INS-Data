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
    is_benign_redirect,
    rule_client_error,
    rule_server_error,
    rule_cross_domain,
    rule_homepage_bounce,
    rule_long_redirect_chain,
    rule_slow_response,
    classify_result,
    build_report_row,
    write_report,
    split_reason,
    count_by_file,
    group_flagged_by_file,
    print_summary,
    check_url_health,
    REPORT_COLUMNS,
)


A = "https://a.gov"


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


def make_redirect(original, final, **kwargs):
    """Build a result for a single 301 redirect ending in 200."""
    chain = [(301, original), (200, final)]
    return make_result(original, final, chain=chain, **kwargs)


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

    @pytest.mark.parametrize("value, expected", [
        (" https://a.gov ;https://b.gov ", ["https://a.gov", "https://b.gov"]),
        ("https://a.gov;; ;", ["https://a.gov"]),
        ("   ", []),
        (None, []),
    ])
    def test_splits_trims_and_skips_empties(self, value, expected):
        assert split_url_field(value) == expected


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

    def test_undecodable_file_warns_instead_of_crashing(
            self, tmp_path, capsys):
        path = tmp_path / "bad.tsv"
        path.write_bytes(b"url\n\xff\xfe\xfa\n")
        assert read_target_urls(str(path), ["url"]) == []
        assert "Could not read" in capsys.readouterr().out


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

    @pytest.mark.parametrize("head_status, uses_get", [
        (200, False),
        (301, False),
        (403, True),
        (405, True),
    ])
    def test_get_fallback_only_when_head_rejected(
            self, head_status, uses_get):
        session = MagicMock()
        session.head.return_value = make_response(head_status)
        session.get.return_value = make_response(200)

        request_once(session, A, 5)

        if uses_get:
            session.get.assert_called_once_with(
                A, timeout=5, allow_redirects=False, stream=True)
        else:
            session.get.assert_not_called()


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
        (requests.exceptions.ChunkedEncodingError("x"),
         "request error (ChunkedEncodingError)"),
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

    def test_redirect_without_location_stops(self):
        session = MagicMock()
        session.head.return_value = make_response(302)
        result = check_url(session, A, 5, 5)
        assert result['status_code'] == 302
        assert result['too_many_redirects'] is False
        assert session.head.call_count == 1

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

class TestIsBenignRedirect:

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

    @pytest.mark.parametrize("status, client, server", [
        (None, None, None),
        (399, None, None),
        (400, "400 bad request", None),
        (499, "499 error", None),
        (500, None, "500 internal server error"),
        (599, None, "599 error"),
    ])
    def test_status_range_boundaries(self, status, client, server):
        result = make_result(A, status=status)
        assert rule_client_error(result) == client
        assert rule_server_error(result) == server

    def test_redirect_ending_in_error(self):
        result = make_result(
            "https://a.gov/x", "https://a.gov/y", status=404,
            chain=[(301, "https://a.gov/x"), (404, "https://a.gov/y")])
        assert rule_client_error(result) == "redirect ends in 404 not found"


# ============================================================
# Yellow rules
# ============================================================

class TestYellowRules:

    @pytest.mark.parametrize("final, expected", [
        ("https://b.org/x", "cross-domain redirect (b.org)"),
        ("https://www.a.gov/x", None),
        ("https://sub.a.gov/x", "cross-domain redirect (sub.a.gov)"),
    ])
    def test_cross_domain(self, final, expected):
        result = make_result("https://a.gov/x", final)
        assert rule_cross_domain(result) == expected

    @pytest.mark.parametrize("original, final, flagged", [
        ("https://a.gov/study/1", "https://a.gov/", True),
        ("https://a.gov/study/1", "https://a.gov", True),
        ("https://a.gov/study/1", "https://a.gov/?ref=x", True),
        ("https://a.gov", "https://www.a.gov/", False),
        ("https://a.gov/x", "https://a.gov/y", False),
    ])
    def test_homepage_bounce(self, original, final, flagged):
        reason = rule_homepage_bounce(make_result(original, final))
        assert (reason is not None) is flagged

    def test_long_redirect_chain_threshold_is_exclusive(self):
        result = make_result(A, chain=[(301, A)] * 3 + [(200, A)])
        assert rule_long_redirect_chain(result, 2) == \
            "long redirect chain (3 hops)"
        assert rule_long_redirect_chain(result, 3) is None

    def test_slow_response_threshold_is_exclusive(self):
        assert rule_slow_response(
            make_result(A, response_time=2.5), 5, 0.5) is None
        assert rule_slow_response(
            make_result(A, response_time=2.6), 5, 0.5) == \
            "slow response (2.6s)"


# ============================================================
# classify_result
# ============================================================

class TestClassifyResult:

    @pytest.mark.parametrize("result, label, reason", [
        pytest.param(
            make_result(A, status=None, chain=[], error="timeout"),
            'R', "timeout", id="connection-error"),
        pytest.param(
            make_result(A, status=302, chain=[(302, A)] * 6, too_many=True),
            'R', "too many redirects (>5 hops)", id="too-many-redirects"),
        pytest.param(
            make_result("https://a.gov/x", "https://b.org/", status=404,
                        response_time=4.0),
            'R', "404 not found", id="red-beats-yellow"),
        pytest.param(
            make_result(A, status=204),
            'Y', "unexpected status 204", id="non-200-success"),
        pytest.param(
            make_redirect("https://a.gov/x", "https://a.gov/y"),
            'Y', "redirected to different page", id="same-host-move"),
        pytest.param(
            make_redirect("https://a.gov/x", "https://b.org/y"),
            'Y', "cross-domain redirect (b.org)",
            id="cross-domain-suppresses-generic"),
        pytest.param(
            make_redirect("https://a.gov/x", "https://a.gov/",
                          response_time=3.0),
            'Y', "redirect to homepage; slow response (3.0s)",
            id="multiple-yellow-reasons"),
        pytest.param(
            make_result("https://a.gov/x"),
            'G', "ok", id="no-redirect"),
        pytest.param(
            make_redirect("http://a.gov/x", "https://www.a.gov/x/"),
            'G', "ok (benign redirect)", id="benign-redirect"),
    ])
    def test_precedence_and_reasons(self, mock_config, result, label,
                                    reason):
        assert classify_result(result) == (label, reason)

    def test_uses_config_thresholds(self, mock_config):
        mock_config.URL_HEALTH_SLOW_THRESHOLD_FRACTION = 0.1
        result = make_result(A, response_time=1.0)
        assert classify_result(result)[0] == 'Y'


# ============================================================
# Reporting
# ============================================================

class TestReporting:

    def test_build_report_row_for_redirect(self, mock_config):
        result = make_redirect("http://a.gov", "https://a.gov",
                               response_time=0.12345)
        row = build_report_row(result, ["f | url | row 2",
                                        "f | url | row 5"])
        assert set(row) == set(REPORT_COLUMNS)
        assert row['url_changed'] == 'changed'
        assert row['redirect_chain'] == \
            "301 http://a.gov -> 200 https://a.gov"
        assert row['response_time_s'] == "0.123"
        assert row['occurrence_count'] == 2
        assert row['source_locations'] == "f | url | row 2; f | url | row 5"

    def test_build_report_row_for_connection_error(self, mock_config):
        result = make_result(A, status=None, chain=[], error="DNS failure")
        row = build_report_row(result, ["loc"])
        assert row['final_status_code'] == ''
        assert row['url_changed'] == 'not_changed'
        assert (row['classification'], row['reason']) == \
            ('R', "DNS failure")

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

    def test_write_report_to_bare_filename(self, tmp_path, monkeypatch):
        monkeypatch.chdir(tmp_path)
        write_report([], "report.csv")
        assert (tmp_path / "report.csv").exists()


# ============================================================
# Console summary grouping
# ============================================================

@pytest.fixture
def summary_rows():
    """Report rows spanning three files and all three classes."""
    return [
        {'original_url': 'https://a.gov', 'classification': 'R',
         'reason': 'timeout',
         'source_locations': 'p.tsv | link | row 2; p.tsv | link | row 3'},
        {'original_url': 'https://b.gov', 'classification': 'Y',
         'reason': 'redirect to homepage; slow response (3.2s)',
         'source_locations': 'p.tsv | link | row 4; d.tsv | url | row 2'},
        {'original_url': 'https://c.gov', 'classification': 'G',
         'reason': 'ok',
         'source_locations': 'd.tsv | url | row 3'},
        {'original_url': 'https://e.gov', 'classification': 'G',
         'reason': 'ok',
         'source_locations': 'g.tsv | url | row 2'},
    ]


class TestSummaryGrouping:

    @pytest.mark.parametrize("reason, expected", [
        ("slow response (3.2s)", ("slow response", "3.2s")),
        ("cross-domain redirect (b.org)", ("cross-domain redirect", "b.org")),
        ("404 not found", ("404 not found", "")),
    ])
    def test_split_reason(self, reason, expected):
        assert split_reason(reason) == expected

    def test_count_by_file(self, summary_rows):
        counts = count_by_file(summary_rows)
        assert counts == {
            'p.tsv': {'R': 1, 'Y': 1, 'G': 0},
            'd.tsv': {'R': 0, 'Y': 1, 'G': 1},
            'g.tsv': {'R': 0, 'Y': 0, 'G': 1},
        }

    def test_group_flagged_by_file(self, summary_rows):
        groups = group_flagged_by_file(summary_rows)
        assert groups['p.tsv']['R'] == {'timeout': [('https://a.gov', '')]}
        assert groups['d.tsv']['Y'] == {
            'redirect to homepage': [('https://b.gov', '')],
            'slow response': [('https://b.gov', '3.2s')],
        }
        assert 'G' not in groups['d.tsv']

    def test_print_summary_hierarchy(self, summary_rows, capsys):
        print_summary(summary_rows)
        out = capsys.readouterr().out
        assert out.index("URL HEALTH SUMMARY") < out.index("  p.tsv\n")
        assert out.index("  p.tsv\n") < out.index("  d.tsv\n")
        assert "RED (1 URLs)" in out
        assert "      - https://b.gov  [3.2s]" in out
        assert "No Red or Yellow URLs." in out.split("  g.tsv\n")[1]


# ============================================================
# check_url_health (end to end)
# ============================================================

class TestCheckUrlHealth:

    def test_checks_each_unique_url_once_and_writes_report(
            self, sample_tsv, tmp_path, mock_config):
        report = tmp_path / "reports" / "url_health.csv"
        mock_config.URL_HEALTH_REPORT_PATH = str(report)
        responses = {
            "https://a.gov": make_response(200),
            "https://b.gov": make_response(404),
            "https://c.gov": make_response(301, "https://other.org/"),
            "https://other.org/": make_response(200),
            "https://d.gov": make_response(200),
        }
        session = MagicMock()
        session.__enter__.return_value = session
        session.head.side_effect = lambda url, **kw: responses[url]
        session.get.side_effect = lambda url, **kw: responses[url]

        with patch("modules.check_url_health.requests.Session",
                   return_value=session):
            rows = check_url_health({sample_tsv: ["url", "links"]})

        head_urls = [c.args[0] for c in session.head.call_args_list]
        assert head_urls.count("https://a.gov") == 1
        session.headers.__setitem__.assert_called_once_with(
            'User-Agent', 'test-agent')
        assert len(rows) == 4

        with open(report, newline='', encoding='utf-8') as f:
            by_url = {r['original_url']: r for r in csv.DictReader(f)}
        assert {u: r['classification'] for u, r in by_url.items()} == {
            "https://a.gov": 'G',
            "https://b.gov": 'R',
            "https://c.gov": 'Y',
            "https://d.gov": 'G',
        }
        assert by_url["https://a.gov"]['occurrence_count'] == '2'
