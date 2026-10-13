"""
check_url_health.py
2026-10-13 ZD

Manually-triggered QA step that scans URL fields in TSV files and reports
link health as Red / Yellow / Green before files are loaded to the site.

For each unique URL across all configured targets, this script:
    1. Sends HEAD, falling back to GET if HEAD is rejected (4xx/5xx status)
    2. Follows redirects one hop at a time, capturing the full chain
    3. Classifies the result as R/Y/G using individually testable rules
       (Red first, then Yellow, else Green)

Requests run sequentially by design.

Input:
    - config.URL_HEALTH_TARGETS (dict of TSV path: list of URL columns)

Output:
    - Console summary of R/Y/G counts with all Red and Yellow URLs and reasons
    - config.URL_HEALTH_REPORT_PATH (CSV)
"""

import csv
import os
import re
import time
from http import HTTPStatus
from urllib.parse import urljoin, urlsplit

import requests
from tqdm import tqdm

import config


REDIRECT_CODES = {301, 302, 303, 307, 308}

# Large free-text fields can exceed csv default
CSV_FIELD_SIZE_LIMIT = 2**31 - 1

# Substrings that identify DNS failures
DNS_ERROR_MARKERS = [
    'nameresolutionerror',
    'name resolution',
    'getaddrinfo failed',
    'name or service not known',
    'nodename nor servname',
]

SEVERITY_ORDER = {'R': 0, 'Y': 1, 'G': 2}
LABEL_NAMES = {'R': 'Red', 'Y': 'Yellow', 'G': 'Green'}

REPORT_COLUMNS = [
    'original_url',
    'final_url',
    'url_changed',
    'final_status_code',
    'classification',
    'reason',
    'response_time_s',
    'redirect_chain',
    'occurrence_count',
    'source_locations',
]


# ---
# URL EXTRACTION

def split_url_field(value):
    """Split a URL field into individual URLs.

    Args:
        value (str | None): Raw cell value. May hold one URL or a
            semicolon-separated list of URLs.

    Returns:
        list[str]: Trimmed, non-empty URLs.
    """

    if not value:
        return []

    return [p.strip() for p in value.split(';') if p.strip()]


def read_target_urls(path, columns):
    """Extract URLs and locations from a TSV.

    Missing files, missing columns, and unreadable rows are logged and
    skipped.

    Args:
        path (str): Path to the TSV file.
        columns (list[str]): Columns that contain URLs.

    Returns:
        list[tuple[str, str]]: (url, location) pairs. Location is formatted
            as 'filename | column | row N', where row 1 is the header row.
    """

    if not os.path.exists(path):
        print(f"  WARNING: File not found, skipping: {path}")
        return []

    csv.field_size_limit(CSV_FIELD_SIZE_LIMIT)
    filename = os.path.basename(path)
    found = []

    try:
        with open(path, newline='', encoding='utf-8-sig') as f:
            reader = csv.DictReader(f, delimiter='\t')
            header = reader.fieldnames or []

            missing = [c for c in columns if c not in header]
            if missing:
                print(f"  WARNING: Columns not found in {path}: {missing}")
            present = [c for c in columns if c in header]

            row_num = 1
            while True:
                row_num += 1
                try:
                    row = next(reader)
                except StopIteration:
                    break
                except csv.Error as exc:
                    print(f"  WARNING: Skipping unreadable row {row_num} "
                          f"in {path}: {exc}")
                    continue

                for col in present:
                    for url in split_url_field(row.get(col)):
                        loc = f"{filename} | {col} | row {row_num}"
                        found.append((url, loc))

    except (OSError, UnicodeDecodeError) as exc:
        print(f"  WARNING: Could not read {path}: {exc}")

    return found


def collect_urls(targets):
    """Deduplicate URLs across all targets.

    Args:
        targets (dict): TSV path: list of URL column names.

    Returns:
        dict[str, list[str]]: Each unique URL mapped to every location where
            it occurs. Occurrence count is the length of the location list.
    """

    url_locations = {}

    for path, columns in targets.items():
        print(f"\n  Reading {path}")
        pairs = read_target_urls(path, columns)
        print(f"    {len(pairs)} URL values found")

        for url, loc in pairs:
            url_locations.setdefault(url, []).append(loc)

    return url_locations


# ---
# HTTP CHECKS

def build_session():
    """Create a requests session for checks.

    Returns:
        requests.Session: Session with the configured User-Agent header.
    """

    session = requests.Session()
    session.headers['User-Agent'] = config.URL_HEALTH_USER_AGENT
    return session


def request_once(session, url, timeout):
    """Request a URL without following redirects.

    Sends HEAD first. Falls back to GET if HEAD is rejected with a 4xx/5xx
    status (e.g. 405/501), since some servers mishandle HEAD.

    Args:
        session (requests.Session): Session.
        url (str): URL to request.
        timeout (float): Timeout in seconds.

    Returns:
        requests.Response: Final response.

    Raises:
        requests.RequestException: On a connection-level failure.
    """

    resp = session.head(url, timeout=timeout, allow_redirects=False)
    if resp.status_code < 400:
        return resp

    # Stream so the body is never downloaded
    resp = session.get(url, timeout=timeout, allow_redirects=False,
                       stream=True)
    resp.close()
    return resp


def describe_error(exc):
    """Map a requests exception to a reason.

    Args:
        exc (requests.RequestException): The raised exception.

    Returns:
        str: Short reason, e.g. 'timeout'.
    """

    errors = requests.exceptions
    text = str(exc).lower()

    if isinstance(exc, errors.SSLError):
        return 'SSL/TLS error'
    # ConnectTimeout is also a ConnectionError
    if isinstance(exc, errors.Timeout):
        return 'timeout'
    if isinstance(exc, (errors.MissingSchema, errors.InvalidSchema,
                        errors.InvalidURL)):
        return 'invalid or unsupported URL'
    if isinstance(exc, errors.ConnectionError):
        if any(m in text for m in DNS_ERROR_MARKERS):
            return 'DNS failure'
        if 'refused' in text:
            return 'connection refused'
        return 'connection error'

    return f"request error ({type(exc).__name__})"


def check_url(session, url, timeout, max_hops):
    """Check a single URL, following redirects.

    Never raises. Connection-level failures and unexpected errors are
    captured in the result's 'error' field.

    Args:
        session (requests.Session): Session.
        url (str): URL to check.
        timeout (float): Per-request timeout in seconds.
        max_hops (int): Maximum redirects to follow before flagging.

    Returns:
        dict: Result with keys 'original_url', 'final_url', 'status_code'
            (int or None), 'chain' (list of (status, url) tuples),
            'response_time' (seconds), 'error' (str or None), and
            'too_many_redirects' (bool).
    """

    result = {
        'original_url': url,
        'final_url': url,
        'status_code': None,
        'chain': [],
        'response_time': 0.0,
        'error': None,
        'too_many_redirects': False,
    }

    current = url
    start = time.perf_counter()

    try:
        for hop in range(max_hops + 1):
            resp = request_once(session, current, timeout)
            code = resp.status_code
            result['chain'].append((code, current))
            result['final_url'] = current
            result['status_code'] = code

            location = resp.headers.get('Location')
            if code not in REDIRECT_CODES or not location:
                break
            if hop == max_hops:
                result['too_many_redirects'] = True
                break
            current = urljoin(current, location)

    except requests.RequestException as exc:
        result['error'] = describe_error(exc)

    except Exception as exc:
        tqdm.write(f"  WARNING: Unexpected error checking {url}: {exc}")
        result['error'] = f"unexpected error ({type(exc).__name__})"

    result['response_time'] = time.perf_counter() - start
    return result


# ---
# URL COMPARISON HELPERS

def count_hops(result):
    """Count redirects followed for a result.

    Args:
        result (dict): Output of check_url.

    Returns:
        int: Number of redirect hops.
    """

    return max(len(result['chain']) - 1, 0)


def normalize_host(url):
    """Lowercase host with any 'www.' removed.

    Args:
        url (str): URL to parse.

    Returns:
        str: Normalized host, or '' if none.
    """

    host = urlsplit(url).hostname or ''
    return host.removeprefix('www.')


def is_root_path(url):
    """Check if a URL points to the site root.

    Args:
        url (str): URL to parse.

    Returns:
        bool: True if the path is '' or '/'.
    """

    return urlsplit(url).path in ('', '/')


def is_benign_redirect(original, final):
    """Check if a redirect is only cosmetic.

    Benign differences are http -> https, trailing slash added/removed, and
    'www.' added/removed. Anything else is not benign.

    Args:
        original (str): Original URL.
        final (str): Final URL.

    Returns:
        bool: True if the URLs differ only in benign ways (or not at all).
    """

    a = urlsplit(original)
    b = urlsplit(final)

    scheme_ok = (a.scheme == b.scheme
                 or (a.scheme == 'http' and b.scheme == 'https'))
    same_host = normalize_host(original) == normalize_host(final)
    same_path = a.path.rstrip('/') == b.path.rstrip('/')

    return scheme_ok and same_host and same_path and a.query == b.query


def describe_status(result):
    """Build a short reason for an HTTP status.

    Args:
        result (dict): Output of check_url.

    Returns:
        str: e.g. '404 not found', or 'redirect ends in 404 not found' when
            redirects were followed.
    """

    code = result['status_code']
    try:
        phrase = HTTPStatus(code).phrase.lower()
    except ValueError:
        phrase = 'error'

    text = f"{code} {phrase}"
    if count_hops(result) > 0:
        return f"redirect ends in {text}"
    return text


# ---
# RED RULES
# Each rule returns a reason string, or None

def rule_connection_error(result):
    """Red: DNS, timeout, refused, SSL, etc.

    Args:
        result (dict): Output of check_url.

    Returns:
        str | None: Reason, or None if OK.
    """

    return result['error']


def rule_too_many_redirects(result):
    """Red: redirects exceeded the hop limit.

    Args:
        result (dict): Output of check_url.

    Returns:
        str | None: Reason, or None if OK.
    """

    if result['too_many_redirects']:
        return f"too many redirects (>{count_hops(result)} hops)"
    return None


def rule_client_error(result):
    """Red: final status is 4xx.

    Args:
        result (dict): Output of check_url.

    Returns:
        str | None: Reason, or None if OK.
    """

    code = result['status_code']
    if code is not None and 400 <= code < 500:
        return describe_status(result)
    return None


def rule_server_error(result):
    """Red: final status is 5xx.

    Args:
        result (dict): Output of check_url.

    Returns:
        str | None: Reason, or None if OK.
    """

    code = result['status_code']
    if code is not None and 500 <= code < 600:
        return describe_status(result)
    return None


# ---
# YELLOW RULES
# Each rule returns a reason string, or None

def rule_unexpected_status(result):
    """Yellow: final status is not 200.

    Covers non-error statuses such as 204 or a 3xx without a Location header.

    Args:
        result (dict): Output of check_url.

    Returns:
        str | None: Reason, or None if OK.
    """

    code = result['status_code']
    if code != 200:
        return f"unexpected status {code}"
    return None


def rule_cross_domain(result):
    """Yellow: final domain differs.

    Hosts are compared after lowercasing and removing 'www.'. Any other host
    change, including subdomains, is cross-domain.

    Args:
        result (dict): Output of check_url.

    Returns:
        str | None: Reason, or None if OK.
    """

    orig = normalize_host(result['original_url'])
    final = normalize_host(result['final_url'])
    if orig != final:
        return f"cross-domain redirect ({final})"
    return None


def rule_homepage_bounce(result):
    """Yellow: non-root URL lands on root.

    A redirect from a specific page to '/' is likely a silent bounce to the
    homepage.

    Args:
        result (dict): Output of check_url.

    Returns:
        str | None: Reason, or None if OK.
    """

    orig_root = is_root_path(result['original_url'])
    final_root = is_root_path(result['final_url'])
    if not orig_root and final_root:
        return "redirect to homepage"
    return None


def rule_non_benign_redirect(result):
    """Yellow: redirect is not cosmetic.

    Flags same-domain redirects that change the path or query (e.g. a moved
    page).

    Args:
        result (dict): Output of check_url.

    Returns:
        str | None: Reason, or None if OK.
    """

    if count_hops(result) == 0:
        return None
    if is_benign_redirect(result['original_url'], result['final_url']):
        return None
    return "redirected to different page"


def rule_long_redirect_chain(result, max_hops):
    """Yellow: redirect chain is long.

    Args:
        result (dict): Output of check_url.
        max_hops (int): Hops allowed before flagging.

    Returns:
        str | None: Reason, or None if OK.
    """

    hops = count_hops(result)
    if hops > max_hops:
        return f"long redirect chain ({hops} hops)"
    return None


def rule_slow_response(result, timeout, fraction):
    """Yellow: response time is slow.

    Args:
        result (dict): Output of check_url.
        timeout (float): Timeout in seconds.
        fraction (float): Fraction of timeout above which a response is
            slow.

    Returns:
        str | None: Reason, or None if OK.
    """

    elapsed = result['response_time']
    if elapsed > timeout * fraction:
        return f"slow response ({elapsed:.1f}s)"
    return None


# ---
# CLASSIFICATION

def classify_result(result):
    """Classify a result as Red/Yellow/Green.

    Red rules are evaluated first, then Yellow, else Green. All matching
    reasons within the winning tier are reported.

    Args:
        result (dict): Output of check_url.

    Returns:
        tuple[str, str]: ('R'|'Y'|'G', reason)
    """

    red = [
        rule_connection_error(result),
        rule_too_many_redirects(result),
        rule_client_error(result),
        rule_server_error(result),
    ]
    red = [r for r in red if r]
    if red:
        return 'R', '; '.join(red)

    redirect_flags = [
        rule_cross_domain(result),
        rule_homepage_bounce(result),
    ]
    # Generic flag only if no specific one
    if not any(redirect_flags):
        redirect_flags.append(rule_non_benign_redirect(result))

    yellow = [
        rule_unexpected_status(result),
        *redirect_flags,
        rule_long_redirect_chain(
            result, config.URL_HEALTH_YELLOW_REDIRECT_HOPS),
        rule_slow_response(
            result, config.URL_HEALTH_TIMEOUT_SECONDS,
            config.URL_HEALTH_SLOW_THRESHOLD_FRACTION),
    ]
    yellow = [r for r in yellow if r]
    if yellow:
        return 'Y', '; '.join(yellow)

    if count_hops(result) > 0:
        return 'G', 'ok (benign redirect)'
    return 'G', 'ok'


# ---
# REPORTING

def format_chain(chain):
    """Format a redirect chain compactly.

    Args:
        chain (list[tuple[int, str]]): (status, url) for each hop.

    Returns:
        str: e.g. '301 http://a -> 200 https://a'
    """

    return ' -> '.join(f"{code} {url}" for code, url in chain)


def build_report_row(result, locations):
    """Build one CSV report row for a URL.

    Args:
        result (dict): Output of check_url.
        locations (list[str]): Every location where the URL occurs.

    Returns:
        dict: Row keyed by REPORT_COLUMNS.
    """

    label, reason = classify_result(result)
    code = result['status_code']
    changed = result['original_url'] != result['final_url']

    return {
        'original_url': result['original_url'],
        'final_url': result['final_url'],
        'url_changed': 'changed' if changed else 'not_changed',
        'final_status_code': '' if code is None else code,
        'classification': label,
        'reason': reason,
        'response_time_s': f"{result['response_time']:.3f}",
        'redirect_chain': format_chain(result['chain']),
        'occurrence_count': len(locations),
        'source_locations': '; '.join(locations),
    }


def location_files(row):
    """Get unique file names from a report row's source locations.

    Args:
        row (dict): Report row.

    Returns:
        list[str]: File names in order of first appearance.
    """

    locs = row['source_locations'].split('; ')
    return list(dict.fromkeys(loc.split(' | ')[0] for loc in locs))


def split_reason(reason):
    """Split a single reason into its type and optional detail.

    Args:
        reason (str): e.g. 'slow response (3.2s)'.

    Returns:
        tuple[str, str]: (type, detail), e.g. ('slow response', '3.2s').
            Detail is '' when the reason has no parenthetical.
    """

    match = re.fullmatch(r'(.+?) \((.+)\)', reason)
    if match:
        return match.group(1), match.group(2)
    return reason, ''


def count_by_file(rows):
    """Count unique URLs per classification for each source file.

    Args:
        rows (list[dict]): Report rows.

    Returns:
        dict[str, dict[str, int]]: File name to {'R': n, 'Y': n, 'G': n}.
    """

    counts = {}
    for row in rows:
        for name in location_files(row):
            file_counts = counts.setdefault(
                name, dict.fromkeys(LABEL_NAMES, 0))
            file_counts[row['classification']] += 1
    return counts


def group_flagged_by_file(rows):
    """Group Red and Yellow URLs by file, class, then reason type.

    A URL with multiple reasons is listed under each reason type.

    Args:
        rows (list[dict]): Report rows.

    Returns:
        dict: file -> label -> type -> list of (url, detail) tuples.
    """

    groups = {}
    for row in rows:
        label = row['classification']
        if label == 'G':
            continue

        for part in row['reason'].split('; '):
            issue, detail = split_reason(part)
            for name in location_files(row):
                by_type = groups.setdefault(name, {}).setdefault(label, {})
                by_type.setdefault(issue, []).append(
                    (row['original_url'], detail))
    return groups


def format_count_row(name, counts, width):
    """Format one row of the R/Y/G summary table.

    Args:
        name (str): Row label (file name or total).
        counts (dict[str, int]): {'R': n, 'Y': n, 'G': n}.
        width (int): Width of the file name column.

    Returns:
        str: Formatted table row.
    """

    cells = ''.join(f"{counts[label]:>8}" for label in LABEL_NAMES)
    return f"  {name:<{width}}{cells}"


def print_count_table(rows, file_counts):
    """Print the R/Y/G summary table with per-file and total rows.

    Args:
        rows (list[dict]): Report rows.
        file_counts (dict): Output of count_by_file.
    """

    totals = dict.fromkeys(LABEL_NAMES, 0)
    for row in rows:
        totals[row['classification']] += 1

    total_label = 'All files (unique URLs)'
    width = max([len(total_label), *map(len, file_counts)]) + 2
    divider = f"  {'-' * (width + 8 * len(LABEL_NAMES))}"
    header = ''.join(f"{name:>8}" for name in LABEL_NAMES.values())

    print(f"\n{'=' * 79}")
    print("  URL HEALTH SUMMARY")
    print(f"{'=' * 79}")
    print(f"  {'File':<{width}}{header}")
    print(divider)
    for name, counts in file_counts.items():
        print(format_count_row(name, counts, width))
    print(divider)
    print(format_count_row(total_label, totals, width))


def print_file_section(name, by_label):
    """Print flagged URLs for one file, grouped by class then type.

    Args:
        name (str): Source file name.
        by_label (dict): label -> type -> list of (url, detail) tuples.
    """

    print(f"\n{'=' * 79}")
    print(f"  {name}")
    print(f"{'=' * 79}")

    if not by_label:
        print("  No Red or Yellow URLs.")
        return

    for label in ('R', 'Y'):
        by_type = by_label.get(label)
        if not by_type:
            continue

        unique = {url for items in by_type.values() for url, _ in items}
        print(f"\n  {LABEL_NAMES[label].upper()} ({len(unique)} URLs)")

        ordered = sorted(by_type.items(),
                         key=lambda kv: (-len(kv[1]), kv[0]))
        for issue, items in ordered:
            print(f"    {issue} ({len(items)})")
            for url, detail in items:
                suffix = f"  [{detail}]" if detail else ''
                print(f"      - {url}{suffix}")


def print_summary(rows):
    """Print the summary table, then flagged URLs grouped by file.

    Args:
        rows (list[dict]): Report rows.
    """

    file_counts = count_by_file(rows)
    groups = group_flagged_by_file(rows)

    print_count_table(rows, file_counts)
    for name in file_counts:
        print_file_section(name, groups.get(name, {}))
    print(f"\n{'=' * 79}")


def write_report(rows, path):
    """Write report rows to CSV, Red first.

    Args:
        rows (list[dict]): Report rows.
        path (str): Output CSV path.
    """

    output_dir = os.path.dirname(path)
    if output_dir:
        os.makedirs(output_dir, exist_ok=True)

    ordered = sorted(rows,
                     key=lambda r: SEVERITY_ORDER[r['classification']])

    with open(path, 'w', newline='', encoding='utf-8') as f:
        writer = csv.DictWriter(f, fieldnames=REPORT_COLUMNS)
        writer.writeheader()
        writer.writerows(ordered)


# ---
# MAIN

def check_url_health(targets=None):
    """Check health of all URLs in TSV targets.

    Reads:
        config.URL_HEALTH_TARGETS (default)

    Writes:
        config.URL_HEALTH_REPORT_PATH

    Args:
        targets (dict, optional): TSV path: list of URL columns. Defaults
            to config.URL_HEALTH_TARGETS.

    Returns:
        list[dict]: One report row per unique URL.
    """

    print(f"\n---\nQA - URL Health Check:\n"
          f"Checking URLs in TSV targets...\n---")

    if targets is None:
        targets = config.URL_HEALTH_TARGETS

    url_locations = collect_urls(targets)
    total = sum(len(locs) for locs in url_locations.values())
    print(f"\nFound {total} URL occurrences "
          f"({len(url_locations)} unique)")

    rows = []
    with build_session() as session:
        for url in tqdm(url_locations, desc="Checking URLs", ncols=80):
            result = check_url(session, url,
                               config.URL_HEALTH_TIMEOUT_SECONDS,
                               config.URL_HEALTH_MAX_REDIRECT_HOPS)
            rows.append(build_report_row(result, url_locations[url]))

    print_summary(rows)

    output_path = config.URL_HEALTH_REPORT_PATH
    write_report(rows, output_path)
    print(f"\nSuccess! URL health report saved "
          f"to {output_path}\n")

    return rows


# Run module as a standalone script when called directly
if __name__ == "__main__":

    print(f"Running {os.path.basename(__file__)} "
          f"as standalone module...")

    check_url_health()
