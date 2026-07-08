"""
test_gather_grant_data.py
2024-12-19 ZD

Pytest test suite for the `gather_grant_data.py` module.
"""

import os
import pandas as pd
import pytest
import requests
from unittest.mock import patch, MagicMock

from modules.gather_grant_data import (
    get_nih_reporter_grants,
    concatenate_full_names,
    format_name_column,
    extract_total_cost,
    clean_abstract,
    format_organization_columns,
    clean_grants_data,
    gather_grant_data,
)


# ============================================================
# Shared Fixtures
# ============================================================

@pytest.fixture
def mock_grant_response():
    """A mock NIH RePORTER grant result with all relevant fields."""
    return {
        "meta": {
            "search_id": "R01CA123456",
            "total": 1,
            "offset": 0,
            "limit": 500,
            "sort_field": "appl_id",
            "sort_order": "desc",
        },
        "results": [
            {
                "appl_id": 1234567,
                "fiscal_year": 2010,
                "project_num": "R01CA123456-01",
                "organization": {
                    "org_name": "Test University",
                    "org_city": "Test City",
                    "org_state": "TS",
                    "org_country": "USA",
                },
                "award_amount": 600000,
                "principal_investigators": [
                    {"full_name": "John Doe"},
                    {"full_name": "Jane Smith"},
                ],
                "program_officers": [
                    {"full_name": "Alice Johnson"},
                ],
                "agency_ic_fundings": [
                    {"code": "CA", "total_cost": 500000},
                    {"code": "Other", "total_cost": 100000},
                ],
                "project_start_date": "2011-01-01T12:01:00Z",
                "project_end_date": "2012-01-01T12:01:00Z",
                "opportunity_number": "RFA-CA-01-123",
                "award_notice_date": "2010-01-01T12:01:00Z",
                "core_project_num": "R01CA123456",
                "pref_terms": "keyword1;keyword2",
                "abstract_text": "Sample\xad abstract  with\xa0spaces.",
                "project_title": "Sample Project Title",
                "api_source_search": "award_R01CA123456",
            }
        ],
    }


@pytest.fixture
def mock_cleaned_grant():
    """Expected cleaned output of mock_grant_response after clean_grants_data."""
    return pd.DataFrame({
        "grant_id": ["R01CA123456-01"],
        "queried_project_id": ["R01CA123456"],
        "application_id": [1234567],
        "fiscal_year": [2010],
        "project_title": ["Sample Project Title"],
        "abstract_text": ["Sample abstract with spaces."],
        "keywords": ["keyword1;keyword2"],
        "principal_investigators": ["John Doe; Jane Smith"],
        "program_officers": ["Alice Johnson"],
        "award_amount": [600000],
        "nci_funded_amount": [500000],
        "award_notice_date": ["2010-01-01T12:01:00Z"],
        "project_start_date": ["2011-01-01T12:01:00Z"],
        "project_end_date": ["2012-01-01T12:01:00Z"],
        "opportunity_number": ["RFA-CA-01-123"],
        "api_source_search": ["award_R01CA123456"],
        "org_name": ["Test University"],
        "org_city": ["Test City"],
        "org_state": ["TS"],
        "org_country": ["USA"],
    })


@pytest.fixture
def mock_api_response():
    """Factory for minimal mock API responses."""
    def _create_response(status_code, results=None, total=0):
        response = MagicMock()
        response.status_code = status_code
        response.json.return_value = {
            'results': results or [],
            'meta': {'total': total},
        }
        return response
    return _create_response


# ============================================================
# Data Cleaning Helpers
# ============================================================

class TestConcatenateFullNames:

    def test_joins_names_with_semicolons(self):
        names = [{"full_name": "John Doe"}, {"full_name": "Jane Smith"}]
        assert concatenate_full_names(names) == "John Doe; Jane Smith"

    def test_empty_and_none(self):
        assert concatenate_full_names([]) == ''
        assert concatenate_full_names(None) == ''


class TestFormatNameColumn:

    def test_title_cases_and_collapses_whitespace(self):
        assert format_name_column("john  doe") == "John Doe"
        assert format_name_column("JANE DOE") == "Jane Doe"
        assert format_name_column("mary ann smith") == "Mary Ann Smith"


class TestExtractTotalCost:

    def test_extracts_nci_ca_funding(self):
        fundings = [
            {"code": "Other_IC", "total_cost": 100000},
            {"code": "CA", "total_cost": 500000},
            {"code": "ABC", "total_cost": 200000},
        ]
        assert extract_total_cost(fundings) == 500000

    def test_no_nci_funding_returns_zero(self):
        fundings = [
            {"code": "Other", "total_cost": 100000},
            {"code": "ABC", "total_cost": 200000},
        ]
        assert extract_total_cost(fundings) == 0

    def test_none_returns_zero(self):
        assert extract_total_cost(None) == 0


class TestCleanAbstract:

    @pytest.mark.parametrize("input_text,expected", [
        ("Example \xad two  spaces", "Example two spaces"),
        ("Example    three spaces", "Example three spaces"),
        ("Example\xa0non-breaking space", "Example non-breaking space"),
        ("Example \nnewline", "Example newline"),
        (None, None),
        ('', ''),
        (123, 123),
    ])
    def test_cleaning(self, input_text, expected):
        assert clean_abstract(input_text) == expected


class TestFormatOrganizationColumns:

    def test_extracts_subfields_and_drops_original(self):
        """Uses the same subfield names as the real config
        (org_name, org_city, org_state, org_country)."""
        df = pd.DataFrame({
            'organization': [{
                'org_name': 'Test University',
                'org_city': 'Test City',
                'org_state': 'TS',
                'org_country': 'USA',
            }],
        })

        result = format_organization_columns(
            df, 'organization',
            ['org_name', 'org_city', 'org_state', 'org_country'])

        assert 'organization' not in result.columns
        assert result['org_name'].iloc[0] == 'Test University'
        assert result['org_city'].iloc[0] == 'Test City'
        assert result['org_state'].iloc[0] == 'TS'
        assert result['org_country'].iloc[0] == 'USA'


# ============================================================
# NIH RePORTER API (mocked)
# ============================================================

class TestGetNihReporterGrants:

    @pytest.mark.parametrize("search_type,search_value", [
        ('award', 'R01CA123456'),
        ('nofo', 'RFA-CA-01-123'),
    ])
    def test_single_page_response(self, mock_api_response,
                                   search_type, search_value):
        with patch('modules.gather_grant_data.requests.post') as mock_post:
            mock_post.return_value = mock_api_response(
                200, [{'project_num': search_value}], 1)

            grants, failures = get_nih_reporter_grants(
                [search_value], search_type)

            assert len(grants) == 1
            assert len(failures) == 0
            assert grants[0]['api_source_search'] == (
                f'{search_type}_{search_value}')

    @pytest.mark.parametrize("search_type,search_value", [
        ('award', 'R01CA123456'),
        ('nofo', 'RFA-CA-01-123'),
    ])
    def test_pagination(self, mock_api_response, search_type, search_value):
        with patch('modules.gather_grant_data.requests.post') as mock_post:
            mock_post.side_effect = [
                mock_api_response(200, [{'id': '1'}] * 500, 750),
                mock_api_response(200, [{'id': '2'}] * 250, 750),
            ]

            grants, failures = get_nih_reporter_grants(
                [search_value], search_type)

            assert len(grants) == 750
            assert mock_post.call_count == 2

    @pytest.mark.parametrize("search_type,search_value", [
        ('award', 'R01CA123456'),
        ('nofo', 'RFA-CA-01-123'),
    ])
    def test_non_500_error(self, mock_api_response,
                            search_type, search_value):
        with patch('modules.gather_grant_data.requests.post') as mock_post:
            mock_post.return_value = mock_api_response(123)
            grants, failures = get_nih_reporter_grants(
                [search_value], search_type)

            assert failures[search_value]['failure_type'] == '123 API Error'

    @pytest.mark.parametrize("search_type,search_value", [
        ('award', 'R01CA123456'),
        ('nofo', 'RFA-CA-01-123'),
    ])
    def test_request_exception(self, mock_api_response, search_type,
                                search_value, capsys):
        with patch('modules.gather_grant_data.requests.post') as mock_post:
            mock_post.side_effect = requests.exceptions.RequestException()
            grants, failures = get_nih_reporter_grants(
                [search_value], search_type)

            assert "An error occurred" in capsys.readouterr().out
            assert len(grants) == 0

    @pytest.mark.parametrize("search_type,search_value", [
        ('award', 'R01CA123456'),
        ('nofo', 'RFA-CA-01-123'),
    ])
    def test_500_retry_success(self, mock_api_response,
                                search_type, search_value, capsys):
        with (patch('modules.gather_grant_data.requests.post') as mock_post,
              patch('modules.gather_grant_data.sleep') as mock_sleep):
            mock_sleep.side_effect = lambda x: None
            mock_post.side_effect = [
                mock_api_response(500),
                mock_api_response(200, [{'id': 1}], 1),
            ]

            grants, failures = get_nih_reporter_grants(
                [search_value], search_type)

            assert len(grants) == 1
            assert len(failures) == 0
            assert mock_post.call_count == 2
            assert mock_sleep.call_count == 1

            captured = capsys.readouterr().out
            assert "Attempt 1/5" in captured
            assert "Attempt 2/5" not in captured

    @pytest.mark.parametrize("search_type,search_value", [
        ('award', 'R01CA123456'),
        ('nofo', 'RFA-CA-01-123'),
    ])
    def test_500_retry_exhausted(self, mock_api_response,
                                  search_type, search_value, capsys):
        with (patch('modules.gather_grant_data.requests.post') as mock_post,
              patch('modules.gather_grant_data.sleep') as mock_sleep):
            mock_sleep.side_effect = lambda x: None
            mock_post.side_effect = [mock_api_response(500)] * 5 + [
                mock_api_response(200, [{'id': 1}], 1),  # never reached
            ]

            grants, failures = get_nih_reporter_grants(
                [search_value], search_type)

            assert len(grants) == 0
            assert len(failures) == 1
            assert "500 Error: Likely too many results" in (
                failures[search_value]['failure_type'])
            assert mock_post.call_count == 5
            assert mock_sleep.call_count == 5
            assert "Attempt 5/5" in capsys.readouterr().out

    @pytest.mark.parametrize("search_type,search_value", [
        ('award', ['', 'R01CA123456']),
        ('nofo', ['', 'RFA-CA-01-123']),
    ])
    def test_skips_empty_search_values(self, mock_api_response,
                                        search_type, search_value, capsys):
        with patch('modules.gather_grant_data.requests.post') as mock_post:
            mock_post.return_value = mock_api_response(200)
            grants, failures = get_nih_reporter_grants(
                search_value, search_type, print_meta=True)

            assert f"Blank {search_type} value" in capsys.readouterr().out
            assert mock_post.call_count == 1

    def test_invalid_search_type_raises(self, mock_api_response):
        with patch('modules.gather_grant_data.requests.post') as mock_post:
            mock_post.return_value = mock_api_response(200)

            with pytest.raises(ValueError, match="Invalid search type"):
                get_nih_reporter_grants(['R01CA123456'], 'wrong_type')

            assert mock_post.call_count == 0

    @pytest.mark.parametrize("search_type,search_value", [
        ('award', 'R01CA123456'),
        ('nofo', 'RFA-CA-01-123'),
    ])
    def test_blank_results_reported_as_failure(self, mock_api_response,
                                                search_type, search_value):
        with patch('modules.gather_grant_data.requests.post') as mock_post:
            mock_post.return_value = mock_api_response(200, results=[], total=0)
            grants, failures = get_nih_reporter_grants(
                [search_value], search_type)

            assert len(grants) == 0
            assert len(failures) == 1
            assert failures[search_value]['failure_type'] == "No results found"


# ============================================================
# clean_grants_data
# ============================================================

class TestCleanGrantsData:

    def test_against_expected_output(self, mock_grant_response,
                                     mock_cleaned_grant):
        """Verify full cleaning pipeline produces the expected DataFrame."""
        result = clean_grants_data(mock_grant_response['results'])

        assert isinstance(result, pd.DataFrame)
        assert len(result) == 1
        pd.testing.assert_frame_equal(result, mock_cleaned_grant)


# ============================================================
# Integration: gather_grant_data
# ============================================================

class TestGatherGrantData:

    def test_orchestrator_produces_output(self, mock_grant_response, tmp_path):
        """Mocks only the API call; lets clean_grants_data run with real
        config to verify the full orchestration loop."""
        programs_df = pd.DataFrame({
            'program_name': ['Test Program'],
            'program_id': ['test_program'],
            'nofo': ['RFA-CA-01-123'],
            'award': [None],
        })

        grant_path = tmp_path / 'grant.csv'
        failed_path = tmp_path / 'failedSearches.csv'

        with (
            patch('modules.gather_grant_data.get_nih_reporter_grants')
                as mock_api,
            patch('modules.gather_grant_data.config.GRANTS_INTERMED_PATH',
                  str(grant_path)),
            patch('modules.gather_grant_data.config.FAILED_GRANT_SEARCH_REPORT',
                  str(failed_path)),
        ):
            # Award search returns nothing (award is None → empty list)
            # NOFO search returns one grant
            mock_api.side_effect = [
                ([], {}),                                       # awards
                (mock_grant_response['results'], {}),           # nofos
            ]

            result_df = gather_grant_data(programs_df)

            assert not result_df.empty
            assert grant_path.exists()
            assert failed_path.exists()
            assert 'program.program_id' in result_df.columns
            assert result_df.iloc[0]['program.program_id'] == 'test_program'


# ============================================================
# Live API smoke tests — run by default, excluded in CI with -m "not live_api"
# ============================================================

@pytest.mark.live_api
@pytest.mark.xfail(reason="Live NIH RePORTER API — may be slow or unavailable",
                   raises=AssertionError)
class TestNIHReporterGrantsAPILive:
    """Live smoke tests for the NIH RePORTER Projects/Grants API.
    Excluded in CI with -m "not live_api". To skip locally, use the same flag.
    """

    def test_reporter_grants_api_reachable_and_returns_expected_schema(self):
        """Verify the Projects API is up and the response schema matches
        what get_nih_reporter_grants expects."""
        response = requests.post(
            "https://api.reporter.nih.gov/v2/projects/search",
            json={
                "criteria": {
                    "opportunity_numbers": ["RFA-CA-21-038"],
                    "exclude_subprojects": True,
                    "fiscal_years": ["2022"],
                },
                "offset": 0,
                "limit": 5,
            },
            headers={"accept": "application/json",
                     "Content-Type": "application/json"},
            timeout=30,
        )

        assert response.status_code == 200, (
            f"NIH RePORTER Grants API returned {response.status_code}")

        data = response.json()

        assert "meta" in data, "Response missing 'meta' key"
        assert "results" in data, "Response missing 'results' key"

        meta = data["meta"]
        for key in ("total", "offset", "limit"):
            assert key in meta, f"Meta missing '{key}' field"

        assert meta["total"] > 0, (
            "Known NOFO RFA-CA-21-038 returned 0 results — "
            "API may have changed or data may have been removed")

        result = data["results"][0]
        for field in ("project_num", "fiscal_year", "organization",
                      "principal_investigators", "program_officers",
                      "agency_ic_fundings", "project_start_date",
                      "project_end_date", "opportunity_number",
                      "award_notice_date", "core_project_num",
                      "abstract_text", "project_title"):
            assert field in result, (
                f"Result missing '{field}' — API schema may have changed")

    def test_reporter_grants_api_returns_results_for_known_nofo(self):
        """Verify a known NOFO returns at least one grant."""
        results, failed = get_nih_reporter_grants(
            search_type="nofo",
            search_values=["RFA-CA-21-038"],
        )
        assert len(results) > 0, (
            "Known NOFO RFA-CA-21-038 returned no grants from RePORTER API")
        assert not failed, (
            f"Unexpected failures for known NOFO: {failed}")

