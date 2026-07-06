"""
test_gather_publication_data.py
2026-07-06 ZD

Pytest test suite for the gather_publication_data.py module.
Tests cover data formatting helpers, mocked API interactions, and live
API smoke tests that verify external services are reachable and returning
expected schemas.
"""

import os
import sys
from datetime import datetime

import pandas as pd
import pytest
from unittest.mock import patch, MagicMock

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from modules.gather_publication_data import (
    get_pmids_from_nih_reporter_api,
    get_pubmed_info_from_pmid,
    format_authors,
    format_publication_date,
    convert_month,
    extract_medline_date_components,
    merge_pubmed_icite_pmid_data,
    merge_and_clean_project_pmid_info,
)


# ============================================================
# Fixtures
# ============================================================

@pytest.fixture
def sample_pubmed_df():
    """Minimal PubMed-sourced publication DataFrame."""
    return pd.DataFrame({
        "pmid": [11111111, 22222222, 33333333],
        "title": ["Title A", "Title B", None],
        "authors": ["Author A", None, None],
        "publication_date": pd.to_datetime(
            ["2020-06-15", "2019-03-01", None]),
    })


@pytest.fixture
def sample_icite_df():
    """Minimal iCite-sourced DataFrame."""
    return pd.DataFrame({
        "pmid": [11111111, 22222222, 33333333],
        "title": ["Title A icite", "Title B icite", "Title C icite"],
        "authors": ["Author A icite", "Author B icite", "Author C icite"],
        "year": [2020, 2019, 2021],
        "citation_count": [10, 5, 0],
        "relative_citation_ratio": [1.5, 0.8, 0.0],
    })


@pytest.fixture
def sample_pmids_df():
    """DataFrame mapping projects to PMIDs."""
    return pd.DataFrame({
        "coreproject": ["R01CA111111", "R01CA111111", "U01CA222222"],
        "pmid": [11111111, 22222222, 33333333],
    })


# ============================================================
# format_authors
# ============================================================

class TestFormatAuthors:
    """Tests for author name formatting."""

    def test_standard_authors(self):
        authors = [
            {"LastName": "Smith", "ForeName": "John"},
            {"LastName": "Doe", "ForeName": "Jane"},
        ]
        result = format_authors(authors)
        assert result == "John Smith;Jane Doe"

    def test_collective_name(self):
        authors = [
            {"LastName": "Smith", "ForeName": "John"},
            {"CollectiveName": "TCGA Research Network"},
        ]
        result = format_authors(authors)
        assert "John Smith" in result
        assert "TCGA Research Network" in result

    def test_empty_author_list(self):
        assert format_authors([]) == ""


# ============================================================
# convert_month
# ============================================================

class TestConvertMonth:
    """Tests for month string/int conversion."""

    def test_integer_passthrough(self):
        assert convert_month(8) == 8

    def test_numeric_string(self):
        assert convert_month("08") == 8

    def test_abbreviation(self):
        assert convert_month("Aug") == 8

    def test_unknown_string_defaults_to_1(self):
        assert convert_month("Unknown") == 1


# ============================================================
# extract_medline_date_components
# ============================================================

class TestExtractMedlineDateComponents:
    """Tests for parsing non-standard Medline date strings."""

    def test_full_date(self):
        result = extract_medline_date_components("2020 Jun 15")
        assert result == datetime(2020, 6, 15)

    def test_year_and_month_only(self):
        result = extract_medline_date_components("2020 Mar")
        assert result.year == 2020
        assert result.month == 3

    def test_year_only(self):
        result = extract_medline_date_components("2020")
        assert result.year == 2020

    def test_season_string(self):
        result = extract_medline_date_components("2019 Winter")
        assert result == datetime(2019, 12, 1)

    def test_range_format(self):
        """Medline sometimes uses ranges like '2020 Jan-Feb'."""
        result = extract_medline_date_components("2020 Jan-Feb")
        assert result.year == 2020
        assert result.month == 1


# ============================================================
# format_publication_date
# ============================================================

class TestFormatPublicationDate:
    """Tests for standardizing publication date dicts."""

    def test_full_date_dict(self):
        pub_date = {"Year": "2021", "Month": "Aug", "Day": "15"}
        result = format_publication_date(pub_date)
        assert result == datetime(2021, 8, 15)

    def test_missing_month_and_day(self):
        pub_date = {"Year": "2021"}
        result = format_publication_date(pub_date)
        assert result == datetime(2021, 1, 1)

    def test_medline_date_fallback(self):
        pub_date = {"MedlineDate": "2020 Spring"}
        result = format_publication_date(pub_date)
        assert result == datetime(2020, 3, 1)

    def test_empty_returns_none(self):
        assert format_publication_date({}) is None
        assert format_publication_date(None) is None


# ============================================================
# merge_pubmed_icite_pmid_data
# ============================================================

class TestMergePubmedIciteData:
    """Tests for combining PubMed and iCite data."""

    def test_icite_fills_missing_pubmed_values(
            self, sample_pubmed_df, sample_icite_df):
        result = merge_pubmed_icite_pmid_data(
            sample_pubmed_df, sample_icite_df)

        # PubMed title should take precedence when present
        row_a = result[result["pmid"] == 11111111].iloc[0]
        assert row_a["title"] == "Title A"

        # iCite should fill in where PubMed is missing
        row_b = result[result["pmid"] == 22222222].iloc[0]
        assert row_b["authors"] == "Author B icite"

        # iCite metrics should always be present
        assert "citation_count" in result.columns
        assert "relative_citation_ratio" in result.columns

    def test_output_columns_ordered(
            self, sample_pubmed_df, sample_icite_df):
        result = merge_pubmed_icite_pmid_data(
            sample_pubmed_df, sample_icite_df)
        expected = ["pmid", "title", "authors", "publication_date",
                     "citation_count", "relative_citation_ratio"]
        assert list(result.columns) == expected


# ============================================================
# merge_and_clean_project_pmid_info
# ============================================================

class TestMergeAndCleanProjectPmidInfo:
    """Tests for the final merge-and-clean step."""

    def test_removes_pre_2000_publications(self, sample_pmids_df):
        pub_info = pd.DataFrame({
            "pmid": [11111111, 22222222, 33333333],
            "title": ["Title A", "Title B", "Title C"],
            "authors": ["Auth A", "Auth B", "Auth C"],
            "publication_date": pd.to_datetime(
                ["1999-06-15", "2020-03-01", "2021-01-01"]),
            "citation_count": [10, 5, 0],
            "relative_citation_ratio": [1.5, 0.8, 0.0],
        })
        result, removed = merge_and_clean_project_pmid_info(
            sample_pmids_df, pub_info)

        assert 11111111 not in result["pmid"].values
        assert any(removed["reason"] == "Published before 2000")

    def test_removes_erratum_titles(self, sample_pmids_df):
        pub_info = pd.DataFrame({
            "pmid": [11111111, 22222222, 33333333],
            "title": ["Erratum: correction", "Title B", "Title C"],
            "authors": ["Auth A", "Auth B", "Auth C"],
            "publication_date": pd.to_datetime(
                ["2020-06-15", "2020-03-01", "2021-01-01"]),
            "citation_count": [10, 5, 0],
            "relative_citation_ratio": [1.5, 0.8, 0.0],
        })
        result, removed = merge_and_clean_project_pmid_info(
            sample_pmids_df, pub_info)

        assert 11111111 not in result["pmid"].values
        assert any(removed["reason"] == 'Title starts with "Erratum"')


# ============================================================
# NIH RePORTER API — mocked
# ============================================================

class TestGetPmidsFromReporterAPI:
    """Tests for NIH RePORTER PMID gathering with mocked responses."""

    def test_successful_single_page(self):
        mock_response = MagicMock()
        mock_response.status_code = 200
        mock_response.json.return_value = {
            "meta": {"total": 1},
            "results": [
                {"coreproject": "R01CA111111", "pmid": 12345678,
                 "applid": 99999}
            ],
        }

        with patch("modules.gather_publication_data.requests.post",
                    return_value=mock_response):
            results, error = get_pmids_from_nih_reporter_api("R01CA111111")

        assert error is False
        assert len(results) == 1
        assert results[0]["pmid"] == 12345678

    def test_500_error_retries_and_fails(self):
        mock_response = MagicMock()
        mock_response.status_code = 500

        with patch("modules.gather_publication_data.requests.post",
                    return_value=mock_response):
            with patch("modules.gather_publication_data.sleep"):
                results, error = get_pmids_from_nih_reporter_api(
                    "R01CA111111")

        assert error is True
        assert results == []

    def test_request_exception_returns_error(self):
        import requests as req
        with patch("modules.gather_publication_data.requests.post",
                    side_effect=req.exceptions.ConnectionError(
                        "Connection timeout")):
            results, error = get_pmids_from_nih_reporter_api("R01CA111111")

        assert error is True
        assert results == []


# ============================================================
# PubMed / Entrez API — mocked
# ============================================================

class TestGetPubmedInfoMocked:
    """Tests for PubMed info retrieval with mocked Entrez responses."""

    def test_successful_fetch(self):
        mock_record = {
            "PubmedArticle": [{
                "MedlineCitation": {
                    "Article": {
                        "ArticleTitle": "Test Title",
                        "AuthorList": [
                            {"LastName": "Smith", "ForeName": "John"},
                        ],
                        "ArticleDate": [
                            {"Year": "2022", "Month": "05", "Day": "10"}
                        ],
                        "Journal": {
                            "JournalIssue": {"PubDate": {}}
                        },
                    }
                }
            }]
        }

        mock_handle = MagicMock()
        with patch("modules.gather_publication_data.Entrez.efetch",
                    return_value=mock_handle):
            with patch("modules.gather_publication_data.Entrez.read",
                        return_value=mock_record):
                result = get_pubmed_info_from_pmid(12345678)

        assert result is not None
        assert result["title"] == "Test Title"
        assert result["authors"] == "John Smith"
        assert result["publication_date"] == datetime(2022, 5, 10)

    def test_fetch_error_returns_none(self):
        with patch("modules.gather_publication_data.Entrez.efetch",
                    side_effect=Exception("API unavailable")):
            result = get_pubmed_info_from_pmid(99999999)

        assert result is None


# ============================================================
# Environment configuration
# ============================================================

class TestEnvironmentConfig:
    """Tests that required environment variables are configured."""

    def test_ncbi_api_key_is_set(self):
        """Verify NCBI_API_KEY is set. Without a key, PubMed API calls are
        rate-limited to 3/sec instead of 10/sec, tripling the runtime of
        publication gathering."""
        assert os.environ.get("NCBI_API_KEY"), (
            "NCBI_API_KEY is not set. Add it to your .env file. "
            "See README step 5 for instructions.")

    def test_ncbi_email_is_set(self):
        """Verify NCBI_EMAIL is set."""
        assert os.environ.get("NCBI_EMAIL"), (
            "NCBI_EMAIL is not set. Add it to your .env file. "
            "See README step 5 for instructions.")


# ============================================================
# Live API smoke tests — run with: pytest -m live_api
# ============================================================

@pytest.mark.live_api
class TestNIHReporterAPILive:
    """Live smoke tests for the NIH RePORTER Publications API.
    Run with: pytest -m live_api -v
    """

    def test_reporter_api_reachable_and_returns_expected_schema(self):
        """Verify the API is up and the response schema hasn't changed."""
        import requests as req

        response = req.post(
            "https://api.reporter.nih.gov/v2/publications/search",
            json={
                "criteria": {"core_project_nums": ["R01CA263500"]},
                "offset": 0, "limit": 5,
            },
            headers={"accept": "application/json",
                     "Content-Type": "application/json"},
            timeout=30,
        )

        assert response.status_code == 200, (
            f"NIH RePORTER API returned {response.status_code}")

        data = response.json()

        # Verify top-level schema keys
        assert "meta" in data, "Response missing 'meta' key"
        assert "results" in data, "Response missing 'results' key"

        # Verify meta contains expected fields
        meta = data["meta"]
        assert "total" in meta, "Meta missing 'total' field"

        # Verify at least one result with expected fields
        if data["meta"]["total"] > 0:
            result = data["results"][0]
            assert "pmid" in result, "Result missing 'pmid' field"
            assert "coreproject" in result, (
                "Result missing 'coreproject' field")


@pytest.mark.live_api
class TestPubMedEntrezLive:
    """Live smoke tests for PubMed via BioPython Entrez.
    Run with: pytest -m live_api -v
    """

    def test_entrez_efetch_returns_valid_record(self):
        """Verify BioPython Entrez can fetch and parse a known PMID."""
        from Bio import Entrez as ent

        ent.email = os.environ.get("NCBI_EMAIL", "test@example.com")
        ent.api_key = os.environ.get("NCBI_API_KEY", "")

        # PMID 36400004 is a stable, known publication
        handle = ent.efetch(db="pubmed", id="36400004", retmode="xml")
        records = ent.read(handle)
        handle.close()

        assert "PubmedArticle" in records, (
            "Entrez response missing 'PubmedArticle' key")
        assert len(records["PubmedArticle"]) > 0, (
            "Entrez returned empty PubmedArticle list")

        article = (records["PubmedArticle"][0]
                   ["MedlineCitation"]["Article"])
        assert "ArticleTitle" in article, (
            "Article missing 'ArticleTitle' field")
        assert "AuthorList" in article, (
            "Article missing 'AuthorList' field")

    def test_get_pubmed_info_integration(self):
        """Verify our wrapper function works end-to-end with a real PMID."""
        result = get_pubmed_info_from_pmid(36400004)

        assert result is not None, (
            "get_pubmed_info_from_pmid returned None for known PMID")
        assert result["publication_id"] == 36400004
        assert isinstance(result["title"], str) and len(result["title"]) > 0
        assert isinstance(result["authors"], str) and len(result["authors"]) > 0
        assert isinstance(result["publication_date"], datetime)
