"""
test_gather_dbgap_data.py
2026-07-06 ZD

Pytest test suite for the gather_dbgap_data.py module.
"""

import os
import sys

import pandas as pd
import pytest
from unittest.mock import patch, MagicMock

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from modules.gather_dbgap_data import (
    get_dbgap_api_data,
    detect_attribute_fieldnames,
    get_principal_investigators,
    get_funding_attributions,
    get_cited_publications,
    get_external_study_urls,
    get_sstr_count,
    get_consent_codes,
    get_assay_method,
    get_dbgap_url,
    join_list,
    clean_dbgap_sstr_metadata,
)


# ============================================================
# Fixtures
# ============================================================

@pytest.fixture
def sample_study_metadata_record():
    """Minimal dbGaP Study Metadata API record."""
    return {
        "full_phs": "phs000001.v1.p1",
        "attribution": [
            {"title": "Principal Investigator", "name": "John Smith"},
            {"title": "Co-Investigator", "name": "Jane Doe"},
            {"title": "Funding Source", "name": "NCI Grant R01CA123456"},
        ],
        "reference": [
            {"pmid": "12345678"},
            {"title": "A Study Title", "authors": "Smith J, Doe J"},
        ],
        "study_url": [
            {"url": "https://example.com/study1"},
            {"url": "https://example.com/study2"},
        ],
    }


@pytest.fixture
def sample_sstr_record():
    """Minimal dbGaP SSTR API record."""
    return {
        "full_phs": "phs000001.v1.p1",
        "study": {
            "consent_groups": [
                {"short_name": "GRU", "subject_count": 100,
                 "sample_count": 200},
                {"short_name": "HMB", "subject_count": 50,
                 "sample_count": 75},
            ]
        },
        "study_stats": {
            "cnt_by_consent_and_sample_use": [
                {"sample_use": "WGS"},
                {"sample_use": "RNA-Seq"},
                {"sample_use": "WGS"},
            ]
        },
    }


# ============================================================
# get_dbgap_url
# ============================================================

class TestGetDbgapUrl:
    """Tests for building dbGaP study page URLs."""

    def test_valid_phs(self):
        url = get_dbgap_url("phs000001.v1.p1")
        assert url == ("https://www.ncbi.nlm.nih.gov/projects/gap/"
                        "cgi-bin/study.cgi?study_id=phs000001")

    def test_short_phs(self):
        url = get_dbgap_url("phs000001")
        assert "phs000001" in url

    def test_invalid_phs_raises(self):
        with pytest.raises(ValueError, match="Invalid phs"):
            get_dbgap_url("GSE12345")


# ============================================================
# join_list
# ============================================================

class TestJoinList:
    """Tests for joining lists into delimited strings."""

    def test_joins_with_semicolons(self):
        assert join_list(["A", "B", "C"]) == "A;B;C"

    def test_custom_separator(self):
        assert join_list(["A", "B"], sep=", ") == "A, B"

    def test_empty_list(self):
        assert join_list([]) == ""


# ============================================================
# get_principal_investigators
# ============================================================

class TestGetPrincipalInvestigators:
    """Tests for extracting PI names from study metadata."""

    def test_finds_pi_by_title(self, sample_study_metadata_record):
        pi_set = {"Principal Investigator"}
        result = get_principal_investigators(
            sample_study_metadata_record, pi_set)
        assert "John Smith" in result

    def test_excludes_non_pi_titles(self, sample_study_metadata_record):
        pi_set = {"Principal Investigator"}
        result = get_principal_investigators(
            sample_study_metadata_record, pi_set)
        assert "Jane Doe" not in result


# ============================================================
# get_cited_publications
# ============================================================

class TestGetCitedPublications:
    """Tests for extracting publication citations."""

    def test_extracts_direct_pmid(self):
        record = {
            "reference": [{"pmid": "99999999"}]
        }
        result = get_cited_publications(record)
        assert "99999999" in result

    def test_handles_missing_pmid(self):
        """When no PMID is available, should try Entrez lookup."""
        record = {
            "reference": [{"title": "Test Title", "authors": "Smith J"}]
        }
        with patch("modules.gather_dbgap_data.get_pmid_from_reference",
                    return_value=None):
            result = get_cited_publications(record)
        assert "No PMID" in result


# ============================================================
# get_external_study_urls
# ============================================================

class TestGetExternalStudyUrls:
    """Tests for extracting study URLs from metadata."""

    def test_joins_multiple_urls(self, sample_study_metadata_record):
        result = get_external_study_urls(sample_study_metadata_record)
        assert "https://example.com/study1" in result
        assert "https://example.com/study2" in result
        assert ";" in result

    def test_empty_urls(self):
        record = {"study_url": []}
        result = get_external_study_urls(record)
        assert result == ""


# ============================================================
# detect_attribute_fieldnames
# ============================================================

class TestDetectAttributeFieldnames:
    """Tests for regex-based detection of attribution field names."""

    def test_detects_pi_titles(self):
        records = [{
            "attribution": [
                {"title": "Principal Investigator"},
                {"title": "Co-PI / Lead Investigator"},
                {"title": "Funding Source"},
            ]
        }]
        result = detect_attribute_fieldnames(records, "pi")
        assert "Principal Investigator" in result
        assert "Funding Source" not in result

    def test_detects_funding_titles(self):
        records = [{
            "attribution": [
                {"title": "Funding Source"},
                {"title": "Grant Number"},
                {"title": "Principal Investigator"},
            ]
        }]
        result = detect_attribute_fieldnames(records, "funding")
        assert "Funding Source" in result
        assert "Grant Number" in result
        assert "Principal Investigator" not in result

    def test_invalid_title_type_raises(self):
        with pytest.raises(ValueError, match="Invalid 'title_type'"):
            detect_attribute_fieldnames([], "invalid")


# ============================================================
# SSTR metadata functions
# ============================================================

class TestGetSstrCount:
    """Tests for summing participant/sample counts across consent groups."""

    def test_participant_count(self, sample_sstr_record):
        assert get_sstr_count(sample_sstr_record, "participant") == 150

    def test_sample_count(self, sample_sstr_record):
        assert get_sstr_count(sample_sstr_record, "sample") == 275

    def test_invalid_count_type_raises(self, sample_sstr_record):
        with pytest.raises(ValueError, match="Invalid 'count_type'"):
            get_sstr_count(sample_sstr_record, "invalid")


class TestGetConsentCodes:
    """Tests for extracting consent codes from SSTR data."""

    def test_extracts_unique_codes(self, sample_sstr_record):
        result = get_consent_codes(sample_sstr_record)
        assert "GRU" in result
        assert "HMB" in result


class TestGetAssayMethod:
    """Tests for extracting assay methods from SSTR data."""

    def test_extracts_unique_methods(self, sample_sstr_record):
        result = get_assay_method(sample_sstr_record)
        # WGS appears twice but should be deduplicated (set)
        assert "WGS" in result
        assert "RNA-Seq" in result


class TestCleanDbgapSstrMetadata:
    """Tests for the full SSTR metadata cleaning function."""

    def test_produces_expected_fields(self, sample_sstr_record):
        result = clean_dbgap_sstr_metadata(sample_sstr_record)
        assert result["full_phs"] == "phs000001.v1.p1"
        assert result["participant_count"] == 150
        assert result["sample_count"] == 275
        assert "GRU" in result["limitations_for_reuse"]
        assert "WGS" in result["assay_method"]


# ============================================================
# get_dbgap_api_data — mocked
# ============================================================

class TestGetDbgapApiDataMocked:
    """Tests for dbGaP API calls with mocked responses."""

    def test_successful_study_metadata_call(self):
        mock_response = MagicMock()
        mock_response.status_code = 200
        mock_response.text = '{"data": {"study_name": "Test Study"}}'

        with patch("modules.gather_dbgap_data.requests.get",
                    return_value=mock_response):
            result = get_dbgap_api_data("phs000001", "study_metadata")

        assert "data" in result
        assert result["data"]["study_name"] == "Test Study"

    def test_502_error_returns_error_record(self):
        mock_response = MagicMock()
        mock_response.status_code = 502

        with patch("modules.gather_dbgap_data.requests.get",
                    return_value=mock_response):
            result = get_dbgap_api_data("phs000001", "study_metadata")

        assert "error" in result
        assert result["error"]["response_code"] == 502

    def test_invalid_api_type_raises(self):
        with pytest.raises(ValueError, match="Invalid 'api_type'"):
            get_dbgap_api_data("phs000001", "invalid_type")


# ============================================================
# Live API smoke tests — run with: pytest -m live_api
# ============================================================

@pytest.mark.live_api
@pytest.mark.xfail(reason="Live dbGaP API — may be slow or unavailable",
                   raises=AssertionError)
class TestDbgapAPILive:
    """Live smoke tests for the dbGaP APIs.
    Skippable offline with: pytest -m "not live_api"
    """

    def test_study_metadata_api_reachable(self):
        """Verify the dbGaP Study Metadata API returns data for a known phs."""
        result = get_dbgap_api_data("phs001115", "study_metadata")

        assert "data" in result, (
            "dbGaP Study Metadata API did not return 'data' key "
            "for known phs001115")
        assert "attribution" in result["data"], (
            "Study Metadata response missing 'attribution' field")

    def test_sstr_summary_api_reachable(self):
        """Verify the dbGaP SSTR Summary API returns data for a known phs."""
        result = get_dbgap_api_data("phs001115", "sstr_summary")

        assert "study" in result, (
            "dbGaP SSTR API did not return 'study' key for known phs001115")
        assert "consent_groups" in result["study"], (
            "SSTR response missing 'consent_groups' field")

    def test_sstr_subjects_api_reachable(self):
        """Verify the dbGaP SSTR Subjects API returns data for a known phs."""
        result = get_dbgap_api_data("phs001115", "sstr_subjects")

        assert "subjects" in result, (
            "dbGaP SSTR Subjects API did not return 'subjects' key "
            "for known phs001115")
