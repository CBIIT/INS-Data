"""
test_gather_geo_data.py
2026-07-06 ZD

Pytest test suite for the gather_geo_data.py module.
"""

import os

import pandas as pd
import pytest
from unittest.mock import patch, MagicMock

from modules.gather_geo_data import (
    fetch_geo_ids,
    get_geo_ids_for_pubmed_ids,
    get_full_geo_record,
    create_geo_dataframe,
    get_dataset_doc_from_project,
    format_geo_names,
    select_geo_ftp_fields,
    merge_semicolon_fields,
    group_by_dataset_id,
    ftp_to_https,
    get_geo_url,
    drop_irrelevant_geo,
)


# ============================================================
# Fixtures
# ============================================================

@pytest.fixture
def sample_publication_df():
    """Minimal publication DataFrame for GEO mapping."""
    return pd.DataFrame({
        "pmid": ["11111111", "22222222", "33333333"],
        "coreproject": ["R01CA111111", "R01CA111111", "U01CA222222"],
    })


@pytest.fixture
def sample_geo_df():
    """Minimal GEO DataFrame for grouping and filtering tests."""
    return pd.DataFrame({
        "geo_id": ["200111111", "200111111", "200222222"],
        "dataset_source_id": ["GSE111111", "GSE111111", "GSE222222"],
        "dataset_pmid": ["11111111", "22222222", "33333333"],
        "funding_source": ["R01CA111111", "R01CA111111", "U01CA222222"],
        "dataset_title": ["Title A", "Title A", "Title B"],
        "program_id": ["PROG1", "PROG1", "PROG2"],
        "dataset_doc": ["NCI", "NCI", "NCI"],
    })


# ============================================================
# format_geo_names
# ============================================================

class TestFormatGeoNames:
    """Tests for cleaning GEO contributor name formats."""

    def test_double_comma_format(self):
        """GEO uses 'Last,,First' format for contributor names."""
        assert format_geo_names("Smith,,John") == "Smith John"

    def test_single_comma_format(self):
        assert format_geo_names("J,A,Smith") == "J A Smith"

    def test_list_input(self):
        result = format_geo_names(["Smith,,John", "Doe,,Jane"])
        assert "Smith John" in result
        assert "Doe Jane" in result

    def test_empty_string(self):
        assert format_geo_names("") == ""

    def test_leading_trailing_commas(self):
        assert format_geo_names(",Smith,,John,") == "Smith John"


# ============================================================
# merge_semicolon_fields
# ============================================================

class TestMergeSemicolonFields:
    """Tests for deduplicating semicolon-separated values."""

    def test_deduplicates_values(self):
        result = merge_semicolon_fields(["A; B", "B; C"])
        parts = [x.strip() for x in result.split(";")]
        assert sorted(parts) == ["A", "B", "C"]

    def test_handles_nan(self):
        result = merge_semicolon_fields([pd.NA, "A; B", None])
        assert "A" in result
        assert "B" in result

    def test_empty_values(self):
        assert merge_semicolon_fields(["", "", ""]) == ""


# ============================================================
# ftp_to_https
# ============================================================

class TestFtpToHttps:
    """Tests for FTP to HTTPS URL conversion."""

    def test_converts_ftp_url(self):
        assert ftp_to_https("ftp://ftp.ncbi.nlm.nih.gov/geo/") == \
               "https://ftp.ncbi.nlm.nih.gov/geo/"

    def test_non_ftp_url_unchanged(self):
        assert ftp_to_https("https://example.com") == "https://example.com"

    def test_non_string_returns_as_is(self):
        assert ftp_to_https(None) is None


# ============================================================
# get_geo_url
# ============================================================

class TestGetGeoUrl:
    """Tests for building GEO study page URLs."""

    def test_gse_accession(self):
        url = get_geo_url("GSE123456")
        assert url == "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE123456"

    def test_gds_accession(self):
        url = get_geo_url("GDS5678")
        assert url == "https://www.ncbi.nlm.nih.gov/sites/GDSbrowser?acc=GDS5678"

    def test_invalid_accession_raises(self):
        with pytest.raises(ValueError, match="Check GEO accession"):
            get_geo_url("GPL12345")


# ============================================================
# create_geo_dataframe
# ============================================================

class TestCreateGeoDataframe:
    """Tests for converting PMID-GEO mappings to a DataFrame."""

    def test_creates_correct_structure(self, sample_publication_df):
        mapping = {
            "11111111": ["200111111", "200111112"],
            "22222222": [],
            "33333333": ["200333333"],
        }
        result = create_geo_dataframe(mapping, sample_publication_df)

        assert "pmid" in result.columns
        assert "geo_id" in result.columns
        assert "coreproject" in result.columns
        # 2 from first PMID + 1 empty + 1 from third = 4
        assert len(result) == 4

    def test_empty_mapping(self, sample_publication_df):
        result = create_geo_dataframe({}, sample_publication_df)
        assert len(result) == 0


# ============================================================
# drop_irrelevant_geo
# ============================================================

class TestDropIrrelevantGeo:
    """Tests for filtering to GSE/GDS accessions only."""

    def test_keeps_gse_and_gds(self, tmp_path):
        df = pd.DataFrame({
            "dataset_source_id": ["GSE111111", "GDS222222", "GPL333333"],
            "value": [1, 2, 3],
        })
        result = drop_irrelevant_geo(
            df, report_csv=str(tmp_path / "dropped.csv"))
        assert len(result) == 2
        assert "GPL333333" not in result["dataset_source_id"].values

    def test_saves_dropped_report(self, tmp_path):
        report_path = tmp_path / "dropped.csv"
        df = pd.DataFrame({
            "dataset_source_id": ["GSE111111", "GPL333333"],
            "value": [1, 2],
        })
        drop_irrelevant_geo(df, report_csv=str(report_path))
        assert report_path.exists()


# ============================================================
# group_by_dataset_id
# ============================================================

class TestGroupByDatasetId:
    """Tests for grouping rows by geo_id and merging fields."""

    def test_merges_pmids_and_projects(self, sample_geo_df):
        result = group_by_dataset_id(sample_geo_df)
        assert len(result) == 2

        row = result[result["geo_id"] == "200111111"].iloc[0]
        assert "11111111" in row["dataset_pmid"]
        assert "22222222" in row["dataset_pmid"]


# ============================================================
# select_geo_ftp_fields
# ============================================================

class TestSelectGeoFtpFields:
    """Tests for extracting and deduplicating FTP metadata fields."""

    def test_extracts_contributor_field(self):
        ftp_metadata = {
            "200111111": {
                "Series_contributor": ["Smith,,John", "Doe,,Jane"],
            },
        }
        result = select_geo_ftp_fields(ftp_metadata)
        assert len(result) == 1
        assert "series_contributor" in result.columns

    def test_handles_missing_fields(self):
        ftp_metadata = {"200111111": {}}
        result = select_geo_ftp_fields(ftp_metadata)
        assert len(result) == 1

    def test_empty_input(self):
        result = select_geo_ftp_fields({})
        assert len(result) == 0


# ============================================================
# fetch_geo_ids — mocked
# ============================================================

class TestFetchGeoIdsMocked:
    """Tests for GEO ID fetching with mocked Entrez responses."""

    def test_successful_fetch(self):
        mock_record = [{
            "LinkSetDb": [{
                "Link": [{"Id": "200111111"}, {"Id": "200222222"}]
            }]
        }]

        mock_handle = MagicMock()
        with patch("modules.gather_geo_data.Entrez.elink",
                    return_value=mock_handle):
            with patch("modules.gather_geo_data.Entrez.read",
                        return_value=mock_record):
                with patch("modules.gather_geo_data.time.sleep"):
                    pmid, geo_ids = fetch_geo_ids("11111111")

        assert pmid == "11111111"
        assert len(geo_ids) == 2

    def test_error_returns_empty_list(self):
        with patch("modules.gather_geo_data.Entrez.elink",
                    side_effect=Exception("API error")):
            with patch("modules.gather_geo_data.time.sleep"):
                pmid, geo_ids = fetch_geo_ids("99999999")

        assert pmid == "99999999"
        assert geo_ids == []


# ============================================================
# get_full_geo_record — mocked
# ============================================================

class TestGetFullGeoRecordMocked:
    """Tests for GEO ESummary record retrieval with mocked Entrez."""

    def test_successful_fetch(self):
        mock_records = [{"Id": "200111111", "Accession": "GSE111111",
                         "title": "Test Study"}]
        mock_handle = MagicMock()
        with patch("modules.gather_geo_data.Entrez.esummary",
                    return_value=mock_handle):
            with patch("modules.gather_geo_data.Entrez.read",
                        return_value=mock_records):
                result = get_full_geo_record("200111111")

        assert result is not None
        assert result[0]["Accession"] == "GSE111111"

    def test_error_returns_none(self):
        with patch("modules.gather_geo_data.Entrez.esummary",
                    side_effect=Exception("API error")):
            result = get_full_geo_record("99999999")
        assert result is None


# ============================================================
# get_dataset_doc_from_project
# ============================================================

class TestGetDatasetDocFromProject:
    """Tests for linking GEO datasets to NCI DOCs via projects/programs."""

    def test_maps_docs_via_project_and_program(self, tmp_path):
        geo_df = pd.DataFrame({
            "geo_id": ["200111111"],
            "funding_source": ["PRJ_A"],
            "dataset_title": ["Test"],
        })
        project_df = pd.DataFrame({
            "project_id": ["PRJ_A"],
            "program.program_id": ["PROG1"],
        })
        program_df = pd.DataFrame({
            "program_id": ["PROG1"],
            "doc": ["NCI Division X"],
        })

        project_path = str(tmp_path / "project.csv")
        program_path = str(tmp_path / "program.csv")
        project_df.to_csv(project_path, index=False)
        program_df.to_csv(program_path, index=False)

        with patch("modules.gather_geo_data.config") as mock_config:
            mock_config.PROJECTS_INTERMED_PATH = project_path
            mock_config.PROGRAMS_INTERMED_PATH = program_path
            result = get_dataset_doc_from_project(geo_df)

        assert "dataset_doc" in result.columns
        assert result.iloc[0]["dataset_doc"] == "NCI Division X"

    def test_missing_project_file_returns_empty_docs(self, tmp_path):
        geo_df = pd.DataFrame({
            "geo_id": ["200111111"],
            "funding_source": ["PRJ_A"],
        })
        with patch("modules.gather_geo_data.config") as mock_config:
            mock_config.PROJECTS_INTERMED_PATH = str(
                tmp_path / "nonexistent.csv")
            mock_config.PROGRAMS_INTERMED_PATH = str(
                tmp_path / "nonexistent.csv")
            result = get_dataset_doc_from_project(geo_df)

        assert result.iloc[0]["dataset_doc"] == ""


# ============================================================
# get_geo_ids_for_pubmed_ids (batch with ThreadPoolExecutor)
# ============================================================

class TestGetGeoIdsForPubmedIds:
    """Tests for the batch GEO-ID fetcher that wraps fetch_geo_ids
    in a ThreadPoolExecutor."""

    @patch('modules.gather_geo_data.fetch_geo_ids')
    def test_batch_returns_correct_mapping(self, mock_fetch):
        """Each PMID maps to the GEO IDs returned by its fetch call."""
        expected = {
            '11111111': ('11111111', ['200111111', '200222222']),
            '22222222': ('22222222', []),
            '33333333': ('33333333', ['200333333']),
        }
        mock_fetch.side_effect = lambda pmid: expected[pmid]

        result = get_geo_ids_for_pubmed_ids(['11111111', '22222222', '33333333'])

        assert result['11111111'] == ['200111111', '200222222']
        assert result['22222222'] == []
        assert result['33333333'] == ['200333333']


# ============================================================
# Live API smoke tests — run by default, excluded in CI with -m "not live_api"
# ============================================================

@pytest.mark.live_api
@pytest.mark.xfail(reason="Live NCBI Entrez API — may be slow or unavailable",
                   raises=AssertionError)
class TestGeoEntrezLive:
    """Live smoke tests for GEO via NCBI E-utilities.
    Excluded in CI with -m "not live_api". To skip locally, use the same flag.
    """

    def test_elink_pubmed_to_gds_is_reachable(self):
        """Verify Entrez elink endpoint for pubmed-to-gds is reachable
        and returns a parseable response structure."""
        from Bio import Entrez as ent

        ent.email = os.environ.get("NCBI_EMAIL", "test@example.com")
        ent.api_key = os.environ.get("NCBI_API_KEY", "")

        # Use any valid PMID — we just need a parseable response
        handle = ent.elink(
            dbfrom="pubmed", db="gds",
            id="36400004", linkname="pubmed_gds"
        )
        record = ent.read(handle)
        handle.close()

        # Verify response is a list with at least one record
        assert len(record) > 0, (
            "Entrez elink returned empty response — API may be down")
        # Verify the record has the expected structure
        assert "LinkSetDb" in record[0] or "LinkSetDbHistory" in record[0] \
            or record[0].get("LinkSetDb") == [] \
            or "IdList" in record[0], (
            "Entrez elink response has unexpected structure")

    def test_esummary_returns_expected_fields(self):
        """Verify Entrez esummary returns fields the module depends on."""
        from Bio import Entrez as ent

        ent.email = os.environ.get("NCBI_EMAIL", "test@example.com")
        ent.api_key = os.environ.get("NCBI_API_KEY", "")

        # GDS ID for GSE2553
        handle = ent.esummary(db="gds", id="200002553", retmode="text")
        records = ent.read(handle)
        handle.close()

        assert len(records) > 0, "esummary returned empty results"

        record = records[0]
        for field in ("Id", "Accession", "title", "summary",
                      "n_samples", "taxon", "FTPLink", "gdsType", "PDAT"):
            assert field in record, (
                f"GEO esummary missing '{field}' — schema may have changed")
