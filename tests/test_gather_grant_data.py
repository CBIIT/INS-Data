"""
test_gather_grant_data.py
2024-12-19 ZD

Pytest test suite for the `gather_grant_data.py` module.
"""

import os
import sys
import pandas as pd
import pytest
import requests
from unittest.mock import patch, MagicMock


sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

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



@pytest.fixture
def mock_grant_response():
    """Create a simple mock NIH RePORTER grant result."""
    return {
        "meta":{
            "search_id": "R01CA123456",
            "total": 1,
            "offset": 0,
            "limit": 500,
            "sort_field": "appl_id",
            "sort_order": "desc",
        },
        "results":[
            {
                "project_num": "R01CA123456",
                "principal_investigators": [
                    {"full_name": "John Doe"},
                    {"full_name": "Jane Smith"}
                ],
                "program_officers": [
                    {"full_name": "Alice Johnson"}
                ],
                "agency_funding": [
                    {"code": "CA", "total_cost": 500000},
                    {"code": "Other", "total_cost": 100000}
                ],
                "organization": {
                    "org_name": "Test University",
                    "city": "Test City",
                    "state": "TS"
                },
                "abstract_text": "This is a sample\xad abstract  with\xa0multiple spaces."
            }
        ]
    }


def test_concatenate_full_names():
    """Test concatenation of full names from first and last"""

    # Test simple name list
    test_names = [
        {"full_name": "John Doe"},
        {"full_name": "Jane Smith"}
    ]
    result = concatenate_full_names(test_names)

    assert result == "John Doe; Jane Smith"
    
    # Test blank or missing values
    assert concatenate_full_names([]) == ''
    assert concatenate_full_names(None) == ''


def test_format_name_column():
    """Test formatting capitalization and whitespace of names"""

    assert format_name_column("john  doe") == "John Doe"
    assert format_name_column("JANE DOE") == "Jane Doe"
    assert format_name_column("mary ann smith") == "Mary Ann Smith"


def test_extract_total_cost():
    """Test extraction of NCI funding across award totals."""

    # Test with NCI (CA) funding present
    test_fundings = [
        {"code": "Other_IC", "total_cost": 100000},
        {"code": "CA", "total_cost": 500000},
        {"code": "ABC", "total_cost": 200000}
    ]
    assert extract_total_cost(test_fundings) == 500000

    # Test with no NCI funding
    test_no_nci = [
        {"code": "Other", "total_cost": 100000},
        {"code": "ABC", "total_cost": 200000}
    ]
    assert extract_total_cost(test_no_nci) == 0


def test_clean_abstract():
    """Test text cleaning and normalizing of abstracts."""

    # Test removal of special characters and normalization
    test_cases = [
        ("Example \xad two  spaces", "Example two spaces"),
        ("Example    three spaces", "Example three spaces"),
        ("Example\xa0non-breaking space", "Example non-breaking space"),
        ("Example \nnewline", "Example newline"),
        (None, None),
        ('', ''),
        (123, 123)
    ]

    for input_text, expected in test_cases:
        result = clean_abstract(input_text)
        
        assert result == expected


def test_format_organization_columns():
    """Test handling of organization subfields into single column"""

    # Create a sample DataFrame with nested organization data
    df = pd.DataFrame({
        'organization': [
            {
                'org_name': 'Test University',
                'city': 'Test City',
                'state': 'TS'
            }
        ]
    })

    # Apply the function
    org_subfields = ['org_name', 'city', 'state']
    result_df = format_organization_columns(df, 'organization', org_subfields)

    # Check columns
    assert 'org_name' in result_df.columns
    assert 'city' in result_df.columns
    assert 'state' in result_df.columns
    assert 'organization' not in result_df.columns

    # Check values
    assert result_df['org_name'].iloc[0] == 'Test University'
    assert result_df['city'].iloc[0] == 'Test City'
    assert result_df['state'].iloc[0] == 'TS'
    

@pytest.fixture
def mock_api_response():
    """Create series of minimal mock API responses."""

    def _create_response(status_code, results=None, total=0):
        response = MagicMock()
        response.status_code = status_code
        response.json.return_value = {
            'results': results or [],
            'meta': {'total': total}
        }

        return response
    return _create_response


@pytest.mark.parametrize("search_type,search_value", [
                        ('award', 'R01CA123456'),
                        ('nofo', 'RFA-CA-01-123')])
def test_grants_api_successful_single_page_response(mock_api_response, 
                                                    search_type, search_value):
    """Test successful API call with single page of results."""

    with patch('modules.gather_grant_data.requests.post') as mock_post:
        mock_post.return_value = mock_api_response(
            status_code=200,
            results=[{'project_num': search_value}],
            total=1
        )
        
        grants, failures = get_nih_reporter_grants([search_value], search_type)
        
        assert len(grants) == 1
        assert len(failures) == 0
        assert grants[0]['api_source_search'] == f'{search_type}_{search_value}'


@pytest.mark.parametrize("search_type,search_value", [
                        ('award', 'R01CA123456'),
                        ('nofo', 'RFA-CA-01-123')])
def test_grants_api_pagination(mock_api_response, search_type, search_value):
    """Test handling multiple pages of results."""

    with patch('modules.gather_grant_data.requests.post') as mock_post:
        mock_post.side_effect = [
            mock_api_response(200, [{'id': '1'}] * 500, 750),
            mock_api_response(200, [{'id': '2'}] * 250, 750)
        ]
        
        grants, failures = get_nih_reporter_grants([search_value], search_type)
        
        assert len(grants) == 750
        assert mock_post.call_count == 2


@pytest.mark.parametrize("search_type,search_value", [
                        ('award', 'R01CA123456'),
                        ('nofo', 'RFA-CA-01-123')])
def test_grants_api_non_500_error_handling(mock_api_response, search_type, 
                                           search_value, capsys):
    """Test common error scenarios other than 500 errors."""

    with patch('modules.gather_grant_data.requests.post') as mock_post:
        # Test unspecified error handling
        mock_post.return_value = mock_api_response(123)
        grants, failures = get_nih_reporter_grants([search_value], search_type)

        assert failures[search_value]['failure_type'] == '123 API Error'

    with patch('modules.gather_grant_data.requests.post') as mock_post:
        # Test exception handling (RequestException)
        mock_post.return_value = mock_api_response(200)
        mock_post.side_effect = requests.exceptions.RequestException()
        grants, failures = get_nih_reporter_grants([search_value], search_type)

        assert "An error occurred" in capsys.readouterr().out
        assert len(grants) == 0


@pytest.mark.parametrize("search_type,search_value", [
                        ('award', 'R01CA123456'),
                        ('nofo', 'RFA-CA-01-123')])
def test_grants_api_500_error_retry_success(mock_api_response, search_type, 
                                             search_value, capsys):
    """Test 500 error handling with multiple attempts."""
    
    with (patch('modules.gather_grant_data.requests.post') as mock_post, 
          patch('modules.gather_grant_data.sleep') as mock_sleep):

        mock_sleep.side_effect = lambda x: None
        mock_post.side_effect = [
            mock_api_response(500),
            mock_api_response(200, [{'id': 1}], 1)
        ]
        
        grants, failures = get_nih_reporter_grants([search_value], search_type)
        
        assert len(grants) == 1
        assert len(failures) == 0
        assert mock_post.call_count == 2
        assert mock_sleep.call_count == 1
        assert "Attempt 1/5" in capsys.readouterr().out
        assert "Attempt 2/5" not in capsys.readouterr().out


@pytest.mark.parametrize("search_type,search_value", [
                        ('award', 'R01CA123456'),
                        ('nofo', 'RFA-CA-01-123')])
def test_grants_api_500_error_retry_failure(mock_api_response, search_type, 
                                             search_value, capsys):
    """Test 500 error handling when max attempts reached."""
    
    with (patch('modules.gather_grant_data.requests.post') as mock_post, 
          patch('modules.gather_grant_data.sleep') as mock_sleep):

        mock_sleep.side_effect = lambda x: None
        mock_post.side_effect = [
            mock_api_response(500),
            mock_api_response(500),
            mock_api_response(500),
            mock_api_response(500),
            mock_api_response(500),
            mock_api_response(200, [{'id': 1}], 1) # Should not actually run
        ]
        
        grants, failures = get_nih_reporter_grants([search_value], search_type)
        
        assert len(grants) == 0
        assert len(failures) == 1
        assert "500 Error: Likely too many results" in (
                                         failures[search_value]['failure_type'])
        assert mock_post.call_count == 5
        assert mock_sleep.call_count == 5
        assert "Attempt 5/5" in capsys.readouterr().out


@pytest.mark.parametrize("search_type,search_value", [
                        ('award', ['', 'R01CA123456']),
                        ('nofo', ['','RFA-CA-01-123'])])
def test_grants_api_empty_search_values(mock_api_response, search_type, 
                                        search_value, capsys):
    """Test handling and skipping of empty search values."""

    with patch('modules.gather_grant_data.requests.post') as mock_post:
        mock_post.return_value = mock_api_response(200)
        grants, failures = get_nih_reporter_grants(search_value, search_type, 
                                                   print_meta=True)

        assert f"Blank {search_type} value" in capsys.readouterr().out
        assert mock_post.call_count == 1
        

def test_grants_api_invalid_search_type(mock_api_response):
    """Test handling of search types other than 'award' or 'nofo'."""

    with patch('modules.gather_grant_data.requests.post') as mock_post:
        mock_post.return_value = mock_api_response(200)

        with pytest.raises(ValueError, match="Invalid search type"):
            grants, failures = get_nih_reporter_grants(search_values=['R01CA123456'],
                                                        search_type='wrong_type')
        assert mock_post.call_count == 0


@pytest.mark.parametrize("search_type,search_value", [
                        ('award', 'R01CA123456'),
                        ('nofo', 'RFA-CA-01-123')])
def test_grants_api_blank_search_results(mock_api_response,search_type, 
                                        search_value):
    """Test handling of empty search results from valid search types."""  

    with patch('modules.gather_grant_data.requests.post') as mock_post:
        mock_post.return_value = mock_api_response(200, results=[], total=0)
        grants, failures = get_nih_reporter_grants([search_value], search_type)

        assert len(grants) == 0
        assert len(failures) == 1
        assert failures[search_value]['failure_type'] == "No results found"

