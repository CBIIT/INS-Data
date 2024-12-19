"""
test_gather_grant_data.py
2024-12-19 ZD

Pytest test suite for the `gather_grant_data.py` module.
"""

import os
import sys
import pandas as pd
import pytest
from unittest.mock import patch

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from modules.gather_grant_data import (
    get_nih_reporter_grants,
    concatenate_full_names,
    format_name_column,
    extract_total_cost,
    clean_abstract,
    format_organization_columns,
    clean_grants_data,
    gather_grant_data
)



@pytest.fixture
def mock_grants_json():
    """Create a simple list of mock grant JSON data for testing."""
    return [
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