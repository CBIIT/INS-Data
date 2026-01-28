"""
gather_sra_data.py
2026-01-28 ZD

This script defines functions to gather metadata for Sequence Read Archive (SRA)
studies supported by the NCI. The main function accepts a DataFrame containing 
publication information with PMIDs, finds all associated SRA IDs, maps them to 
study-level SRP/ERP IDs, gathers metadata for each study, and outputs a single
enriched dataset CSV for further analysis.

The workflow follows these steps:
1. Read input DataFrame with PubMed IDs
2. Map PMIDs to SRA IDs using NCBI's E-utilities (batched with resumption)
3. Map SRA IDs to SRP/ERP IDs (study-level identifiers)
4. Aggregate batches to create SRP-centric view
5. Gather metadata for each unique SRP/ERP ID using NCBI E-utilities
6. Enrich SRP data with metadata and link to NCI DOCs via projects/programs
7. Output two files:
   - PMID → SRA → SRP/ERP mapping file (for reference)
   - SRP/ERP-centric enriched dataset (main output)

Batch Processing Features:
- Processes PMIDs in configurable batch sizes (default 2000)
- Automatically resumes from last completed batch if interrupted
- Minimal console output for large batch runs
- Intermediate files saved at each step and reused if available

Metadata Gathering:
- Fetches XML metadata from NCBI E-utilities (efetch with db=sra)
- Efficient: Only one API call per SRP (samples first SRA in each study)
- Parses STUDY section for title, abstract, BioProject, PMIDs, etc.
- Rate-limited to 3 requests/second with exponential backoff retry logic
- Returns blank strings for missing fields (never NaN/None/null)

Intermediate files are saved at each step and reused if they already exist.
"""

import os
import sys
import json
import re
import time
import uuid
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Any
from collections import defaultdict
import xml.etree.ElementTree as ET

import pandas as pd
import concurrent.futures as cf
from tqdm import tqdm  # for progress bars
from Bio import Entrez  # for e-Utils API

# Append the project's root directory to the Python path
# This allows for importing config when running as part of main.py or alone
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import config

# Load .env into the process environment so os.environ.get(...) works
# Use an explicit path to the repository root .env for reliability
try:
    from dotenv import load_dotenv
    dotenv_path = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), '.env')
    if os.path.exists(dotenv_path):
        load_dotenv(dotenv_path)
    else:
        # Fall back to default behavior (load from CWD) if explicit file not found
        load_dotenv()
except Exception:
    # If python-dotenv is not available, assume environment variables are
    # already present in the running process (e.g., set by the shell or IDE).
    pass



def get_composite_uuid5(
    df: pd.DataFrame, 
    fields: list, 
    uuid_col: str = 'dataset_uuid'
) -> pd.DataFrame:
    """
    Generate a UUID5 for each row in the DataFrame based on a composite of 
    specified fields. The UUID5 is deterministic and reproducible for the 
    same field values. Raises an error if duplicate UUIDs are generated.

    Args:
        df: Input DataFrame
        fields: List of column names to combine for the UUID5 'name'
        uuid_col: Name of the output column for the UUIDs

    Returns:
        DataFrame with a new column containing the UUID5 values
    """
    # Use a fixed namespace UUID for reproducibility (can be any valid UUID)
    NAMESPACE = uuid.UUID('12345678-1234-5678-1234-567812345678')

    def row_to_uuid(row):
        name = '||'.join(str(row[field]) for field in fields)
        return uuid.uuid5(NAMESPACE, name)

    df[uuid_col] = df.apply(row_to_uuid, axis=1)

    # Check for duplicate UUIDs and raise error if found
    duplicate_uuids = df[uuid_col][df[uuid_col].duplicated()]
    if not duplicate_uuids.empty:
        raise ValueError(
            f"Duplicate UUID5 values found in {uuid_col} column! "
            f"Duplicates:\n{duplicate_uuids.to_list()}\n"
            f"Check your input data for non-unique combinations "
            f"of {fields}."
        )

    return df



def fetch_sra_ids(pmid: str) -> Tuple[str, List[str]]:
    """
    Fetch SRA IDs for a single PubMed ID with retry logic for rate limiting
    
    Args:
        pmid: PubMed ID to query
    
    Returns:
        Tuple of (pmid, list_of_sra_ids)
    """
    
    # Configure Entrez for this thread
    Entrez.email = os.environ.get('NCBI_EMAIL', 'your-email@example.com')
    Entrez.api_key = os.environ.get('NCBI_API_KEY', '')
    Entrez.max_tries = 3
    Entrez.sleep_between_tries = 2
    
    # Retry logic for rate limiting
    max_retries = 5
    retry_delay = 1.0
    
    for attempt in range(max_retries):
        try:
            # Small delay to control API rate
            time.sleep(0.2)
            
            link_handle = Entrez.elink(
                dbfrom="pubmed",
                db="sra",
                id=pmid,
                linkname="pubmed_sra"
            )
            
            link_record = Entrez.read(link_handle)
            link_handle.close()
            
            sra_ids = [
                link['Id']
                for link_set in link_record
                for link in link_set.get('LinkSetDb', [])
                for link in link.get('Link', [])
            ]
            
            return (pmid, sra_ids)
        
        except Exception as e:
            error_msg = str(e)
            # Check if it's a rate limiting error (429)
            if "429" in error_msg or "Too Many Requests" in error_msg:
                if attempt < max_retries - 1:
                    # Exponential backoff: 1s, 2s, 4s, 8s, 16s
                    wait_time = retry_delay * (2 ** attempt)
                    time.sleep(wait_time)
                    continue
                else:
                    print(f"Error processing PMID {pmid}: {e}")
                    return (pmid, [])
            else:
                # For non-rate-limiting errors, fail immediately
                print(f"Error processing PMID {pmid}: {e}")
                return (pmid, [])
    
    return (pmid, [])



def get_sra_ids_for_pubmed_ids(pubmed_ids: List[str]) -> Dict[str, List[str]]:
    """
    Retrieve SRA dataset IDs associated with each PubMed ID in a list of PMIDs. 
    
    Args:
        pubmed_ids: List of PubMed IDs to query
    
    Returns:
        Dictionary mapping PubMed IDs to their associated SRA IDs
    """

    # Configure Entrez
    # Email and api key from hidden local env file. Use default if not defined
    Entrez.email = os.environ.get('NCBI_EMAIL', 'your-email@example.com')
    Entrez.api_key = os.environ.get('NCBI_API_KEY', '')
    if not Entrez.api_key: 
        print(f"WARNING: No NCBI API key in use. Check readme and local .env file."
              f"\nNCBI E-Utilities rate will be limited and may cause errors.")
    
    # Configure API rate limiting and max parallel threads
    Entrez.max_tries = 3
    Entrez.sleep_between_tries = 2
    max_workers = 10
    
    # Get counts for progress tracking
    pmid_count = len(pubmed_ids)
    
    # Use ThreadPoolExecutor with rate-limited concurrency
    # This will run multiple API-calling threads while waiting for responses
    pmid_sra_links = {}
    with cf.ThreadPoolExecutor(max_workers=max_workers) as executor:
        # Submit tasks (PMIDs) to the executor and store futures in dict
        futures = {executor.submit(fetch_sra_ids, pmid): pmid for pmid in pubmed_ids}
        
        # Iterate through futures (PMIDs) as they become available
        for future in tqdm(cf.as_completed(futures),
                          unit="PMID", total=pmid_count, ncols=80,
                          desc="Fetching SRA IDs"):
            pmid, sra_ids = future.result()
            pmid_sra_links[pmid] = sra_ids
    
    return pmid_sra_links



def get_srp_ids_for_sra_id(sra_id: str) -> List[str]:
    """
    Get all SRP/ERP (study-level) IDs associated with a SRA (run-level) ID.

    This parses the XML structure returned by NCBI SRA EFetch and searches for
    Study accessions (SRP numbers from NCBI or ERP numbers from ENA/EBI).

    Args:
        sra_id: String of a Sequence Read Archive run ID (e.g. '31172113')

    Returns:
        List of SRP/ERP accessions e.g. ['SRP123456', 'ERP001234', ...]
        Returns [] if none found or on error.
    """
    # Entrez setup (matches your typical pattern)
    Entrez.email = os.environ.get("NCBI_EMAIL", "your-email@example.com")
    Entrez.api_key = os.environ.get("NCBI_API_KEY", "")
    Entrez.max_tries = 3
    Entrez.sleep_between_tries = 2

    # Retry logic for rate limiting
    max_retries = 5
    retry_delay = 1.0
    
    for attempt in range(max_retries):
        try:
            # EFetch returns full XML for this SRA ID
            handle = Entrez.efetch(db="sra", id=str(sra_id), rettype="xml", retmode="text")
            xml_data = handle.read()
            handle.close()

            if not xml_data.strip():
                return []

            # Use ElementTree (ET) to parse XML for SRP/ERP
            root = ET.fromstring(xml_data)

            study_ids: List[str] = []
            for elem in root.iter():
                # Check attributes for SRP or ERP accessions
                for attr_val in elem.attrib.values():
                    if isinstance(attr_val, str):
                        # Match SRP (NCBI) or ERP (ENA/EBI) study accessions
                        if re.match(r"^(SRP|ERP)\d+$", attr_val.strip()):
                            study_ids.append(attr_val.strip())
                # Occasionally appears in text content too
                if elem.text:
                    if re.match(r"^(SRP|ERP)\d+$", elem.text.strip()):
                        study_ids.append(elem.text.strip())

            # Deduplicate while preserving order
            seen = set()
            out = []
            for s in study_ids:
                if s not in seen:
                    seen.add(s)
                    out.append(s)
            return out

        except Exception as e:
            error_msg = str(e)
            # Check if it's a rate limiting error (429)
            if "429" in error_msg or "Too Many Requests" in error_msg:
                if attempt < max_retries - 1:
                    # Exponential backoff: 1s, 2s, 4s, 8s, 16s
                    wait_time = retry_delay * (2 ** attempt)
                    time.sleep(wait_time)
                    continue
                else:
                    print(f"Error fetching SRP/ERP for SRA ID {sra_id}: {e}")
                    return []
            else:
                # For non-rate-limiting errors, fail immediately
                print(f"Error fetching SRP/ERP for SRA ID {sra_id}: {e}")
                return []
    
    return []



def get_srp_ids_for_sra_ids(
    sra_ids: List[str],
    max_workers: int = 5,
    sleep_between_calls: float = 0.2,
) -> Dict[str, List[str]]:
    """
    Batch map SRA IDs to SRP/ERP accessions using EFetch(db='sra').
    
    Args:
        sra_ids: List of SRA IDs to query
        max_workers: Number of concurrent threads (default 5)
        sleep_between_calls: Delay between API calls in seconds (default 0.2)
    
    Returns:
        Dictionary mapping SRA IDs to their SRP/ERP accessions
        e.g. { '31172113': ['SRP123456'], ... }
    """
    # one-time Entrez init
    Entrez.email = os.environ.get("NCBI_EMAIL", "your-email@example.com")
    Entrez.api_key = os.environ.get("NCBI_API_KEY", "")
    Entrez.max_tries = 3
    Entrez.sleep_between_tries = 2
    if not Entrez.api_key:
        print("WARNING: No NCBI_API_KEY set; requests may be rate limited.")

    sra_ids = [str(u).strip() for u in sra_ids if str(u).strip()]
    result: Dict[str, List[str]] = {}

    def worker(sra_id: str) -> Tuple[str, List[str]]:
        if sleep_between_calls > 0:
            time.sleep(sleep_between_calls)
        srps = get_srp_ids_for_sra_id(sra_id)
        return sra_id, srps

    with cf.ThreadPoolExecutor(max_workers=max_workers) as ex:
        futures = {ex.submit(worker, sra_id): sra_id for sra_id in sra_ids}
        for fut in tqdm(
            cf.as_completed(futures),
            total=len(futures),
            unit="ID",
            ncols=80,
            desc="Fetching SRP/ERP IDs for SRA IDs",
        ):
            uid, srps = fut.result()
            result[uid] = srps

    return result



def create_sra_mapping_dataframe(
    pmid_to_sra_ids: Dict[str, List[str]],
    sra_to_srp_ids: Dict[str, List[str]]
) -> pd.DataFrame:
    """
    Create PMID to SRA to SRP/ERP mapping DataFrame for reference.
    
    Args:
        pmid_to_sra_ids: Dictionary mapping PMIDs to SRA IDs
        sra_to_srp_ids: Dictionary mapping SRA IDs to SRP/ERP IDs
    
    Returns:
        DataFrame with columns: PMID, SRA_IDs, SRP_ERP_IDs
    """
    mapping_rows = []
    
    for pmid, sra_ids in pmid_to_sra_ids.items():
        # Gather all SRP/ERP IDs associated with this PMID's SRA IDs
        srp_ids_set = set()
        for sra_id in sra_ids:
            srp_list = sra_to_srp_ids.get(str(sra_id), [])
            srp_ids_set.update(srp_list)
        
        # Convert to sorted lists for consistent output
        sra_ids_sorted = sorted([str(x) for x in sra_ids])
        srp_ids_sorted = sorted(list(srp_ids_set))
        
        mapping_rows.append({
            'PMID': pmid,
            'SRA_IDs': '; '.join(sra_ids_sorted) if sra_ids_sorted else '',
            'SRP_ERP_IDs': '; '.join(srp_ids_sorted) if srp_ids_sorted else ''
        })
    
    # Create DataFrame with explicit columns even if empty
    if mapping_rows:
        return pd.DataFrame(mapping_rows)
    else:
        return pd.DataFrame(columns=['PMID', 'SRA_IDs', 'SRP_ERP_IDs'])



def create_srp_centric_dataframe(
    pmid_to_sra_ids: Dict[str, List[str]],
    sra_to_srp_ids: Dict[str, List[str]]
) -> pd.DataFrame:
    """
    Create SRP/ERP-centric DataFrame (inverted view).
    
    Args:
        pmid_to_sra_ids: Dictionary mapping PMIDs to SRA IDs
        sra_to_srp_ids: Dictionary mapping SRA IDs to SRP/ERP IDs
    
    Returns:
        DataFrame with columns: SRP_ERP_ID, PMIDs, PMID_Count
    """
    srp_to_pmids = defaultdict(set)
    
    for pmid, sra_ids in pmid_to_sra_ids.items():
        for sra_id in sra_ids:
            srp_list = sra_to_srp_ids.get(str(sra_id), [])
            for srp_id in srp_list:
                if srp_id:
                    srp_to_pmids[srp_id].add(pmid)
    
    # Create SRP-centric rows
    srp_rows = []
    for srp_id in sorted(srp_to_pmids.keys()):
        pmids = sorted(list(srp_to_pmids[srp_id]))
        srp_rows.append({
            'SRP_ERP_ID': srp_id,
            'PMIDs': '; '.join(pmids),
            'PMID_Count': len(pmids)
        })
    
    # Create DataFrame with explicit columns even if empty
    if srp_rows:
        return pd.DataFrame(srp_rows)
    else:
        return pd.DataFrame(columns=['SRP_ERP_ID', 'PMIDs', 'PMID_Count'])



def merge_semicolon_fields(values):
    """Deduplicate values within semicolon-separated list-like strings."""

    # Create an empty set to collect all unique items
    unique_items = set()
    
    # Process each string value
    for val in values:
        if pd.notna(val) and val != '':
            # Split by semicolon and add each item to set
            items = [item.strip() for item in str(val).split(';')]
            unique_items.update(item for item in items if item)
    
    # Join the unique items back with semicolons
    return '; '.join(sorted(unique_items)) if unique_items else ''



def aggregate_batch_mappings(mapping_dfs: List[pd.DataFrame]) -> pd.DataFrame:
    """
    Aggregate mapping DataFrames from multiple batches.
    
    Args:
        mapping_dfs: List of mapping DataFrames from batches
    
    Returns:
        Combined mapping DataFrame
    """
    return pd.concat(mapping_dfs, ignore_index=True)



def aggregate_batch_srp_data(srp_dfs: List[pd.DataFrame]) -> pd.DataFrame:
    """
    Aggregate SRP-centric DataFrames from multiple batches.
    Combines PMIDs for each unique SRP/ERP ID.
    
    Args:
        srp_dfs: List of SRP-centric DataFrames from batches
    
    Returns:
        Combined SRP-centric DataFrame with merged PMIDs
    """
    srp_pmid_map = defaultdict(set)
    
    for batch_df in srp_dfs:
        for _, row in batch_df.iterrows():
            srp_id = row['SRP_ERP_ID']
            # Convert PMIDs to string to handle both string and int types
            pmid_value = row['PMIDs']
            if pd.notna(pmid_value):
                # Convert to string first if it's an integer
                pmid_str = str(pmid_value) if not isinstance(pmid_value, str) else pmid_value
                pmids = pmid_str.split('; ') if pmid_str else []
            else:
                pmids = []
            srp_pmid_map[srp_id].update(pmids)
    
    # Create final SRP-centric rows
    final_srp_rows = []
    for srp_id in sorted(srp_pmid_map.keys()):
        pmids = sorted(list(srp_pmid_map[srp_id]))
        final_srp_rows.append({
            'SRP_ERP_ID': srp_id,
            'PMIDs': '; '.join(pmids),
            'PMID_Count': len(pmids)
        })
    
    return pd.DataFrame(final_srp_rows)



def get_max_batch_number(batch_dir: str) -> int:
    """Get the maximum batch number from existing batch filenames."""
    
    if not os.path.exists(batch_dir):
        return 0
    
    existing_batch_files = [
        file for file in os.listdir(batch_dir) 
        if file.startswith('batch_') and file.endswith('_mapping.csv')
    ]
    
    if not existing_batch_files:
        return 0
    
    # Extract batch numbers from filenames and return the maximum
    batch_numbers = []
    for file in existing_batch_files:
        match = re.search(r'batch_(\d+)_', file)
        if match:
            batch_numbers.append(int(match.group(1)))
    
    return max(batch_numbers) if batch_numbers else 0



def load_all_batch_files(batch_dir: str, file_suffix: str) -> pd.DataFrame:
    """Load all batch files with a given suffix into a single DataFrame."""
    
    if not os.path.exists(batch_dir):
        return pd.DataFrame()
    
    filenames = [
        file for file in os.listdir(batch_dir) 
        if file.endswith(file_suffix)
    ]
    
    if not filenames:
        return pd.DataFrame()
    
    filepaths = [os.path.join(batch_dir, file) for file in filenames]
    df = pd.concat(map(pd.read_csv, filepaths), ignore_index=True)
    
    return df



def fetch_sra_metadata_for_srp(
    srp_id: str, 
    sra_id_sample: str
) -> Dict[str, Any]:
    """
    Fetch metadata for a single SRP/ERP study using one sample SRA ID.
    
    Fetches XML metadata from NCBI E-utilities and parses the STUDY section.
    The STUDY section is identical across all SRAs within the same SRP, so
    sampling one SRA per SRP is sufficient.
    
    Args:
        srp_id: SRP/ERP study ID (e.g., 'SRP472047')
        sra_id_sample: One SRA ID from this study (used to fetch metadata)
    
    Returns:
        Dictionary with metadata fields (blank strings for missing values)
    """
    # Configure Entrez for this thread
    Entrez.email = os.environ.get('NCBI_EMAIL', 'your-email@example.com')
    Entrez.api_key = os.environ.get('NCBI_API_KEY', '')
    Entrez.max_tries = 3
    Entrez.sleep_between_tries = 2
    
    # Rate limiting: max 3 requests/second
    time.sleep(0.34)
    
    # Retry logic for rate limiting
    max_retries = 5
    retry_delay = 1.0
    
    for attempt in range(max_retries):
        try:
            # Fetch XML metadata for the sample SRA ID
            handle = Entrez.efetch(
                db="sra", 
                id=str(sra_id_sample), 
                rettype="xml", 
                retmode="text"
            )
            xml_data = handle.read()
            handle.close()
            
            if not xml_data.strip():
                return _create_blank_metadata_dict(srp_id)
            
            # Parse XML
            root = ET.fromstring(xml_data)
            
            # Initialize metadata dictionary with blank defaults
            metadata = {
                'type': 'dataset',
                'dataset_uuid': '',
                'dataset_source_repo': 'SRA',
                'dataset_title': '',
                'description': '',
                'dataset_source_id': srp_id,
                'dataset_source_url': (
                    f'https://trace.ncbi.nlm.nih.gov/Traces/?study={srp_id}'
                ),
                'PI_name': '',
                'GPA': '',
                'dataset_doc': '',
                'dataset_pmid': '',
                'funding_source': '',
                'release_date': '',
                'limitations_for_reuse': '',
                'assay_method': '',
                'study_type': 'Genomic sequencing',
                'primary_disease': 'Not Reported',
                'participant_count': '',
                'sample_count': '',
                'study_links': '',
                'related_genes': '',
                'related_diseases': '',
                'related_terms': ''
            }
            
            # Extract STUDY section fields
            study_elem = root.find('.//STUDY')
            if study_elem is not None:
                # Study title
                title_elem = study_elem.find('.//STUDY_TITLE')
                if title_elem is not None and title_elem.text:
                    metadata['dataset_title'] = title_elem.text.strip()
                
                # Study abstract
                abstract_elem = study_elem.find('.//STUDY_ABSTRACT')
                if abstract_elem is not None and abstract_elem.text:
                    metadata['description'] = abstract_elem.text.strip()
                
                # Center name (sometimes contains PI information)
                center_name = study_elem.get('center_name', '')
                if center_name:
                    # Filter out "GEO Curators" and "GEO"
                    if center_name not in ['GEO Curators', 'GEO']:
                        metadata['PI_name'] = center_name.strip()
                
                # Study type
                study_type_elem = study_elem.find(
                    './/STUDY_TYPE'
                )
                if study_type_elem is not None:
                    existing_type = study_type_elem.get(
                        'existing_study_type', 
                        ''
                    )
                    if existing_type:
                        metadata['assay_method'] = existing_type.strip()
                
                # BioProject ID
                bioproject_elem = study_elem.find(
                    ".//EXTERNAL_ID[@namespace='BioProject']"
                )
                bioproject_id = ''
                if bioproject_elem is not None and bioproject_elem.text:
                    bioproject_id = bioproject_elem.text.strip()
                
                # GEO ID
                geo_elem = study_elem.find(
                    ".//EXTERNAL_ID[@namespace='GEO']"
                )
                geo_id = ''
                if geo_elem is not None and geo_elem.text:
                    geo_id = geo_elem.text.strip()
                
                # Build study_links
                links = []
                if geo_id:
                    geo_url = (
                        f'https://www.ncbi.nlm.nih.gov/geo/query/'
                        f'acc.cgi?acc={geo_id}'
                    )
                    links.append(geo_url)
                if bioproject_id:
                    bioproject_url = (
                        f'https://www.ncbi.nlm.nih.gov/bioproject/'
                        f'{bioproject_id}'
                    )
                    links.append(bioproject_url)
                metadata['study_links'] = '; '.join(links)
                
                # PubMed IDs from STUDY_LINKS
                pmids = []
                for study_link in study_elem.findall('.//STUDY_LINK'):
                    xref_link = study_link.find('.//XREF_LINK')
                    if xref_link is not None:
                        db_elem = xref_link.find('.//DB')
                        id_elem = xref_link.find('.//ID')
                        if (db_elem is not None and 
                            db_elem.text and 
                            db_elem.text.strip().lower() == 'pubmed' and
                            id_elem is not None and 
                            id_elem.text):
                            pmids.append(id_elem.text.strip())
                if pmids:
                    metadata['dataset_pmid'] = '; '.join(pmids)
            
            # Try to get PI name from Organization/Contact
            # Collect all collaborators (semicolon-separated list)
            collaborators = []
            for contact in root.findall('.//Organization/Contact'):
                name_elem = contact.find('.//Name')
                if name_elem is not None:
                    first = name_elem.find('.//First')
                    last = name_elem.find('.//Last')
                    if first is not None and last is not None:
                        first_text = first.text.strip() if first.text else ''
                        last_text = last.text.strip() if last.text else ''
                        if first_text and last_text:
                            full_name = f'{first_text} {last_text}'
                            # Filter out "GEO Curators" and "GEO"
                            if full_name not in ['GEO Curators', 'GEO']:
                                collaborators.append(full_name)
            
            # If we found collaborators, use them (semicolon-separated)
            # Otherwise keep the center_name if it was set
            if collaborators:
                metadata['PI_name'] = '; '.join(collaborators)
            
            # Release date from first RUN published date
            run_elem = root.find('.//RUN')
            if run_elem is not None:
                published = run_elem.get('published', '')
                if published:
                    # Format as yyyy-mm-dd (remove time if present)
                    date_part = published.strip().split()[0]
                    metadata['release_date'] = date_part
            
            # Library strategy (e.g., RNA-Seq)
            library_strategy_elem = root.find('.//LIBRARY_STRATEGY')
            if library_strategy_elem is not None and library_strategy_elem.text:
                # Only override if assay_method not already set from study type
                if not metadata['assay_method']:
                    metadata['assay_method'] = (
                        library_strategy_elem.text.strip()
                    )
            
            # NOTE: Sample count is unreliable when fetching from a single 
            # SRA per SRP. Leaving blank for now.
            
            return metadata
            
        except Exception as e:
            error_msg = str(e)
            # Check if it's a rate limiting error (429)
            if "429" in error_msg or "Too Many Requests" in error_msg:
                if attempt < max_retries - 1:
                    # Exponential backoff
                    wait_time = retry_delay * (2 ** attempt)
                    time.sleep(wait_time)
                    continue
                else:
                    print(
                        f"Error fetching metadata for SRP {srp_id}: {e}"
                    )
                    return _create_blank_metadata_dict(srp_id)
            else:
                # For non-rate-limiting errors, fail immediately
                print(f"Error fetching metadata for SRP {srp_id}: {e}")
                return _create_blank_metadata_dict(srp_id)
    
    return _create_blank_metadata_dict(srp_id)


def _create_blank_metadata_dict(srp_id: str) -> Dict[str, Any]:
    """
    Create a blank metadata dictionary with all required fields.
    
    Args:
        srp_id: SRP/ERP study ID
    
    Returns:
        Dictionary with all fields set to appropriate blank values
    """
    return {
        'type': 'dataset',
        'dataset_uuid': '',
        'dataset_source_repo': 'SRA',
        'dataset_title': '',
        'description': '',
        'dataset_source_id': srp_id,
        'dataset_source_url': (
            f'https://trace.ncbi.nlm.nih.gov/Traces/?study={srp_id}'
        ),
        'PI_name': '',
        'GPA': '',
        'dataset_doc': '',
        'dataset_pmid': '',
        'funding_source': '',
        'release_date': '',
        'limitations_for_reuse': '',
        'assay_method': '',
        'study_type': 'Genomic sequencing',
        'primary_disease': 'Not Reported',
        'participant_count': '',
        'sample_count': '',
        'study_links': '',
        'related_genes': '',
        'related_diseases': '',
        'related_terms': ''
    }


def create_srp_to_sra_mapping(
    pmid_to_sra_ids: Dict[str, List[str]],
    sra_to_srp_ids: Dict[str, List[str]]
) -> Dict[str, str]:
    """
    Create mapping of SRP/ERP IDs to one sample SRA ID per study.
    
    For each SRP, selects the first SRA ID found in the mapping. This is
    used to efficiently fetch metadata (one API call per SRP instead of
    one per SRA).
    
    Args:
        pmid_to_sra_ids: Dictionary mapping PMIDs to SRA IDs
        sra_to_srp_ids: Dictionary mapping SRA IDs to SRP/ERP IDs
    
    Returns:
        Dictionary mapping SRP/ERP IDs to one sample SRA ID
    """
    srp_to_sample_sra = {}
    
    for pmid, sra_ids in pmid_to_sra_ids.items():
        for sra_id in sra_ids:
            srp_list = sra_to_srp_ids.get(str(sra_id), [])
            for srp_id in srp_list:
                if srp_id and srp_id not in srp_to_sample_sra:
                    srp_to_sample_sra[srp_id] = str(sra_id)
    
    return srp_to_sample_sra


def gather_srp_metadata(
    srp_ids: List[str],
    srp_to_sra_mapping: Dict[str, str],
    max_workers: int = 3
) -> pd.DataFrame:
    """
    Gather metadata for SRP/ERP IDs from NCBI SRA database.
    
    Uses NCBI E-utilities to fetch XML metadata for each SRP/ERP study.
    For efficiency, fetches metadata from only one SRA ID per SRP (the
    STUDY section is identical across all SRAs in the same SRP).
    
    Args:
        srp_ids: List of SRP/ERP study IDs
        srp_to_sra_mapping: Dictionary mapping SRP IDs to sample SRA IDs
        max_workers: Number of concurrent threads (default 3 for rate limit)
    
    Returns:
        DataFrame with SRP_ERP_ID and metadata columns
    """
    # Configure Entrez
    Entrez.email = os.environ.get('NCBI_EMAIL', 'your-email@example.com')
    Entrez.api_key = os.environ.get('NCBI_API_KEY', '')
    if not Entrez.api_key:
        print(
            "WARNING: No NCBI API key in use. "
            "Check readme and local .env file."
        )
    
    # Build list of (srp_id, sra_id_sample) tuples
    srp_sra_pairs = []
    for srp_id in srp_ids:
        sra_sample = srp_to_sra_mapping.get(srp_id, '')
        if sra_sample:
            srp_sra_pairs.append((srp_id, sra_sample))
        else:
            # No SRA sample found - add blank metadata
            srp_sra_pairs.append((srp_id, ''))
    
    # Fetch metadata for each SRP using thread pool
    metadata_records = []
    
    def worker(srp_sra_pair):
        srp_id, sra_sample = srp_sra_pair
        if sra_sample:
            return fetch_sra_metadata_for_srp(srp_id, sra_sample)
        else:
            return _create_blank_metadata_dict(srp_id)
    
    with cf.ThreadPoolExecutor(max_workers=max_workers) as executor:
        futures = {
            executor.submit(worker, pair): pair 
            for pair in srp_sra_pairs
        }
        
        for future in tqdm(
            cf.as_completed(futures),
            total=len(futures),
            unit="SRP",
            ncols=80,
            desc="Fetching SRP metadata"
        ):
            metadata = future.result()
            metadata_records.append(metadata)
    
    # Convert to DataFrame
    metadata_df = pd.DataFrame(metadata_records)
    
    # Ensure dataset_source_id is set correctly
    # (it's used as the join key with SRP_ERP_ID)
    if 'dataset_source_id' in metadata_df.columns:
        metadata_df['SRP_ERP_ID'] = metadata_df['dataset_source_id']
    
    return metadata_df



def enrich_srp_data_with_metadata(
    srp_df: pd.DataFrame,
    metadata_df: pd.DataFrame
) -> pd.DataFrame:
    """
    Merge SRP-centric data with metadata.
    
    Args:
        srp_df: SRP-centric DataFrame with SRP_ERP_ID, PMIDs, PMID_Count
        metadata_df: Metadata DataFrame with SRP_ERP_ID and metadata columns
    
    Returns:
        Enriched DataFrame combining both sources
    """
    # Merge on SRP_ERP_ID
    enriched_df = pd.merge(
        metadata_df,
        srp_df[['SRP_ERP_ID', 'PMIDs', 'PMID_Count']], 
        on='SRP_ERP_ID', 
        how='left'
    )
    
    # Merge PMIDs from the mapping with PMIDs from metadata
    # (metadata may have PMIDs from the SRA record itself)
    def merge_pmids(row):
        """Combine PMIDs from mapping and metadata, deduplicate."""
        mapping_pmids = str(row.get('PMIDs', '')) if pd.notna(row.get('PMIDs')) else ''
        metadata_pmids = str(row.get('dataset_pmid', '')) if pd.notna(row.get('dataset_pmid')) else ''
        
        all_pmids = set()
        if mapping_pmids:
            all_pmids.update(p.strip() for p in mapping_pmids.split(';') if p.strip())
        if metadata_pmids:
            all_pmids.update(p.strip() for p in metadata_pmids.split(';') if p.strip())
        
        return '; '.join(sorted(all_pmids)) if all_pmids else ''
    
    # Update dataset_pmid with merged PMIDs
    enriched_df['dataset_pmid'] = enriched_df.apply(merge_pmids, axis=1)
    
    # Drop the PMIDs column from the original mapping (we have dataset_pmid now)
    if 'PMIDs' in enriched_df.columns:
        enriched_df = enriched_df.drop(columns=['PMIDs'])
    
    # Rename dataset_source_id to match expected output format if needed
    # (It should already be set to SRP_ERP_ID in the metadata)
    
    return enriched_df


def get_dataset_doc_from_project(df_sra: pd.DataFrame) -> pd.DataFrame:
    """
    Merge dataset df with project and program dfs to link DOCs with datasets.
    
    Links SRA datasets to NCI DOCs via the publications -> projects -> 
    programs chain. Uses dataset_pmid field to find matching publications
    and their associated projects.
    
    Args:
        df_sra: DataFrame with SRA dataset metadata
    
    Returns:
        DataFrame with program_id and dataset_doc columns added
    """
    print(f"---\nLinking SRA datasets with NCI DOCs via programs and projects...")

    # Create copy with empty fields as backup
    df_null_result = df_sra.copy()
    df_null_result['program_id'] = ''
    df_null_result['dataset_doc'] = ''
    df_null_result['funding_source'] = ''

    # Check and load projects and program intermediate CSVs
    project_path = config.PROJECTS_INTERMED_PATH
    program_path = config.PROGRAMS_INTERMED_PATH
    publication_path = config.PUBLICATIONS_INTERMED_PATH

    # Check if publication CSV file exists
    if os.path.exists(publication_path):
        df_publication = pd.read_csv(publication_path, dtype={'pmid': str})
        print(f"Loaded publications from {publication_path}")
    else:
        print("Warning: No publication file available. DOCs not mapped.")
        return df_null_result

    # Check if project CSV file exists
    if os.path.exists(project_path):
        df_project = pd.read_csv(project_path)
        print(f"Loaded projects from {project_path}")
    else: 
        print("Warning: No project file available. DOCs not mapped.")
        return df_null_result
    
    # Check if program CSV file exists
    if os.path.exists(program_path):
        df_program = pd.read_csv(program_path)
        print(f"Loaded programs from {program_path}")
    else: 
        print("Warning: No program file available. DOCs not mapped.")
        return df_null_result
    
    # Check if right columns exist in publication CSV
    if not all(
        col in df_publication.columns for col in ['pmid', 'coreproject']):
        print(f"Warning: Publication CSV missing required columns. DOCs not mapped.")
        return df_null_result
    
    # Check if right columns exist in project CSV
    if not all(
        col in df_project.columns for col in ['project_id', 'program.program_id']):
        print(f"Warning: Project CSV missing required columns. DOCs not mapped.")
        return df_null_result
    
    # Check if right columns exist in program CSV
    if not all(
        col in df_program.columns for col in ['program_id', 'doc']):
        print("Warning: Program CSV missing required columns. DOCs not mapped.")
        return df_null_result

    # Prepare publication data: create PMID to project mapping
    df_publication['pmid'] = df_publication['pmid'].astype(str)
    pmid_to_project = dict(zip(
        df_publication['pmid'], 
        df_publication['coreproject']
    ))

    # For each SRA dataset, extract PMIDs and find associated projects
    def get_projects_for_dataset(row):
        """Extract all unique projects from semicolon-separated PMIDs."""
        pmid_str = str(row.get('dataset_pmid', ''))
        if not pmid_str or pmid_str == 'nan':
            return ''
        
        pmids = [p.strip() for p in pmid_str.split(';') if p.strip()]
        projects = set()
        for pmid in pmids:
            project = pmid_to_project.get(pmid, '')
            if project:
                projects.add(project)
        
        return '; '.join(sorted(projects)) if projects else ''
    
    df_sra['funding_source'] = df_sra.apply(get_projects_for_dataset, axis=1)

    # Get programs for each SRA via projects
    # First, create a lookup for project to program
    project_to_program = dict(zip(
        df_project['project_id'],
        df_project['program.program_id']
    ))
    
    def get_programs_for_dataset(row):
        """Extract all unique programs from semicolon-separated projects."""
        project_str = str(row.get('funding_source', ''))
        if not project_str or project_str == 'nan':
            return ''
        
        projects = [p.strip() for p in project_str.split(';') if p.strip()]
        programs = set()
        for project in projects:
            program = project_to_program.get(project, '')
            if program:
                programs.add(program)
        
        return '; '.join(sorted(programs)) if programs else ''
    
    df_sra['program_id'] = df_sra.apply(get_programs_for_dataset, axis=1)

    # Get DOCs for each SRA via programs
    program_to_doc = dict(zip(
        df_program['program_id'],
        df_program['doc']
    ))
    
    def get_docs_for_dataset(row):
        """Extract all unique DOCs from semicolon-separated programs."""
        program_str = str(row.get('program_id', ''))
        if not program_str or program_str == 'nan':
            return ''
        
        programs = [p.strip() for p in program_str.split(';') if p.strip()]
        docs = set()
        for program in programs:
            doc = program_to_doc.get(program, '')
            if doc and doc != 'nan':
                docs.add(str(doc))
        
        return '; '.join(sorted(docs)) if docs else ''
    
    df_sra['dataset_doc'] = df_sra.apply(get_docs_for_dataset, axis=1)

    print(f"Done! NCI DOCs mapped to SRA datasets.\n"
          f"Unique DOC combinations:        {df_sra['dataset_doc'].nunique()}\n"
          f"Unique program combinations:    {df_sra['program_id'].nunique()}\n"
          f"Unique project combinations:    {df_sra['funding_source'].nunique()}")

    return df_sra




def gather_sra_data(
    publications_df: pd.DataFrame,
    overwrite_intermeds: bool = False,
    batch_size: int = 2000
) -> pd.DataFrame:
    """
    Main function for the SRA data gathering workflow. Uses NCBI E-Utilities to
    pull information for all SRA datasets associated with provided publications.
    
    This function performs the following steps:
    1. Map PMIDs to SRA IDs (run-level) 
    2. Map SRA IDs to SRP/ERP IDs (study-level)
    3. Create SRP-centric view with associated PMIDs
    4. Gather metadata for each unique SRP/ERP ID
    5. Enrich SRP data with metadata
    6. Export mapping CSV (for reference) and enriched SRP datasets CSV
    
    Input and output filepaths are defined in config.py

    Args:
        publications_df: DataFrame output of publication gathering workflow
        overwrite_intermeds: Rerun all steps, even if intermediate files 
            already exist (default False). False will save time and avoid
            unnecessary repeated API calls.
        batch_size: Number of PMIDs to process per batch (default 2000)
        
    Returns:
        pd.DataFrame: SRP/ERP-centric DataFrame enriched with metadata
            Columns: SRP_ERP_ID, PMIDs, PMID_Count, [metadata columns]
    """

    print(f"\n---\nDATASETS - SRA\nGathering, formatting, and saving SRA data...")

    # Define paths from config
    sra_mapping_path = config.SRA_PMID_MAPPING_PATH
    srp_centric_path = config.SRA_SRP_CENTRIC_PATH
    sra_datasets_path = config.SRA_INTERMED_PATH
    batch_dir = config.SRA_BATCH_DIR

    # Create directories if they don't exist
    for path in [sra_mapping_path, srp_centric_path, sra_datasets_path]:
        Path(path).parent.mkdir(parents=True, exist_ok=True)
    Path(batch_dir).mkdir(parents=True, exist_ok=True)
    
    # Validate required columns
    required_cols = ['pmid']
    missing_cols = [col for col in required_cols if col not in publications_df.columns]
    if missing_cols:
        raise ValueError(f"Missing required columns: {missing_cols}")
    
    # Convert PMIDs to strings for consistency and get unique values
    publications_df['pmid'] = publications_df['pmid'].astype(str)
    unique_pmids = publications_df['pmid'].dropna().unique().tolist()
    
    print(f"\nLoaded {len(publications_df)} publications with {len(unique_pmids)} unique PMIDs")

    # Check if SRP-centric aggregation already exists and we're not overwriting
    # This is the key file that indicates all batch processing is complete
    if os.path.exists(srp_centric_path) and not overwrite_intermeds:
        print(f"\nSRP-centric aggregation file already exists:")
        print(f"  {srp_centric_path}")
        print(f"Skipping batch processing and loading aggregated data...")
        
        # Load the aggregated files
        final_mapping_df = pd.read_csv(sra_mapping_path)
        final_srp_df = pd.read_csv(srp_centric_path)
        
        print(f"  Loaded {len(final_srp_df)} unique SRP/ERP IDs")
        
        # Jump to metadata gathering (Step 5)
        # Continue with the rest of the workflow from metadata gathering onwards
        skip_to_metadata = True
    else:
        skip_to_metadata = False
        
        # If overwriting, clean up existing aggregated files
        if overwrite_intermeds:
            if os.path.exists(sra_mapping_path):
                os.remove(sra_mapping_path)
            if os.path.exists(srp_centric_path):
                os.remove(srp_centric_path)

    # STEP 1-3: Process PMIDs in batches to build PMID to SRA to SRP mapping
    
    if not skip_to_metadata:
        print(f"\nStep 1-3: Mapping PMIDs to SRA datasets")
        
        total_pmids = len(unique_pmids)
        num_batches = (total_pmids + batch_size - 1) // batch_size

        print(f"  Total PMIDs: {total_pmids}")
        print(f"  Batch size: {batch_size}")
        print(f"  Total batches: {num_batches}")
        
        # Check for existing batches to enable resumption
        max_existing_batch = get_max_batch_number(batch_dir)
        
        if max_existing_batch > 0 and not overwrite_intermeds:
            print(f"  Found {max_existing_batch} existing batch file(s)")
            print(f"  Resuming from batch {max_existing_batch + 1}")
        elif overwrite_intermeds:
            print(f"  Overwriting existing intermediate files")
        
        # Process batches
        for batch_num in range(num_batches):
            batch_number = batch_num + 1
            start_idx = batch_num * batch_size
            end_idx = min((batch_num + 1) * batch_size, total_pmids)
            batch_pmids = unique_pmids[start_idx:end_idx]
            
            # Define output paths for this batch
            mapping_path_batch = os.path.join(
                batch_dir, 
                f'batch_{batch_number:04d}_mapping.csv'
            )
            srp_path_batch = os.path.join(
                batch_dir,
                f'batch_{batch_number:04d}_srp.csv'
            )
            
            # Check if batch already processed
            if (os.path.exists(mapping_path_batch) and 
                os.path.exists(srp_path_batch) and 
                not overwrite_intermeds):
                # Skip with minimal output for cleaner logs
                if batch_number % 10 == 0 or batch_number == num_batches:
                    print(f"  Batch {batch_number}/{num_batches}: Loaded from cache (PMIDs {start_idx}-{end_idx-1})")
                continue
            
            # Process new batch
            print(f"\n  Processing batch {batch_number}/{num_batches}")
            print(f"    PMIDs {start_idx} to {end_idx - 1} ({len(batch_pmids)} PMIDs)")
            
            # Step 1: PMID → SRA IDs
            print(f"    [1/3] Fetching SRA IDs...")
            pmid_to_sra_ids = get_sra_ids_for_pubmed_ids(batch_pmids)
            
            pmids_with_sra = len([p for p, sra_list in pmid_to_sra_ids.items() if sra_list])
            print(f"          Found SRA IDs for {pmids_with_sra}/{len(batch_pmids)} PMIDs")
            
            # Step 2: SRA IDs → SRP/ERP IDs
            all_sra_ids = [item for sublist in pmid_to_sra_ids.values() for item in sublist]
            unique_sra_ids = list(set(all_sra_ids))
            
            print(f"    [2/3] Fetching SRP/ERP IDs for {len(unique_sra_ids)} unique SRA IDs...")
            
            if unique_sra_ids:
                sra_to_srp_ids = get_srp_ids_for_sra_ids(unique_sra_ids)
            else:
                sra_to_srp_ids = {}
            
            sra_with_srp = len([sra for sra, srp_list in sra_to_srp_ids.items() if srp_list])
            print(f"          Found SRP/ERP IDs for {sra_with_srp}/{len(unique_sra_ids)} SRA IDs")
            
            # Step 3: Build DataFrames
            print(f"    [3/3] Building DataFrames...")
            batch_mapping_df = create_sra_mapping_dataframe(pmid_to_sra_ids, sra_to_srp_ids)
            batch_srp_df = create_srp_centric_dataframe(pmid_to_sra_ids, sra_to_srp_ids)
            
            # Save batch files
            batch_mapping_df.to_csv(mapping_path_batch, index=False)
            batch_srp_df.to_csv(srp_path_batch, index=False)
            print(f"          Batch saved")

        print(f"\nAggregating batch results...")

        # STEP 4: Aggregate all batch files

        print(f"  Loading and combining {num_batches} batch files...")
        
        # Load all mapping batches
        all_mapping_dfs = []
        all_srp_dfs = []
        
        for batch_num in range(num_batches):
            batch_number = batch_num + 1
            mapping_path_batch = os.path.join(batch_dir, f'batch_{batch_number:04d}_mapping.csv')
            srp_path_batch = os.path.join(batch_dir, f'batch_{batch_number:04d}_srp.csv')
            
            if os.path.exists(mapping_path_batch) and os.path.exists(srp_path_batch):
                all_mapping_dfs.append(pd.read_csv(mapping_path_batch))
                all_srp_dfs.append(pd.read_csv(srp_path_batch))
        
        # Aggregate
        final_mapping_df = aggregate_batch_mappings(all_mapping_dfs)
        final_srp_df = aggregate_batch_srp_data(all_srp_dfs)
        
        # Save aggregated mapping (for reference)
        final_mapping_df.to_csv(sra_mapping_path, index=False)
        print(f"\n  PMID-SRA-SRP mapping saved to:")
        print(f"    {sra_mapping_path}")
        print(f"    Total PMIDs: {len(final_mapping_df)}")
        
        # Save aggregated SRP-centric data
        final_srp_df.to_csv(srp_centric_path, index=False)
        print(f"\n  SRP-centric data saved to:")
        print(f"    {srp_centric_path}")
        print(f"    Total unique SRP/ERP IDs: {len(final_srp_df)}")
    
    # STEP 5: Gather metadata for each unique SRP/ERP ID

    print(f"\nStep 5: Gathering metadata for SRP/ERP studies")
    
    unique_srp_ids = final_srp_df['SRP_ERP_ID'].tolist()
    print(f"  Gathering metadata for {len(unique_srp_ids)} unique SRP/ERP IDs...")
    
    # Create SRP to sample SRA mapping for efficient metadata fetching
    # (Need to reload the full mapping data for this step)
    if skip_to_metadata:
        # If we skipped batch processing, we need to rebuild the mappings
        print(f"  Reconstructing SRA mappings from batch files...")
        all_mapping_dfs = []
        batch_num = 1
        while True:
            mapping_path_batch = os.path.join(
                batch_dir, 
                f'batch_{batch_num:04d}_mapping.csv'
            )
            if os.path.exists(mapping_path_batch):
                all_mapping_dfs.append(pd.read_csv(mapping_path_batch))
                batch_num += 1
            else:
                break
        
        if all_mapping_dfs:
            full_mapping_df = pd.concat(all_mapping_dfs, ignore_index=True)
        else:
            # Fallback: load from aggregated file and parse
            full_mapping_df = pd.read_csv(sra_mapping_path)
        
        # Reconstruct the dictionaries from the mapping DataFrame
        pmid_to_sra_ids = {}
        sra_to_srp_ids = {}
        
        for _, row in full_mapping_df.iterrows():
            pmid = str(row['PMID'])
            sra_str = str(row['SRA_IDs']) if pd.notna(row['SRA_IDs']) else ''
            srp_str = str(row['SRP_ERP_IDs']) if pd.notna(row['SRP_ERP_IDs']) else ''
            
            # Parse semicolon-separated SRA IDs
            if sra_str:
                sra_list = [s.strip() for s in sra_str.split(';') if s.strip()]
                pmid_to_sra_ids[pmid] = sra_list
                
                # Parse semicolon-separated SRP IDs
                if srp_str:
                    srp_list = [s.strip() for s in srp_str.split(';') if s.strip()]
                    for sra_id in sra_list:
                        sra_to_srp_ids[sra_id] = srp_list
    
    # Create SRP to sample SRA mapping
    srp_to_sra_mapping = create_srp_to_sra_mapping(
        pmid_to_sra_ids, 
        sra_to_srp_ids
    )
    
    # Gather metadata using the mapping
    metadata_df = gather_srp_metadata(unique_srp_ids, srp_to_sra_mapping)
    
    # STEP 6: Enrich SRP data with metadata

    print(f"\nStep 6: Enriching SRP data with metadata")
    
    sra_datasets_df = enrich_srp_data_with_metadata(final_srp_df, metadata_df)
    
    # Link datasets with DOCs via publications, projects, and programs
    sra_datasets_df = get_dataset_doc_from_project(sra_datasets_df)
    
    # Generate UUID5 based on dataset_source_id
    print(f"\nGenerating UUID5 for each dataset...")
    uuid_fields = ['dataset_source_id']
    sra_datasets_df = get_composite_uuid5(
        sra_datasets_df, 
        uuid_fields, 
        uuid_col='dataset_uuid'
    )
    
    # Replace NaN with empty strings for consistency
    sra_datasets_df = sra_datasets_df.fillna('')
    
    # Save final enriched dataset
    sra_datasets_df.to_csv(sra_datasets_path, index=False)
    print(f"\n  Enriched SRA datasets saved to:")
    print(f"    {sra_datasets_path}")
    print(f"    Total unique studies: {len(sra_datasets_df)}")
    
    # Final summary
    print(f"\n---\nSRA Gathering Summary")
    
    total_pmids = len(final_mapping_df)
    pmids_with_sra = len(final_mapping_df[final_mapping_df['SRA_IDs'] != ''])
    pmids_with_srp = len(final_mapping_df[final_mapping_df['SRP_ERP_IDs'] != ''])
    
    print(f"  Total PMIDs processed:        {total_pmids:>8}")
    print(f"  PMIDs with SRA IDs:           {pmids_with_sra:>8} ({pmids_with_sra/total_pmids*100:>5.1f}%)")
    print(f"  PMIDs with SRP/ERP IDs:       {pmids_with_srp:>8} ({pmids_with_srp/total_pmids*100:>5.1f}%)")
    print(f"  Total unique SRP/ERP studies: {len(sra_datasets_df):>8}")
    print(f"---")
    
    return sra_datasets_df



# Run module as a standalone script when called directly
if __name__ == "__main__":
    
    print(f"Running {os.path.basename(__file__)} as standalone module...")

    # Load publications
    pubs_path = config.PUBLICATIONS_INTERMED_PATH
    if os.path.exists(pubs_path):
        publications = pd.read_csv(pubs_path)
        print(f"Publications data loaded from {pubs_path}")

        # Run the gathering workflow
        sra_datasets_df = gather_sra_data(publications, overwrite_intermeds=False)

    else:
        print(f"Error: Publications file not found at {pubs_path}")
        print("Please run the publication gathering workflow first.")

    # # TESTING
    # print(f"\n---\nTEST MODE")

    # # Small test params
    # publications_df = pd.DataFrame({'pmid': ['38738472', # Standard
    #                                          '10637239', # No SRA match
    #                                          '26829319', # ERP match
    #                                          '38260414', # Many-to-one (1/2)
    #                                          '38802751', # Many-to-one (2/2)
    #                                          ]}) 
    # batch_size = 2
    # overwrite_intermeds=False

    # print(f"TEST PARAM: PMID list: {publications_df['pmid'].tolist()}")
    # print(f"TEST PARAM: Batch size: {batch_size}")

    # # Run the test workflow
    # sra_datasets_df = gather_sra_data(publications_df, 
    #                                   overwrite_intermeds=overwrite_intermeds, 
    #                                   batch_size=batch_size)
    
    # print(f"\n---\nTEST MODE COMPLETE: ")
    # print(f"SRA intermediate CSV shape: {sra_datasets_df.shape}")
