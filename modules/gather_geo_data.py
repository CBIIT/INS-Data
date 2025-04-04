"""
gather_geo_data.py
2025-03-06 ZD

This script defines functions to gather metadata for Gene Expression Omnibus (GEO)
studies supported by the NCI. The main function accepts a CSV containing publication
information with PMIDs, finds all associated GEO IDs and their metadata, and outputs
a structured CSV file for further analysis.

The workflow follows these steps:
1. Read input CSV with 'coreproject' and 'pmid' columns
2. Map PMIDs to GEO IDs using NCBI's E-utilities
3. Gather Esummary metadata for each GEO ID
4. Gather FTP metadata for each GEO ID
5. Process and combine the metadata
6. Output results as a TSV file

Intermediate files are saved at each step and reused if they already exist.
"""

import io
import gzip
import os
import sys
import json
import re
import time
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Any

from ftplib import FTP
import pandas as pd
import requests
import concurrent.futures
from urllib.parse import urlparse
from tqdm import tqdm  # for progress bars
from Bio import Entrez  # for e-Utils API

# Append the project's root directory to the Python path
# This allows for importing config when running as part of main.py or alone
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import config



def fetch_geo_ids(pmid: str) -> Tuple[str, List[str]]:
    """
    Fetch GEO IDs for a single PubMed ID
    
    Args:
        pmid: PubMed ID to query
    
    Returns:
        Tuple of (pmid, list_of_geo_ids)
    """

    try:
        # Small delay to control API rate
        time.sleep(0.1)
        
        link_handle = Entrez.elink(
            dbfrom="pubmed",
            db="gds",
            id=pmid,
            linkname="pubmed_gds"
        )
        
        link_record = Entrez.read(link_handle)
        link_handle.close()
        
        geo_ids = [
            link['Id']
            for link_set in link_record
            for link in link_set.get('LinkSetDb', [])
            for link in link.get('Link', [])
        ]
        
        return (pmid, geo_ids)
    
    except Exception as e:
        print(f"Error processing PMID {pmid}: {e}")
        return (pmid, [])



def get_geo_ids_for_pubmed_ids(pubmed_ids: List[str]) -> Dict[str, List[str]]:
    """
    Retrieve GEO dataset IDs associated with each PubMed ID in a list of PMIDs. 
    
    Args:
        pubmed_ids: List of PubMed IDs to query
    
    Returns:
        Dictionary mapping PubMed IDs to their associated GEO IDs
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
    pmid_geo_links = {}
    with concurrent.futures.ThreadPoolExecutor(max_workers=max_workers) as executor:
        # Submit tasks (PMIDs) to the executor and store futures in dict
        futures = {executor.submit(fetch_geo_ids, pmid): pmid for pmid in pubmed_ids}
        
        # Iterate through futures (PMIDs) as they become available
        for future in tqdm(concurrent.futures.as_completed(futures),
                          unit="PMID", total=pmid_count, ncols=80,
                          desc="Fetching GEO IDs"):
            pmid, geo_ids = future.result()
            pmid_geo_links[pmid] = geo_ids
    
    return pmid_geo_links



def get_full_geo_record(geo_id: str) -> Optional[Dict]:
    """
    Get the full NCBI ESummary record for a given GEO Accession.
    
    Args:
        geo_id: GEO ID from Accession (e.g. '200252411'). Note that GEO 
            Accession sometimes appears with GSE prefix which is replaced by 
            '200' in the id. 
    
    Returns:
        Full record of data available via Entrez or None if error
    """

    try:
        handle = Entrez.esummary(db="gds", id=geo_id, retmode="text")
        records = Entrez.read(handle)
        handle.close()
        return records
    except Exception as e:
        print(f"Error fetching information for GEO Accession {geo_id}: {e}")
        return None



def get_all_geo_summary_records(geo_id_list: List[str]) -> List[Dict]:
    """
    Get NCBI ESummary metadata for each GEO accession in a list of GEO IDs. 
    
    Args:
        geo_id_list: List of GEO IDs
    
    Returns:
        List of metadata records
    """

    # Configure Entrez
    # Email and api key from hidden local env file. Use default if not defined
    Entrez.email = os.environ.get('NCBI_EMAIL', 'your-email@example.com')
    Entrez.api_key = os.environ.get('NCBI_API_KEY', '')
    if not Entrez.api_key: 
        print(f"WARNING: No NCBI API key in use. Check readme and local .env file."
              f"\nNCBI E-Utilities rate will be limited and may cause errors.")

    # Build empty list to hold gathered records
    all_records = []
    
    # Iterate through list of ids to get metadata for each
    for geo_id in tqdm(geo_id_list, ncols=80, desc="Fetching GEO metadata"):
        record = get_full_geo_record(geo_id)
        
        if record:
            all_records.append(record)
        else:
            print(f"Error retrieving {geo_id} metadata. Skipping.")
    
    return all_records



def extract_ftp_metadata(ftp_link: str) -> Dict[str, Any]:
    """
    Get full GEO metadata from series matrix files given an FTP link.
    
    Args:
        ftp_link: Full FTP link (can start with ftp:// or be a direct path)
        
    Returns:
        Dictionary of extracted metadata fields
    """

    # Clean up the FTP link
    if ftp_link.startswith('ftp://'):
        parsed = urlparse(ftp_link)
        ftp_dir = parsed.path
    else:
        ftp_dir = ftp_link
    
    # Remove any trailing slashes
    ftp_dir = ftp_dir.rstrip('/')
    
    # Check that location is the matrix subdirectory
    if not ftp_dir.endswith('/matrix'):
        ftp_dir = f"{ftp_dir}/matrix"
    
    # Connect to FTP
    ftp = FTP('ftp.ncbi.nlm.nih.gov')
    ftp.login()
    
    try:
        # Remove leading slash if present
        ftp_dir = ftp_dir.lstrip('/')
        ftp.cwd(ftp_dir)
        
        matrix_files = [f for f in ftp.nlst() if f.endswith('series_matrix.txt.gz')]
        
        if not matrix_files:
            raise FileNotFoundError(f"No series matrix files found in {ftp_dir}")
        
        metadata_dict = {}
        
        # Process each matrix file
        for matrix_file in matrix_files:
            # Read and decompress file directly from FTP
            bio = io.BytesIO()
            ftp.retrbinary(f"RETR {matrix_file}", bio.write)
            bio.seek(0)
            
            with gzip.GzipFile(fileobj=bio, mode='rb') as gz:
                content = gz.read().decode('utf-8', errors='replace')
                
                # Extract all metadata fields
                field_pattern = r'!(Series|Sample)_(\w+)\t(.+)(?:\r\n|\r|\n)'
                metadata_fields = re.findall(field_pattern, content)
                
                # Organize fields into categories
                for prefix, field, value in metadata_fields:
                    key = f"{prefix}_{field}"
                    # Clean up the value (remove quotes)
                    value = value.strip('"')
                    
                    # Accumulate values for repeating fields
                    if key in metadata_dict:
                        if isinstance(metadata_dict[key], list):
                            metadata_dict[key].append(value)
                        else:
                            metadata_dict[key] = [metadata_dict[key], value]
                    else:
                        metadata_dict[key] = value
        
        return metadata_dict
        
    except Exception as e:
        print(f"Error extracting FTP metadata from {ftp_link}: {e}")
        return {}
    
    finally:
        ftp.quit()



def create_geo_dataframe(pmid_geo_links: Dict[str, List[str]], 
                         publication_df: pd.DataFrame) -> pd.DataFrame:
    """
    Convert PubMed to GEO ID mapping to a structured DataFrame with project info.
    
    Args:
        pmid_geo_links: Mapping of PubMed IDs to GEO IDs
        publication_df: Original publication DataFrame with project information
    
    Returns:
        DataFrame with columns 'pmid', 'geo_id', and 'coreproject'
    """
    # Create pmid to coreproject mapping for quick lookup
    pmid_to_core = dict(zip(publication_df['pmid'], publication_df['coreproject']))
    
    df_data = []
    for pmid, geo_ids in pmid_geo_links.items():
        coreproject = pmid_to_core.get(pmid, '')
        
        if not geo_ids:
            df_data.append({'pmid': pmid, 'geo_id': pd.NA, 'coreproject': coreproject})
        else:
            for geo_id in geo_ids:
                df_data.append({'pmid': pmid, 'geo_id': geo_id, 'coreproject': coreproject})
    
    return pd.DataFrame(df_data)



def get_all_geo_ftp_metadata(geo_records: List[Dict]) -> Dict[str, Dict]:
    """
    Extract FTP metadata for multiple GEO records.
    
    Args:
        geo_records: List of GEO metadata records
    
    Returns:
        Dictionary mapping GEO IDs to their FTP metadata
    """
    ftp_metadata = {}
    
    for record in tqdm(geo_records, desc="Extracting FTP metadata", ncols=80):
        try:
            # Extract GEO ID and FTP link from record
            record = record[0]
            geo_id = record.get('Id', 'unknown')
            ftp_link = record.get('FTPLink', '')
            
            if ftp_link:
                metadata = extract_ftp_metadata(ftp_link)
                ftp_metadata[geo_id] = metadata
            else:
                tqdm.write(f"No FTP link found for GEO ID {geo_id}")
        except Exception as e:
            tqdm.write(f"Error processing FTP metadata for record: {e}")
    
    return ftp_metadata



# def file_exists_with_content(filepath: str) -> bool:
#     """
#     Check if a file exists and has content.
    
#     Args:
#         filepath: Path to file
        
#     Returns:
#         True if file exists and has content, False otherwise
#     """
#     if not os.path.exists(filepath):
#         return False
    
#     try:
#         if filepath.endswith('.csv') or filepath.endswith('.tsv'):
#             # Check if dataframe has content
#             df = pd.read_csv(filepath, sep='\t' if filepath.endswith('.tsv') else ',')
#             return not df.empty
#         elif filepath.endswith('.json'):
#             # Check if json has content
#             with open(filepath, 'r') as f:
#                 content = json.load(f)
#             return bool(content)
#         else:
#             # For other file types, check if size is > 0
#             return os.path.getsize(filepath) > 0
#     except Exception as e:
#         print(f"Error checking file {filepath}: {e}")
#         return False


def gather_geo_data(input_csv: str, 
                    output_csv: str, 
                    esummary_intermed_path: str, 
                    raw_ftp_path: str,
                    geo_mapping_path: str,
                    overwrite_intermeds: bool=False) -> pd.DataFrame:
    """
    Main function to orchestrate the entire GEO data gathering workflow.
    Implements checkpoint loading to skip completed steps.
    
    Args:
        input_csv: Path to input CSV with 'coreproject' and 'pmid' columns
        output_csv: Path for output CSV with GEO metadata
        raw_metadata_path: Path to save consolidated raw metadata JSON
        raw_ftp_path: Path to save consolidated FTP metadata JSON
        geo_mapping_path: Path to save GEO-PMID-Project mapping CSV
        overwrite_intermeds: Rerun all steps, even if intermediate files 
            already exist (default False). True will save time and avoid
            unnecessary API calls. 
        
    Returns:
        DataFrame with processed GEO metadata
    """
    print(f"Starting GEO data gathering workflow...")
    
    # Create directories if they don't exist
    for path in [esummary_intermed_path, raw_ftp_path, geo_mapping_path, output_csv]:
        os.makedirs(os.path.dirname(path), exist_ok=True)
    
    # Read input publication data
    print(f"Reading publication data from {input_csv}")
    pub_df = pd.read_csv(input_csv)
    
    # Validate required columns
    required_cols = ['coreproject', 'pmid']
    missing_cols = [col for col in required_cols if col not in pub_df.columns]
    if missing_cols:
        raise ValueError(f"Input CSV missing required columns: {', '.join(missing_cols)}")
    
    # Convert PMIDs to strings for consistency
    pub_df['pmid'] = pub_df['pmid'].astype(str)
    print(f"Read {len(pub_df)} publications")


    # Load or create GEO-PMID-Project mapping intermediate
    if os.path.exists(geo_mapping_path) and not overwrite_intermeds:
        print(f"Loading existing GEO-PMID-Project mapping from {geo_mapping_path}")
        geo_mapping_df = pd.read_csv(geo_mapping_path, 
                                     dtype={'geo_id': str, 'pmid':str})
    else:
        print(f"Generating GEO mapping...")
        # Extract unique PMIDs
        unique_pmids = pub_df['pmid'].unique().tolist()
        print(f"Found {len(unique_pmids)} unique PMIDs")
        
        # Get GEO IDs for each PMID
        pmid_geo_mapping = get_geo_ids_for_pubmed_ids(unique_pmids)
        
        # Create mapping DataFrame
        geo_mapping_df = create_geo_dataframe(pmid_geo_mapping, pub_df)
        
        # Save ID mapping as checkpoint file
        print(f"Saving GEO-PMID-Project mapping to {geo_mapping_path}")
        geo_mapping_df.to_csv(geo_mapping_path, index=False)
    
    # Filter to rows with valid GEO IDs
    valid_geo_df = geo_mapping_df.dropna(subset=['geo_id']).copy()
    print(f"Found {len(valid_geo_df)} GEO entries linked to "
          f"{len(valid_geo_df['pmid'].unique())} PMIDs")
    
    # Get unique GEO IDs for metadata gathering
    unique_geo_ids = valid_geo_df['geo_id'].unique().tolist()
    

    # Load or create GEO Esummary API response intermediate
    if os.path.exists(esummary_intermed_path) and not overwrite_intermeds:
        print(f"Loading existing GEO ESummary metadata from {esummary_intermed_path}")
        with open(esummary_intermed_path, 'r') as f:
            raw_records = json.load(f)
    else:
        print(f"Retrieving metadata for {len(unique_geo_ids)} unique GEO datasets")
        raw_records = get_all_geo_summary_records(unique_geo_ids)
        
        # Save summary metadata checkpoint
        print(f"Saving NCBI ESummary GEO metadata to {esummary_intermed_path}")
        with open(esummary_intermed_path, 'w') as f:
            json.dump(raw_records, f, indent=2)


    # Load or create GEO FTP metadata intermediate
    if os.path.exists(raw_ftp_path) and not overwrite_intermeds:
        print(f"Loading existing FTP metadata from {raw_ftp_path}")
        with open(raw_ftp_path, 'r') as f:
            ftp_metadata = json.load(f)
    else:
        print("Extracting FTP metadata from GEO records")
        ftp_metadata = get_all_geo_ftp_metadata(raw_records)
        
        # Save consolidated FTP metadata
        print(f"Saving all FTP metadata to {raw_ftp_path}")
        with open(raw_ftp_path, 'w') as f:
            json.dump(ftp_metadata, f, indent=2)


    # Process records into combined dataframe
    print("Processing metadata into final dataset")
    processed_data = []
    
    for record in tqdm(raw_records, desc="Processing records", ncols=80):
        try:
            # Extract GEO ID from the first item in the record list
            record_data = record[0]  # This is correct - each record is a list with one dict
            geo_id = record_data.get('Id', '')

            if not geo_id:
                print(f"WARNING: No geo_id found in record. Skipping.")
                continue
            
            # Debug print to check for matches
            matching_rows = valid_geo_df[valid_geo_df['geo_id'] == geo_id]
            if matching_rows.empty:
                print(f"No matches found for GEO ID: {geo_id}")
                # Try a string comparison in case of type mismatch
                matching_rows = valid_geo_df[valid_geo_df['geo_id'].astype(str) == str(geo_id)]
                if not matching_rows.empty:
                    print(f"  Found match after string conversion")
            
            for _, row in matching_rows.iterrows():
                entry = {
                    'geo_id': geo_id,
                    'pmid': row['pmid'],
                    'coreproject': row['coreproject'],
                    'title': record_data.get('title', ''),
                    'summary': record_data.get('summary', ''),
                    'geo_accession': record_data.get('Accession', ''),
                    'n_samples': record_data.get('n_samples', ''),
                    'related_terms': record_data.get('taxon', ''),
                    'release_date': record_data.get('PDAT', ''),
                    'assay_method': record_data.get('gdsType', ''),
                    'study_links': record_data.get('FTPLink', ''),
                }
                
                processed_data.append(entry)
        except Exception as e:
            print(f"Error processing record: {e}")
    
    # Create final dataframe
    result_df = pd.DataFrame(processed_data)
    
    # Save results
    print(f"Saving processed results to {output_csv}")
    result_df.to_csv(output_csv, index=False, sep='\t')
    
    print(f"GEO data gathering complete. Found {len(result_df)} records.")
    return result_df


# Run module as a standalone script when called directly
if __name__ == "__main__":
    print(f"Running {os.path.basename(__file__)} as standalone module...")
    
    # Get the directory of the current script
    script_dir = os.path.dirname(os.path.abspath(__file__))

    # Construct absolute paths
    input_path = os.path.join(script_dir, "..", "data", "01_intermediate", "geo_test", "publication.csv")
    output_path = os.path.join(script_dir, "..", "data", "02_output", "geo_test", "geo_datasets.tsv")
    raw_metadata_path = os.path.join(script_dir, "..", "data", "01_intermediate", "geo_test", "geo_metadata.json")
    raw_ftp_path = os.path.join(script_dir, "..", "data", "01_intermediate", "geo_test", "geo_ftp_metadata.json")
    geo_mapping_path = os.path.join(script_dir, "..", "data", "01_intermediate", "geo_test", "geo_pmid_project_map.csv")
    
    # Execute the main workflow
    result = gather_geo_data(input_path, output_path, raw_metadata_path, raw_ftp_path, geo_mapping_path)