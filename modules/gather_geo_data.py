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
3. Gather raw metadata for each GEO ID
4. Process and clean the metadata
5. Output results as a CSV file

Results are stored both as consolidated raw JSON files and as a final 
cleaned CSV for loading into INS.
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
        api_key: Whether an API key is being used
    
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
        api_key: Whether an API key is being used
    
    Returns:
        Dictionary mapping PubMed IDs to their associated GEO IDs
    """

    # Configure Entrez
    # Get user email from hidden local env file. Use default if not defined
    Entrez.email = os.environ.get('NCBI_EMAIL', 'your-email@example.com')
    Entrez.api_key = os.environ.get('NCBI_API_KEY', '')
    if Entrez.api_key == None: 
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
    Get single-GEO metadata from series matrix files given an FTP link.
    
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
        
        metadata = {
            'contributors': []}
        
        # Process each matrix file
        for matrix_file in matrix_files:
            # Read and decompress file directly from FTP
            bio = io.BytesIO()
            ftp.retrbinary(f"RETR {matrix_file}", bio.write)
            bio.seek(0)
            
            with gzip.GzipFile(fileobj=bio, mode='rb') as gz:
                content = gz.read().decode('utf-8')
                
                # Extract fields using regex
                contributors = re.findall(r'!Series_contributor\t(.+)(?:\r\n|\r|\n)', content)
                contributors = [contributor.strip('"') for contributor in contributors]
                metadata['contributors'].extend(contributors)
        
        # Remove duplicate contributors while preserving order
        metadata['contributors'] = list(dict.fromkeys(metadata['contributors']))
        return metadata
        
    except Exception as e:
        print(f"Error extracting FTP metadata from {ftp_link}: {e}")
        return {'contributors': []}
    
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



# def get_all_geo_ftp_metadata(geo_records: List[Dict]) -> Dict[str, Dict]:
#     """
#     Extract FTP metadata for multiple GEO records.
    
#     Args:
#         geo_records: List of GEO metadata records
    
#     Returns:
#         Dictionary mapping GEO IDs to their FTP metadata
#     """
#     ftp_metadata = {}
    
#     for record in tqdm(geo_records, desc="Extracting FTP metadata", ncols=80):
#         try:
#             # Extract GEO ID and FTP link from record
#             # This is a placeholder - adjust based on actual record structure
#             geo_id = record.get('Id', 'unknown')
#             ftp_link = record.get('FTPLink', '')
            
#             if ftp_link:
#                 metadata = extract_ftp_metadata(ftp_link)
#                 ftp_metadata[geo_id] = metadata
#             else:
#                 print(f"No FTP link found for GEO ID {geo_id}")
#         except Exception as e:
#             print(f"Error processing FTP metadata for record: {e}")
    
#     return ftp_metadata



def gather_geo_data(input_csv: str, output_csv: str, raw_metadata_path: str, raw_ftp_path: str) -> pd.DataFrame:
    """
    Main function to orchestrate the entire GEO data gathering workflow.
    
    Args:
        input_csv: Path to input CSV with 'coreproject' and 'pmid' columns
        output_csv: Path for output CSV with GEO metadata
        raw_metadata_path: Path to save consolidated raw metadata JSON
        raw_ftp_path: Path to save consolidated FTP metadata JSON
        
    Returns:
        DataFrame with processed GEO metadata
    """
    print(f"Starting GEO data gathering workflow...")
    
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
    
    # Extract unique PMIDs
    unique_pmids = pub_df['pmid'].unique().tolist()
    print(f"Found {len(unique_pmids)} unique PMIDs")
    
    # Get GEO IDs for each PMID
    pmid_geo_mapping = get_geo_ids_for_pubmed_ids(unique_pmids)
    
    # Create mapping DataFrame
    geo_mapping_df = create_geo_dataframe(pmid_geo_mapping, pub_df)
    
    # Filter to rows with valid GEO IDs
    valid_geo_df = geo_mapping_df.dropna(subset=['geo_id'])
    print(f"Found {len(valid_geo_df)} GEO entries linked to {len(valid_geo_df['pmid'].unique())} PMIDs")
    
    # Get unique GEO IDs for metadata gathering
    unique_geo_ids = valid_geo_df['geo_id'].unique().tolist()
    
    # Gather raw metadata
    print(f"Retrieving metadata for {len(unique_geo_ids)} unique GEO datasets")
    raw_records = get_all_geo_summary_records(unique_geo_ids)
    
    # Save summary metadata checkpoint
    print(f"Saving NCBI ESummary GEO metadata to {raw_metadata_path}")
    os.makedirs(os.path.dirname(raw_metadata_path), exist_ok=True)
    with open(raw_metadata_path, 'w') as f:
        json.dump(raw_records, f, indent=2)
    
    # # Extract and gather FTP metadata
    # print("Extracting FTP metadata from GEO records")
    # ftp_metadata = get_all_geo_ftp_metadata(raw_records)
    
    # # Save consolidated FTP metadata
    # print(f"Saving consolidated FTP metadata to {raw_ftp_path}")
    # os.makedirs(os.path.dirname(raw_ftp_path), exist_ok=True)
    # with open(raw_ftp_path, 'w') as f:
    #     json.dump(ftp_metadata, f, indent=2)
    
    # Process records into a final dataframe
    # This is placeholder code - you'll need to implement the actual processing
    print("Processing metadata into final dataset")
    processed_data = []
    
    for record in tqdm(raw_records, desc="Processing records", ncols=80):
        # Extract GEO ID
        record = record[0]
        geo_id = record.get('Id', '')
        if geo_id == None:
            print(f"WARNING: No geo_id found in {record}. Stopping loop.")
            break
        
        # Find matching rows in geo_mapping_df
        matching_rows = valid_geo_df[valid_geo_df['geo_id'] == geo_id]
        
        # # Get FTP metadata if available
        # ftp_data = ftp_metadata.get(geo_id, {})
        
        for _, row in matching_rows.iterrows():
            processed_data.append({
                'geo_id': geo_id,
                'pmid': row['pmid'],
                'coreproject': row['coreproject'],
                'title': record.get('title', ''),
            })
    
    # Create final dataframe
    result_df = pd.DataFrame(processed_data)
    
    # Save results
    print(f"Saving processed results to {output_csv}")
    os.makedirs(os.path.dirname(output_csv), exist_ok=True)
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

    # # Define paths
    # input_path = "data/01_intermediate/geo_test/publication.csv"
    # output_path = "data/02_output/geo_test/geo_datasets.tsv"
    # raw_metadata_path = "data/01_intermediate/geo_test/geo_metadata.json"
    # raw_ftp_path = "data/01_intermediate/geo_test/geo_ftp_metadata.json"
    
    # Execute the main workflow
    result = gather_geo_data(input_path, output_path, raw_metadata_path, raw_ftp_path)