"""
gather_sra_data.py
2026-01-23 ZD

This script defines functions to gather metadata for Sequence Read Archive (SRA)
studies supported by the NCI. The main function accepts a DataFrame containing 
publication information with PMIDs, finds all associated SRA IDs, maps them to 
study-level SRP/ERP IDs, and outputs structured CSV files for further analysis.

The workflow follows these steps:
1. Read input DataFrame with PubMed IDs
2. Map PMIDs to SRA IDs using NCBI's E-utilities
3. Map SRA IDs to SRP/ERP IDs (study-level identifiers)
4. Create two outputs:
   - PMID → SRA → SRP/ERP mapping file
   - SRP/ERP-centric file (unique studies with associated PMIDs)
5. Process in batches to handle large datasets

Intermediate files are saved at each step and reused if they already exist.
"""

import os
import sys
import json
import re
import time
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



def fetch_sra_ids(pmid: str) -> Tuple[str, List[str]]:
    """
    Fetch SRA IDs for a single PubMed ID
    
    Args:
        pmid: PubMed ID to query
    
    Returns:
        Tuple of (pmid, list_of_sra_ids)
    """

    try:
        # Small delay to control API rate
        time.sleep(0.1)
        
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
        print(f"Error processing PMID {pmid}: {e}")
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
        print(f"Error fetching SRP/ERP for SRA ID {sra_id}: {e}")
        return []



def get_srp_ids_for_sra_ids(
    sra_ids: List[str],
    max_workers: int = 10,
    sleep_between_calls: float = 0.05,
) -> Dict[str, List[str]]:
    """
    Batch map SRA IDs to SRP/ERP accessions using EFetch(db='sra').
    
    Args:
        sra_ids: List of SRA IDs to query
        max_workers: Number of concurrent threads (default 10)
        sleep_between_calls: Delay between API calls in seconds
    
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
    
    return pd.DataFrame(mapping_rows)



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
    
    return pd.DataFrame(srp_rows)



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
            pmids = row['PMIDs'].split('; ') if row['PMIDs'] else []
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



def gather_sra_data(
    publications_df: pd.DataFrame,
    overwrite_intermeds: bool = False,
    batch_size: int = 2000
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Main function for the SRA data gathering workflow. Uses NCBI E-Utilities to
    pull information for all SRA datasets associated with provided publications.
    
    Input and output filepaths are defined in config.py

    Args:
        publications_df: DataFrame output of publication gathering workflow
        overwrite_intermeds: Rerun all steps, even if intermediate files 
            already exist (default False). False will save time and avoid
            unnecessary repeated API calls.
        batch_size: Number of PMIDs to process per batch (default 2000)
        
    Returns:
        Tuple of (mapping_df, srp_df):
        - mapping_df: PMID → SRA IDs → SRP/ERP IDs mapping
        - srp_df: SRP/ERP-centric with associated PMIDs
    """

    print(f"\n---\nDATASETS - SRA:\n"
          f"Gathering, formatting, and saving SRA data...\n---\n")

    # Define paths from config
    sra_mapping_path = config.SRA_PMID_MAPPING_PATH
    srp_centric_path = config.SRA_SRP_CENTRIC_PATH
    batch_dir = config.SRA_BATCH_DIR

    # Create directories if they don't exist
    for path in [sra_mapping_path, srp_centric_path]:
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
    
    print(f"Loaded {len(publications_df)} publications.")
    print(f"Finding SRA IDs associated with {len(unique_pmids)} unique PMIDs...")

    # Check if final outputs exist and we're not overwriting
    if (os.path.exists(sra_mapping_path) and os.path.exists(srp_centric_path) 
        and not overwrite_intermeds):
        print(f"\n---\nFinal SRA outputs already exist. Loading from files...")
        mapping_df = pd.read_csv(sra_mapping_path)
        srp_df = pd.read_csv(srp_centric_path)
        
        print(f"Loaded mapping with {len(mapping_df)} rows")
        print(f"Loaded SRP-centric data with {len(srp_df)} unique studies")
        return mapping_df, srp_df

    # Process in batches
    total_pmids = len(unique_pmids)
    num_batches = (total_pmids + batch_size - 1) // batch_size

    print(f"Processing {total_pmids} PMIDs in {num_batches} batches...")
    print(f"Batch size: {batch_size}")

    all_mapping_dfs = []
    all_srp_dfs = []

    for batch_num in range(num_batches):
        start_idx = batch_num * batch_size
        end_idx = min((batch_num + 1) * batch_size, total_pmids)
        batch_pmids = unique_pmids[start_idx:end_idx]
        
        print(f"Processing batch {batch_num + 1}/{num_batches}")
        print(f"PMIDs {start_idx} to {end_idx - 1} ({len(batch_pmids)} PMIDs)")
        
        # Define output paths for this batch
        mapping_path_batch = os.path.join(
            batch_dir, 
            f'batch_{batch_num+1:04d}_mapping_{start_idx:06d}-{end_idx-1:06d}.csv'
        )
        srp_path_batch = os.path.join(
            batch_dir,
            f'batch_{batch_num+1:04d}_srp_{start_idx:06d}-{end_idx-1:06d}.csv'
        )
        
        # Check if batch files exist and we're not overwriting
        if (os.path.exists(mapping_path_batch) and os.path.exists(srp_path_batch) 
            and not overwrite_intermeds):
            print(f"  Loading existing batch files...")
            batch_mapping_df = pd.read_csv(mapping_path_batch)
            batch_srp_df = pd.read_csv(srp_path_batch)
        else:
            # Step 1: PMID → SRA IDs
            print(f"\n  [Step 1/3] Fetching SRA IDs for PMIDs...")
            pmid_to_sra_ids = get_sra_ids_for_pubmed_ids(batch_pmids)
            
            pmids_with_sra = len([p for p, sra_list in pmid_to_sra_ids.items() if sra_list])
            print(f"    Found SRA IDs for {pmids_with_sra}/{len(batch_pmids)} PMIDs")
            
            # Step 2: Get unique SRA IDs and fetch their SRP/ERP IDs
            all_sra_ids = [item for sublist in pmid_to_sra_ids.values() for item in sublist]
            unique_sra_ids = list(set(all_sra_ids))
            
            print(f"\n  [Step 2/3] Found {len(unique_sra_ids)} unique SRA IDs. Fetching SRP/ERP IDs...")
            
            if unique_sra_ids:
                sra_to_srp_ids = get_srp_ids_for_sra_ids(
                    unique_sra_ids,
                    max_workers=10,
                    sleep_between_calls=0.05
                )
            else:
                sra_to_srp_ids = {}
            
            # Check for SRA IDs that did not return any SRP/ERP IDs
            sra_without_srp = [sra_id for sra_id, srp_list in sra_to_srp_ids.items() if not srp_list]
            if sra_without_srp:
                print(f"\n    Warning: {len(sra_without_srp)} SRA ID(s) did not return any SRP/ERP IDs")
                print(f"       (showing first 5): {sra_without_srp[:5]}")
            
            sra_with_srp = len([sra for sra, srp_list in sra_to_srp_ids.items() if srp_list])
            print(f"    Found SRP/ERP IDs for {sra_with_srp}/{len(unique_sra_ids)} SRA IDs")
            
            # Step 3: Build DataFrames
            print(f"\n  [Step 3/3] Building DataFrames...")
            batch_mapping_df = create_sra_mapping_dataframe(pmid_to_sra_ids, sra_to_srp_ids)
            batch_srp_df = create_srp_centric_dataframe(pmid_to_sra_ids, sra_to_srp_ids)
            
            # Save batch files
            batch_mapping_df.to_csv(mapping_path_batch, index=False)
            batch_srp_df.to_csv(srp_path_batch, index=False)
            print(f"    Batch files saved")
        
        # Store for aggregation
        all_mapping_dfs.append(batch_mapping_df)
        all_srp_dfs.append(batch_srp_df)

    print("\nAll batches processed. Aggregating results...")

    # Aggregate all batches
    print("\nAggregating mapping data from all batches...")
    final_mapping_df = aggregate_batch_mappings(all_mapping_dfs)
    
    print("Aggregating SRP-centric data from all batches...")
    final_srp_df = aggregate_batch_srp_data(all_srp_dfs)
    
    # Save final outputs
    final_mapping_df.to_csv(sra_mapping_path, index=False)
    print(f"Complete mapping saved to: {sra_mapping_path}")
    print(f"Total rows: {len(final_mapping_df)}")
    
    final_srp_df.to_csv(srp_centric_path, index=False)
    print(f"Complete SRP-centric data saved to: {srp_centric_path}")
    print(f"Total unique SRP/ERP IDs: {len(final_srp_df)}")
    
    # Final summary
    total_pmids = len(final_mapping_df)
    pmids_with_sra = len(final_mapping_df[final_mapping_df['SRA_IDs'] != ''])
    pmids_with_srp = len(final_mapping_df[final_mapping_df['SRP_ERP_IDs'] != ''])
    
    print("\nFINAL SUMMARY")
    print(f"Total PMIDs processed: {total_pmids}")
    print(f"PMIDs with SRA IDs: {pmids_with_sra} ({pmids_with_sra/total_pmids*100:.1f}%)")
    print(f"PMIDs with SRP/ERP IDs: {pmids_with_srp} ({pmids_with_srp/total_pmids*100:.1f}%)")
    print(f"Total unique SRP/ERP IDs: {len(final_srp_df)}")
    print()
    
    return final_mapping_df, final_srp_df



# Run module as a standalone script when called directly
if __name__ == "__main__":
    
    print(f"Running {os.path.basename(__file__)} as standalone module...")

    # # Load publications
    # pubs_path = config.PUBLICATIONS_INTERMED_PATH
    # if os.path.exists(pubs_path):
    #     publications = pd.read_csv(pubs_path)
    #     print(f"Publications data loaded from {pubs_path}")

    #     # Run the gathering workflow
    #     mapping_df, srp_df = gather_sra_data(publications, overwrite_intermeds=False)

    # else:
    #     print(f"Error: Publications file not found at {pubs_path}")
    #     print("Please run the publication gathering workflow first.")

    # TESTING
    print(f"\n---\nTEST MODE\n---\n")

    # Small test params
    publications_df = pd.DataFrame({'pmid': ['38738472','38227896','10637239']})
    batch_size = 2
    overwrite_intermeds=False

    print(f"TEST: PMID list: {publications_df['pmid'].tolist()}")
    print(f"TEST: Batch size: {batch_size}")

    # Run the test workflow
    mapping_df, srp_df = gather_sra_data(publications_df, 
                                         overwrite_intermeds=True, 
                                         batch_size=batch_size)
