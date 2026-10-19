"""
package_output_data.py
2023-12-21 ZD

This script defines primary function package_output_data that applies final 
formatting, filtering, etc. to the intermediate data gathering outputs in order
to prepare them for INS ingestion. These outputs are saved as TSVs. 
"""

import os
import unicodedata
import uuid
import html
from html.entities import html5 as HTML5_ENTITIES
import hashlib

import pandas as pd
import re

import config


# Dataset-like datatypes, all mapped to the general 'dataset' node type
DATASET_TYPES = {'dbgap_dataset',
                 'geo_dataset',
                 'sra_dataset',
                 'cedcd_dataset',
                 'ctd2_dataset',
                 'dceg_dataset',
                 'nccr_dataset',
                 'dctd_dataset'}

# Dataset fields allowed to retain HTML entities for rich-text display
DATASET_HTML_ENTITY_EXCLUDE_COLS = ['description', 'experimental_approaches']

# Fields that keep legacy ASCII-only character processing
ASCII_ONLY_COLS = ['description', 'experimental_approaches']

# Non-ASCII symbols kept in output alongside accented Latin and Greek letters
PRESERVED_SYMBOLS = set('°±×µ≤≥¹²³⁰⁴⁵⁶⁷⁸⁹₀₁₂₃₄₅₆₇₈₉')

# Unicode blocks whose letters (accented Latin and Greek) are kept in output
PRESERVED_LETTER_RANGES = [
    (0x00C0, 0x024F),   # Latin-1 Supplement, Latin Extended-A/B (í, ł)
    (0x0370, 0x03FF),   # Greek and Coptic (β, Δ)
    (0x1E00, 0x1EFF),   # Latin Extended Additional (ễ, Ḥ)
    (0x1F00, 0x1FFF),   # Greek Extended (ἄ)
]

# Excel ANSI saves replace unsupported characters with '?' (e.g. RORγt -> ROR?t)
# A '?' followed by key= starts a URL query string, so it is not flagged
SUBSTITUTED_CHAR_PATTERN = re.compile(r'[^\W\d_]\?(?![\w.-]*=)[^\W\d_]')

# UTF-8 bytes misread as cp1252 (e.g. é -> Ã©, ’ -> â€™)
CP1252_CONTINUATION_CHARS = bytes(range(0x80, 0xC0)).decode(
    'cp1252', errors='ignore')
MOJIBAKE_PATTERN = re.compile(
    '[\u00c2-\u00df][' + re.escape(CP1252_CONTINUATION_CHARS) + ']|\u00e2\u20ac')

# URLs legitimately contain '?'; stop at HTML delimiters so link text is scanned
URL_PATTERN = re.compile(r'(?:https?://|www\.)[^\s<>"\']+')

# Numeric (semicolon optional) and complete named references only, so URL
# parameters like &region= are never read as &reg
HTML_REFERENCE_PATTERN = re.compile(
    r'&(?:#[0-9]+;?|#[xX][0-9a-fA-F]+;?|[A-Za-z][A-Za-z0-9]*;)')

# Output QA patterns: HTML codes (semicolon optional for numeric), leftover
# tags, control characters, and any non-ASCII character
HTML_ENTITY_PATTERN = re.compile(r'&(?:#\d+;?|#[xX][0-9a-fA-F]+;?|[a-zA-Z]+;)')
HTML_TAG_PATTERN = re.compile(r'</?\w+>')
CONTROL_CHAR_PATTERN = re.compile(r'[\x00-\x1f\x7f]')
NON_ASCII_PATTERN = re.compile(r'[^\x00-\x7f]')

# '?' placeholder for a lost character (e.g. 5-10 ?l, (?DC1s)), not a query string
PLACEHOLDER_PATTERN = re.compile(r'\?(?![\w.-]*=)(?=[^\W_])')

# Plain-language QA questions and severity for each output check
OUTPUT_QA_CHECKS = {
    'file_encoding': ('FAIL',
        'Is any file not plain UTF-8, or does it start with a byte order mark?'),
    'replacement_char': ('FAIL',
        'Are there any U+FFFD replacement characters (diamond question mark)?'),
    'garbled_text': ('FAIL',
        'Is there garbled text from encoding mix-ups (UTF-8 read as cp1252)?'),
    'html_entity': ('FAIL',
        'Are there HTML codes (e.g. &#237;) outside dataset description '
        'and experimental_approaches?'),
    'html_tag': ('FAIL',
        'Are there leftover HTML tags in columns configured for tag removal?'),
    'unexpected_character': ('FAIL',
        'Are there characters outside the allowed set (ASCII, accented/Greek '
        'letters, approved symbols; ASCII only in description fields)?'),
    'control_character': ('FAIL',
        'Are there hidden control characters (tabs, line breaks) in values?'),
    'placeholder_char': ('INFO',
        "How many '?' placeholders mark characters lost in source data?"),
}


def warn_on_encoding_damage(df, source, max_examples=10):
    """Print warnings for text likely damaged by an Excel save and return
    a list of finding dicts (source_file, issue, excel_row, column, snippet)."""

    checks = {
        "'?' between letters (character likely replaced by Excel)":
            SUBSTITUTED_CHAR_PATTERN,
        "garbled characters (UTF-8 text read as cp1252)": MOJIBAKE_PATTERN,
        "U+FFFD replacement characters": re.compile('\ufffd'),
    }

    findings = []
    for issue, pattern in checks.items():
        issue_findings = []
        for col in df.columns:
            if not (pd.api.types.is_object_dtype(df[col])
                    or pd.api.types.is_string_dtype(df[col])):
                continue
            for idx, value in df[col].items():
                if not isinstance(value, str):
                    continue
                text = URL_PATTERN.sub(' ', value)
                match = pattern.search(text)
                if match:
                    snippet = text[max(0, match.start() - 20):match.end() + 20]
                    issue_findings.append({
                        'source_file': str(source),
                        'issue': issue,
                        # Header is Excel row 1, so first data row is row 2
                        'excel_row': idx + 2,
                        'column': col,
                        'snippet': snippet,
                    })

        if issue_findings:
            print(f"WARNING: {len(issue_findings)} value(s) in {source} "
                  f"contain {issue}. Fix in the curated file and re-run:")
            for finding in issue_findings[:max_examples]:
                print(f"    Excel row {finding['excel_row']}, "
                      f"{finding['column']}: ...{finding['snippet']}...")
        findings.extend(issue_findings)

    return findings



def write_special_character_scan(findings, report_path):
    """Save every special character finding to a CSV report, one per row.
    No report is written when there are no findings."""

    if not findings:
        print(f"Special character scan: 0 finding(s). No report written.")
        return

    columns = ['source_file', 'issue', 'excel_row', 'column', 'snippet']
    os.makedirs(os.path.dirname(report_path), exist_ok=True)
    # BOM so snippets with accents/symbols display correctly in Excel
    pd.DataFrame(findings, columns=columns).to_csv(
        report_path, index=False, encoding='utf-8-sig')
    print(f"Special character scan: {len(findings)} finding(s) saved to "
          f"{report_path}")



def read_csv_with_encoding_fallback(filepath, findings=None, **kwargs):
    """Read a curated CSV/TSV using the encoding Excel may have saved it in,
    then warn about likely Excel encoding damage. Findings are appended to
    the optional findings list for reporting.

    - UTF-16 (has BOM): Excel 'Unicode Text' save. Keeps all characters.
    - UTF-8, with or without BOM: pipeline outputs or Excel 'CSV UTF-8'.
    - cp1252 fallback: Excel 'Text (Tab delimited)' save. Accented letters
      survive, but Excel has already replaced characters outside cp1252
      (e.g. β, ≤) with '?'.
    """
    with open(filepath, 'rb') as f:
        bom = f.read(2)

    if bom in (b'\xff\xfe', b'\xfe\xff'):
        df = pd.read_csv(filepath, encoding='utf-16', **kwargs)
    else:
        try:
            df = pd.read_csv(filepath, encoding='utf-8-sig', **kwargs)
        except UnicodeDecodeError:
            print(f"WARNING: {filepath} is not valid utf-8 (likely saved from "
                  f"Excel). Retrying with cp1252 encoding.")
            df = pd.read_csv(filepath, encoding='cp1252', **kwargs)

    found = warn_on_encoding_damage(df, filepath)
    if findings is not None:
        findings.extend(found)

    return df



def add_type_column(df, datatype):
    """Add a column for datatype and fill with specified datatype."""

    df_with_type = df.copy()
    
    # Map to 'dataset' if in the list, otherwise keep the original value
    df_with_type['type'] = 'dataset' if datatype in DATASET_TYPES else datatype
    
    return df_with_type



def reorder_columns(df, column_configs, datatype):
    """Reorder and validate columns of df using dict from config.py. 
    Report any dropped columns and throw error if defined column not found.
    """
    # Check if the datatype exists in column_configs
    if datatype not in column_configs:
        raise ValueError(f"Invalid datatype: {datatype}")

    # Extract relevant configuration for the datatype
    config = column_configs[datatype]

    # Extract relevant information from the configuration
    keep_and_rename = config.get('keep_and_rename', {})

    # Define column names from old:new dict
    old_cols = set(keep_and_rename.keys())
    new_cols = set(keep_and_rename.values())

    # Get missing columns from current df
    missing_columns = old_cols - set(df.columns)

    # Check if old columns all exist or are already renamed
    if missing_columns and not new_cols.issubset(set(df.columns)):
        raise ValueError(f"Fields missing from intermediate {datatype} data: "
                         f"{', '.join(missing_columns)}")

    # Rename columns using keep_and_rename
    df.rename(columns=keep_and_rename, inplace=True)

    # Define list of columns to keep
    columns_to_keep = list(keep_and_rename.values())

    # Drop any columns not part of keep_and_rename
    dropped_columns = set(df.columns) - set(columns_to_keep)
    if dropped_columns:
        print(f"Columns dropped from {datatype} output: "
              f"{', '.join(dropped_columns)}")

    # Keep and reorder columns
    df_reordered = df[columns_to_keep]

    return df_reordered



def validate_first_columns(df, column_configs, datatype):
    """Validate that type, node_id, and link_id (opt) are first columns."""

    # Check if the datatype exists in column_configs and then get it
    if datatype not in column_configs:
        raise ValueError(f"Invalid datatype: {datatype}")
    config = column_configs.get(datatype)

    # Get values from column configuration
    node_id = config.get('node_id')
    link_id = config.get('link_id', None)

    # If link_id exists, add it to the list of expected columns 
    expected_cols = ['type', node_id]
    if link_id:
        expected_cols.append(link_id)

    # Check the order of columns in df
    if list(df.columns[:len(expected_cols)]) != expected_cols:
        raise ValueError(f"First columns in {datatype} output are not in the "
                         f"expected order. Reorder columns in config.py to fix.\n"
                         f"Expected order: {', '.join(expected_cols)}.\n"
                         f"Actual order:   "
                         f"{', '.join(df.columns[:len(expected_cols)])}")
    
    return None



def is_preserved_character(char):
    """True for non-ASCII characters kept in output: letters in
    PRESERVED_LETTER_RANGES and PRESERVED_SYMBOLS."""

    if char in PRESERVED_SYMBOLS:
        return True
    if not unicodedata.category(char).startswith('L'):
        return False
    code = ord(char)
    return any(start <= code <= end for start, end in PRESERVED_LETTER_RANGES)



# Mapping of specific non-standard characters to replacements
CHARACTER_REPLACEMENTS = {
        '\"'    : '\'',     # STANDARD double quotes - replace with single quote
        '\u2028': ' ',      # Line Separator (LS) - replace with space
        '\u2029': ' ',      # Paragraph Separator (PS) - replace with space
        '\u00a0': ' ',      # Non-breaking space - replace with space
        '\ufeff': '',       # Byte Order Mark (BOM) warning sign - remove
        '\ufffd': '?',      # Unicode replacement char (lost upstream) - placeholder
        '\u2013': '-',      # En Dash - replace with dash
        '\u2014': '-',      # Em Dash - replace with dash
        '\u2010': '-',      # Hyphen - replace with dash
        '\u201c': '\'',     # Left double quotes - replace with single quote
        '\u201d': '\'',     # Right double quotes - replace with single quote
        '\u2018': '\'',     # Left single quote - replace with single quote
        '\u2019': '\'',     # Right single quote - replace with single quote
        '\u0004': '',       # Unused unicode - remove
        '\u0005': '',       # Unused unicode - remove
        '\u0006': '',       # Unused unicode - remove
        '\u0007': '',       # Unused unicode - remove
        '\x02': '-',        # Start of Text control char - seen in place of a lost hyphen
        '\u03b1': 'a',      # Greek alpha - replace with a
        '\u03b2': 'B',      # Greek beta - replace with B
        '\u03b3': 'g',      # Greek gamma - replace with g
        '\u03b4': 'd',      # Greek lower delta
        '\u03F5': 'e',      # Greek lower epsilon
        '\u03ba': 'k',      # Greek kappa - replace with k
        '\u03bc': 'u',      # Greek mu - replace with u
        '\u03BB': 'l',      # Greek lambda - replace with l
        '\u0394': 'D',      # Greek upper delta
        '\u2070': '^0',     # Superscript 0
        '\u00B9': '^1',     # Superscript 1
        '\u00B2': '^2',     # Superscript 2
        '\u00B3': '^3',     # Superscript 3
        '\u2074': '^4',     # Superscript 4
        '\u2075': '^5',     # Superscript 5
        '\u2076': '^6',     # Superscript 6
        '\u2077': '^7',     # Superscript 7
        '\u2078': '^8',     # Superscript 8
        '\u2079': '^9',     # Superscript 9
        '\u2081': '1',      # Subscript 1
        '\u2082': '2',      # Subscript 2
        '\u2083': '3',      # Subscript 3
        '\u2084': '4',      # Subscript 4
        '\u2085': '5',      # Subscript 5
        '\u2086': '6',      # Subscript 6
        '\u2087': '7',      # Subscript 7
        '\u2088': '8',      # Subscript 8
        '\u2089': '9',      # Subscript 9
        '\u00a9': '(C)',    # Copywrite symbol
        '\u2264': '<=',     # Less than or equal to
        '\u2265': '>=',     # Greater than or equal to

        # Add mappings here in the future as needed
        }

# Remove all other control characters; tabs and line breaks are collapsed to
# spaces by TAB_LINE_BREAK_PATTERN instead
CHARACTER_REPLACEMENTS.update({
    chr(code): ''
    for code in [*range(0x00, 0x09), *range(0x0e, 0x20), 0x7f]
    if chr(code) not in CHARACTER_REPLACEMENTS})

# Tabs and line breaks, with any surrounding spaces, become a single space
TAB_LINE_BREAK_PATTERN = re.compile(r' *[\t\n\r\x0b\x0c][\t\n\r\x0b\x0c ]*')

ASCII_TRANSLATION_TABLE = str.maketrans(CHARACTER_REPLACEMENTS)
PRESERVING_TRANSLATION_TABLE = str.maketrans(
    {char: replacement for char, replacement in CHARACTER_REPLACEMENTS.items()
     if not is_preserved_character(char)})



def replace_defined_characters(text, preserve_characters=False):
    """Use a mappping to replace specific non-standard characters in text.
    With preserve_characters, accented/Greek letters and PRESERVED_SYMBOLS
    are left in place instead of being mapped to ASCII."""

    # Return NaN and non-string values as-is
    if pd.isna(text) or not isinstance(text, str):
        return text

    text = TAB_LINE_BREAK_PATTERN.sub(' ', text)
    if preserve_characters:
        return text.translate(PRESERVING_TRANSLATION_TABLE)
    return text.translate(ASCII_TRANSLATION_TABLE)



def normalize_encoding(text, preserve_characters=False):
    """Normalizes text encoding for consistency and compatibility. By default
    reduces text to ASCII. With preserve_characters, keeps accented/Greek
    letters and PRESERVED_SYMBOLS, reducing only other characters to ASCII."""

    # Return NaN and non-string values as-is
    if pd.isna(text) or not isinstance(text, str):
        return text

    if preserve_characters:
        # NFC joins base letter + combining accent into one character (e.g. í)
        text = unicodedata.normalize('NFC', text)
        return ''.join(
            char if ord(char) < 128 or is_preserved_character(char)
            else (unicodedata.normalize('NFD', char)
                  .encode('ascii', 'ignore').decode('utf-8'))
            for char in text)

    # Encode and decode text through an ascii "filter" to normalize
    try:
        # Using NFD will decompose characters into base and combining accent
        formatted_text = (unicodedata.normalize('NFD', text)
                # 'ignore' argument will drop any non-ascii, including accents
                          .encode('ascii', 'ignore').decode('utf-8'))
        return formatted_text
    
    # Second catch to handle non-string values as-is
    except TypeError:
        return text

    

def process_special_characters(df, column_configs, datatype):
    """
    Cleans a DataFrame by normalizing special characters and encoding.
    Excludes fields listed in config as 'exclude_special_char_processing'.
    ASCII_ONLY_COLS are reduced to ASCII; other fields keep accented/Greek
    letters and PRESERVED_SYMBOLS.

    Args:
        df (pd.DataFrame): Input DataFrame.
        column_configs (dict, optional): Configuration dictionary.
        datatype (str, optional): Data type key for column_configs.

    Returns:
        pd.DataFrame: Cleaned DataFrame.
    """

    # Get any predefined excluded fields from config
    exclude_cols = (column_configs[datatype]
                    .get('exclude_special_char_processing', []))

    # Use a copy for cleaner changes
    df_cleaned = df.copy()

    # Iterate through non-excluded fields
    cols_to_process = [col for col in df_cleaned.columns 
                       if col not in exclude_cols]
    
    for col in cols_to_process:
        preserve = col not in ASCII_ONLY_COLS
        df_cleaned[col] = (df_cleaned[col]
                            # Replace specific non-standard characters
                           .map(lambda x: replace_defined_characters(x, preserve))
                           # General normalization to utf-8 encoding
                           .map(lambda x: normalize_encoding(x, preserve)))

    return df_cleaned



def remove_html_tags(text):
    """Removes HTML tags using regular expressions."""

    # Matches opening or closing tags with tag name like <i>italics</i>
    tag_pattern = r'<(/?)(\w+)>'

    return re.sub(tag_pattern, '', text)



def remove_html_tags_from_df(df, column_configs, datatype):
    """Removes HTML tags such as <i>italics</i> or <b>bold</b> from text."""

    # Get predefined columns containing HTML tags from configuration
    html_tag_cols = column_configs[datatype].get('html_tag_cols', None)

    # Iterate through expected HTML-tagged columns
    if html_tag_cols is not None:

        for col in html_tag_cols:
        
            # Check that expected column exists
            if col not in df.columns:
                raise ValueError(f"Expected html-tagged column {col} not found "
                                 f"within data. Check data or redefine config.")
        
            # Use regex to remove HTML tags but keep enclosed text
            df[col] = df[col].apply(remove_html_tags)

    return df



def decode_html_reference(match):
    """Decode one numeric or complete named HTML reference. Unknown names
    are left as-is."""

    reference = match.group()
    if reference[1] == '#':
        return html.unescape(reference)
    return HTML5_ENTITIES.get(reference[1:], reference)



def unescape_html_entities(text):
    """Decode HTML references until unchanged so nested ones like &amp;#945;
    fully decode to α. Only numeric and complete named (;-terminated)
    references are decoded. Non-string values are returned as-is."""

    if not isinstance(text, str):
        return text
    decoded = HTML_REFERENCE_PATTERN.sub(decode_html_reference, text)
    while decoded != text:
        text = decoded
        decoded = HTML_REFERENCE_PATTERN.sub(decode_html_reference, text)
    return decoded



def clean_html_entities(df, datatype):
    """Replaces HTML character entities (&#8217) with unicode characters.
    Must run before process_special_characters so decoded characters get the
    same character handling as other text."""
    
    exclude_cols = []

    # Retain html encoding and characters for dataset display
    if datatype in DATASET_TYPES:
        exclude_cols = DATASET_HTML_ENTITY_EXCLUDE_COLS

    for col in df.columns:
        # pandas 3 text columns use StringDtype rather than object
        if col not in exclude_cols and (
                pd.api.types.is_object_dtype(df[col])
                or pd.api.types.is_string_dtype(df[col])):
            # Replace htmls with unicode characters
            df[col] = df[col].map(unescape_html_entities)

    return df



def format_datetime_columns(df, column_configs, datatype):
    """Converts datetime columns to "yyyy-mm-dd" strings."""

    # Get predefined list-like columns from configuration
    datetime_cols = column_configs[datatype].get('datetime_cols', None)

    # Iterate through expected datetime-like columns
    if datetime_cols is not None:

        for col in datetime_cols:

            # Check that expected column exists
            if col not in df.columns:
                raise ValueError(f"Expected datetime column {col} not found "
                                 f"within data. Check data or redefine config.")
            
            # Convert datetime string values to datetime and back to string
            format_string = '%Y-%m-%d'
            df[col] = (pd.to_datetime(df[col], utc=True)
                       .dt.strftime(format_string))

    return df



def validate_listlike_columns(df, column_configs, datatype):
    """Validate semicolon separation without spaces in list-like columns."""

    # Get predefined list-like columns from configuration
    listlike_cols = column_configs[datatype].get('list_like_cols', None)

    # Iterate through expected list-like columns
    if listlike_cols is not None:

        for col in listlike_cols:

            # Check that expected column exists
            if col not in df.columns:
                raise ValueError(f"Expected list-like column {col} not found "
                                 f"within data. Check data or redefine config.")

            # Skip current iteration if col contains non-string (None) values
            if not pd.api.types.is_string_dtype(df[col]):
                print(f"Column '{col}' contains non-string values. Skipping "
                      f"validation of list-like field.")
                continue

            # Replace any semicolons followed by whitespace with just semicolon
            df[col] = df[col].str.replace('; ', ';')

            # Check if any semicolons exist in any value of the column
            if not df[col].astype(str).str.contains(';').any():
                print(f"WARNING: No semicolons found in list-like column {col}")

    return df



def validate_and_clean_unique_nodes(df, column_configs, datatype):
    """
    Validate that all node_ids are either unique values or valid duplicates.

    Args:
        df (pd.DataFrame): Input DataFrame
        column_configs (dict): Configuration settings for columns, including 
            node_id and link_id
        datatype (str): Type of data in the DataFrame. Indicates which info
            will be accessed from column_configs

    Returns:
        pd.DataFrame: A new DataFrame with invalid duplicate rows removed.

    Raises:
        ValueError: If the specified node_id column is not found in the DataFrame.

    Notes:
        This function checks for valid uniqueness or duplicates in node_ids 
            based on the provided configuration.

        - For many-to-many nodes, it ensures that any duplicates only differ 
            in the link_id column.
        - True duplicates (identical values in all columns) are kept in the 
            first row for production.
        - For mismatched duplicates, the first row is kept for production.
        - All duplicates (including kept first occurrences) are recorded in the 
            invalid_duplicates_df.
        - Rows with the same node_id but different link_ids are not considered 
            duplicates and are all kept.
    """

    # Get node_id and link_id from column configurations
    node_id = column_configs[datatype].get('node_id')
    link_id = column_configs[datatype].get('link_id', None)

    # Check if the specified node_id column exists
    if node_id not in df.columns:
        raise ValueError(f"Node ID column '{node_id}' not found in DataFrame.")

    # Define duplicate check columns
    dup_check_cols = [node_id, link_id] if link_id else [node_id]

    # Get only rows with duplicated node_ids and link_ids (if applicable)
    dup_mask = df.duplicated(subset=dup_check_cols, keep=False)
    dup_nodes = df[dup_mask]
    
    # Get true duplicates where values in all columns are identical
    true_duplicates = dup_nodes[dup_nodes.duplicated(keep=False)].assign(
        error_reason=f"True duplicate in all columns. "
                     f"Kept first row for production.")
    
    # Empty dataframe to store mismatch rows
    mismatch_df = pd.DataFrame()

    # Get list of value columns to check for mismatch
    value_cols = list(set(df.columns) - set(dup_check_cols))

    # Group rows by node_id (and link_id if applicable) for shared value checking
    for node, node_df in dup_nodes.groupby(dup_check_cols):
        
        # Check for mismatched values in other columns within the same group
        value_mismatches = node_df[~node_df.duplicated(
            subset=value_cols, keep=False)
            ].assign(error_reason=f"Mismatched values in columns with same "
                                  f"node_id (and link_id if applicable). "
                                  f"Kept first row for production.")

        # Add all mismatched rows to the mismatch_df
        mismatch_df = pd.concat([mismatch_df, value_mismatches])

    # Combine list of annotated invalid duplicates for export to csv
    invalid_duplicates_df = pd.concat([true_duplicates, mismatch_df])

    # Save any invalid duplicate rows to a report csv for troubleshooting
    if len(invalid_duplicates_df) > 0:
        invalid_duplicates_path = f"{config.REMOVED_DUPLICATES}{datatype}.csv"
        os.makedirs(os.path.dirname(invalid_duplicates_path), exist_ok=True)
        invalid_duplicates_df.to_csv(invalid_duplicates_path, index=False)
        print(f"Invalid duplicate rows (including kept first rows) recorded in "
              f"{invalid_duplicates_path}.")

    # Final df with only unique combos of node_id and link_id (if applicable)
    valid_df = df.drop_duplicates(subset=dup_check_cols, keep='first'
                                  ).reset_index(drop=True)

    return valid_df



def enforce_int_values(df, column_configs, datatype):
    """Enforce integer values in expected columns."""
    
    # Get predefined columns from configuration
    int_cols = column_configs[datatype].get('int_cols', None)
    
    # Iterate through expected columns
    if int_cols is not None:
        for col in int_cols:
            
            # Check that expected column exists
            if col not in df.columns:
                raise ValueError(f"Expected integer column {col} not found "
                               f"within data. Check data or redefine config.")
                
            # Convert column to Int64. Keeps NaN and blank vals
            df[col] = df[col].astype('Int64')

    return df



def remove_nan_strings(df):
    """Remove 'nan' string values and replace with blank space."""

    df = df.map(lambda x: '' if isinstance(x, str) and x.lower() == 'nan' else x)

    return df



def standardize_data(df, column_configs, datatype):
    """Group standardization functions common to all data types"""

    # Edit data to standardize
    df = add_type_column(df, datatype)
    df = reorder_columns(df, column_configs, datatype)
    df = clean_html_entities(df, datatype)
    df = process_special_characters(df, column_configs, datatype)
    df = remove_html_tags_from_df(df, column_configs, datatype)
    df = validate_listlike_columns(df, column_configs, datatype)
    df = format_datetime_columns(df, column_configs, datatype)
    df = validate_and_clean_unique_nodes(df, column_configs, datatype)
    df = remove_nan_strings(df)
    df = enforce_int_values(df, column_configs, datatype)

    # Validate that data meets loading standards
    validate_first_columns(df, column_configs, datatype)
    
    return df



def package_programs(df_programs, column_configs):
    """Package programs data for INS loading."""

    print(f"---\nFinalizing TSV for program data...")

    # Standardize and validate data
    df_programs_output = standardize_data(df_programs, column_configs, 
                                        datatype='program')

    # Export as TSV
    output_filepath = config.PROGRAMS_OUTPUT_PATH
    os.makedirs(os.path.dirname(output_filepath), exist_ok=True)
    df_programs_output.to_csv(output_filepath, sep='\t', index=False, 
                            encoding='utf-8')

    print(f"Done! Final program data saved as {output_filepath}.")

    return df_programs_output



def package_grants(df_grants, column_configs):
    """Package grants data for INS loading."""

    print(f"---\nFinalizing TSV for grant data...")

    # Standardize and validate data
    df_grants_output = standardize_data(df_grants, column_configs, 
                                        datatype='grant')

    # Export as TSV
    output_filepath = config.GRANTS_OUTPUT_PATH
    os.makedirs(os.path.dirname(output_filepath), exist_ok=True)
    df_grants_output.to_csv(output_filepath, sep='\t', index=False, 
                            encoding='utf-8')

    print(f"Done! Final grant data saved as {output_filepath}.")

    return df_grants_output



def package_projects(df_projects, column_configs):
    """Package projects data for INS loading."""

    print(f"---\nFinalizing TSV for project data...") 

    # Standardize and validate data
    df_projects_output = standardize_data(df_projects, column_configs, 
                                        datatype='project')

    # Export as TSV
    output_filepath = config.PROJECTS_OUTPUT_PATH
    os.makedirs(os.path.dirname(output_filepath), exist_ok=True)
    df_projects_output.to_csv(output_filepath, sep='\t', index=False, 
                            encoding='utf-8')

    print(f"Done! Final projects data saved as {output_filepath}.")

    return df_projects_output



def package_publications(df_publications, column_configs):
    """Package publications data for INS loading."""

    print(f"---\nFinalizing TSV for publication data...") 

    # Standardize and validate data
    df_publications_output = standardize_data(df_publications, column_configs, 
                                        datatype='publication')

    # Export as TSV
    output_filepath = config.PUBLICATIONS_OUTPUT_PATH
    os.makedirs(os.path.dirname(output_filepath), exist_ok=True)
    df_publications_output.to_csv(output_filepath, sep='\t', index=False, 
                                    encoding='utf-8')

    print(f"Done! Final publication data saved as {output_filepath}.")

    return df_publications_output



def package_dbgap_datasets(df_dbgap_datasets, column_configs, dbgap_curated=False):
    """Package dbGaP datasets data for INS loading."""

    print(f"---\nFinalizing TSV for dbGaP datasets data...") 

    # Standardize and validate data
    df_dbgap_datasets_output = standardize_data(df_dbgap_datasets, column_configs, 
                                        datatype='dbgap_dataset')

    # Export as TSV
    output_filepath = config.DBGAP_OUTPUT_PATH
    if dbgap_curated: 
        output_filepath = config.DBGAP_OUTPUT_CURATED_CLEANED
    os.makedirs(os.path.dirname(output_filepath), exist_ok=True)
    df_dbgap_datasets_output.to_csv(output_filepath, sep='\t', index=False, 
                                    encoding='utf-8')

    print(f"Done! Final dbGaP Datasets data saved as {output_filepath}.")

    return df_dbgap_datasets_output


def backfill_curated_dbgap_accessions(curated_df, gathered_path):
    """Add full accessions to a legacy curated file from gathered data.
    Likely one-time-use after addition of full accessions to the curated file."""

    accession_col = 'dataset_source_accession'
    if accession_col in curated_df.columns:
        return curated_df

    if not os.path.exists(gathered_path):
        raise FileNotFoundError(
            f"Cannot backfill '{accession_col}': gathered dbGaP file not "
            f"found at {gathered_path}")

    gathered_df = pd.read_csv(
        gathered_path, dtype=str, keep_default_na=False,
        usecols=['accession', 'full_accession'])
    accession_map = (gathered_df.drop_duplicates('accession')
                     .set_index('accession')['full_accession'])

    result_df = curated_df.copy()
    result_df[accession_col] = (
        result_df['dataset_source_id'].map(accession_map).fillna(''))

    missing_count = (result_df[accession_col] == '').sum()
    print(f"Backfilled full dbGaP accessions for "
          f"{len(result_df) - missing_count} curated rows.")
    if missing_count:
        print(f"WARNING: {missing_count} curated row(s) were not present in "
              f"the gathered file and have a blank {accession_col}.")

    return result_df



def package_geo_datasets(df_geo_datasets, column_configs):
    """Package GEO Datasets data for INS loading."""

    print(f"---\nFinalizing TSV for GEO Datasets data...") 

    # Standardize and validate data
    df_geo_datasets_output = standardize_data(df_geo_datasets, column_configs, 
                                        datatype='geo_dataset')

    # Export as TSV
    output_filepath = config.GEO_OUTPUT_PATH
    os.makedirs(os.path.dirname(output_filepath), exist_ok=True)
    df_geo_datasets_output.to_csv(output_filepath, sep='\t', index=False, 
                                    encoding='utf-8')

    print(f"Done! Final GEO Datasets data saved as {output_filepath}.")

    return df_geo_datasets_output


def package_sra_datasets(df_sra_datasets, column_configs, sra_curated=False):
    """Package SRA Datasets data for INS loading."""

    print(f"---\nFinalizing TSV for SRA Datasets data...") 

    # Standardize and validate data
    df_sra_datasets_output = standardize_data(df_sra_datasets, column_configs, 
                                        datatype='sra_dataset')

    # Export as TSV
    output_filepath = config.SRA_OUTPUT_PATH
    if sra_curated: 
        output_filepath = config.SRA_OUTPUT_CURATED_CLEANED
    os.makedirs(os.path.dirname(output_filepath), exist_ok=True)
    df_sra_datasets_output.to_csv(output_filepath, sep='\t', index=False, 
                                    encoding='utf-8')

    print(f"Done! Final SRA Datasets data saved as {output_filepath}.")

    return df_sra_datasets_output



def package_cedcd_datasets(df_cedcd_datasets, column_configs):
    """Package CEDCD Cohorts as datasets for INS loading."""

    print(f"---\nFinalizing TSV for CEDCD Datasets data...") 

    # Standardize and validate data
    df_cedcd_datasets_output = standardize_data(df_cedcd_datasets, column_configs, 
                                        datatype='cedcd_dataset')

    # Export as TSV
    output_filepath = config.CEDCD_OUTPUT_PATH
    os.makedirs(os.path.dirname(output_filepath), exist_ok=True)
    df_cedcd_datasets_output.to_csv(output_filepath, sep='\t', index=False, 
                                    encoding='utf-8')

    print(f"Done! Final cedcd Datasets data saved as {output_filepath}.")

    return df_cedcd_datasets_output



def package_dceg_datasets(df_dceg_datasets, column_configs):
    """Package DCEG Cohorts as datasets for INS loading.

    Generates deterministic UUID5 values for dataset_uuid based on
    dataset_title and dataset_source_id. This is a curated-only source
    with no gathering pipeline.
    """

    print(f"---\nFinalizing TSV for DCEG Cohorts data...")

    # Generate deterministic UUID5 for each row based on title + source ID
    NAMESPACE = uuid.UUID('12345678-1234-5678-1234-567812345678')

    def generate_dceg_uuid(row):
        name = '||'.join(str(row[field]) for field in 
                         ['dataset_title', 'dataset_source_id'])
        return uuid.uuid5(NAMESPACE, name)

    df_dceg_datasets['dataset_uuid'] = df_dceg_datasets.apply(
        generate_dceg_uuid, axis=1)

    # Standardize and validate data
    df_dceg_datasets_output = standardize_data(df_dceg_datasets, column_configs,
                                        datatype='dceg_dataset')

    # Export as TSV
    output_filepath = config.DCEG_OUTPUT_PATH
    os.makedirs(os.path.dirname(output_filepath), exist_ok=True)
    df_dceg_datasets_output.to_csv(output_filepath, sep='\t', index=False,
                                    encoding='utf-8')

    print(f"Done! Final DCEG Cohorts data saved as {output_filepath}.")

    return df_dceg_datasets_output



def package_nccr_datasets(df_nccr_datasets, column_configs):
    """Package NCCR datasets for INS loading.

    Generates deterministic UUID5 values for dataset_uuid based on
    dataset_source_id and dataset_title. This is a curated-only source
    with no gathering pipeline.
    """

    print(f"---\nFinalizing TSV for NCCR Datasets data...")

    # Generate deterministic UUID5 for each row based on source ID + title
    NAMESPACE = uuid.UUID('12345678-1234-5678-1234-567812345678')

    def generate_nccr_uuid(row):
        name = '||'.join(str(row[field]) for field in
                         ['dataset_source_id', 'dataset_title'])
        return uuid.uuid5(NAMESPACE, name)

    df_nccr_datasets['dataset_uuid'] = df_nccr_datasets.apply(
        generate_nccr_uuid, axis=1)

    # Standardize and validate data
    df_nccr_datasets_output = standardize_data(df_nccr_datasets, column_configs,
                                        datatype='nccr_dataset')

    # Export as TSV
    output_filepath = config.NCCR_OUTPUT_PATH
    os.makedirs(os.path.dirname(output_filepath), exist_ok=True)
    df_nccr_datasets_output.to_csv(output_filepath, sep='\t', index=False,
                                    encoding='utf-8')

    print(f"Done! Final NCCR Datasets data saved as {output_filepath}.")

    return df_nccr_datasets_output



def package_dctd_datasets(df_dctd_datasets, column_configs):
    """Package DCTD datasets for INS loading.

    Generates deterministic UUID5 values for dataset_uuid based on
    dataset_source_id and dataset_title. This is a curated-only source
    with no gathering pipeline.
    """

    print(f"---\nFinalizing TSV for DCTD Datasets data...")

    # Generate deterministic UUID5 for each row based on source ID + title
    NAMESPACE = uuid.UUID('12345678-1234-5678-1234-567812345678')

    def generate_dctd_uuid(row):
        name = '||'.join(str(row[field]) for field in
                         ['dataset_source_id', 'dataset_title'])
        return uuid.uuid5(NAMESPACE, name)

    df_dctd_datasets['dataset_uuid'] = df_dctd_datasets.apply(
        generate_dctd_uuid, axis=1)

    # Standardize and validate data
    df_dctd_datasets_output = standardize_data(df_dctd_datasets, column_configs,
                                        datatype='dctd_dataset')

    # Export as TSV
    output_filepath = config.DCTD_OUTPUT_PATH
    os.makedirs(os.path.dirname(output_filepath), exist_ok=True)
    df_dctd_datasets_output.to_csv(output_filepath, sep='\t', index=False,
                                    encoding='utf-8')

    print(f"Done! Final DCTD Datasets data saved as {output_filepath}.")

    return df_dctd_datasets_output



def package_ctd2_datasets(df_ctd2_datasets, column_configs, ctd2_curated=False):
    """Package CTD^2 Datasets data for INS loading."""

    print(f"---\nFinalizing TSV for CTD^2 Datasets data...") 

    if ctd2_curated:
        print(f"Locked curation of CTD^2 datasets found. Keep without changing or resaving.")
        return df_ctd2_datasets

    else: 
        # Standardize and validate data
        df_ctd2_datasets_output = standardize_data(df_ctd2_datasets, column_configs, 
                                            datatype='ctd2_dataset')

        # Export as TSV
        output_filepath = config.CTD2_DATASET_OUTPUT_PATH
        os.makedirs(os.path.dirname(output_filepath), exist_ok=True)
        df_ctd2_datasets_output.to_csv(output_filepath, sep='\t', index=False, 
                                        encoding='utf-8')

        print(f"Done! Final CTD^2 Datasets data saved as {output_filepath}.")

        return df_ctd2_datasets_output


def package_ctd2_filedata(df_ctd2_filedata, column_configs):
    """Package CTD^2 file metadata for INS loading."""

    print(f"---\nFinalizing TSV for CTD^2 file metadata...")

    # Standardize and validate data
    df_ctd2_filedata_output = standardize_data(df_ctd2_filedata, column_configs, 
                                        datatype='file')

    # Export as TSV
    output_filepath = config.CTD2_FILE_OUTPUT_PATH
    os.makedirs(os.path.dirname(output_filepath), exist_ok=True)
    df_ctd2_filedata_output.to_csv(output_filepath, sep='\t', index=False, 
                            encoding='utf-8')

    print(f"Done! Final CTD^2 file metadata saved as {output_filepath}.")

    return df_ctd2_filedata_output



def remove_publications_before_projects(df_publications: pd.DataFrame, 
                                        df_projects: pd.DataFrame, 
                                        day_diff_allowed:int=0) -> pd.DataFrame:
    """Remove publications published before the associated project start date.
    
    Args:
        df_publications (pd.DataFrame): DataFrame containing publication data.
        df_projects (pd.DataFrame): DataFrame containing project information.
        day_diff_allowed (int, optional): Allowable difference in days between
            publication date and project start date. Defaults to 0. 
            If None, then no filtering is applied and original publications
            are returned.

    Returns:
        pd.DataFrame: DataFrame containing publications with valid dates. 
            If day_diff_allowed == None, returns original df_publications.

    Raises:
        ValueError: If column structure of the filtered DataFrame is altered

    Notes:
        Exports the DataFrame of removed publications as CSV to a location
            defined in config.py.
    """

    # Skip this function entirely if None provided for day_diff_allowed 
    if not day_diff_allowed:
        print(f"Filter NOT applied for publications published before project "
              f"start date. Keeping all. Change this in config.py if desired.")
        return df_publications

    # Validate day difference and convert to pd.Timedelta value
    if not isinstance(day_diff_allowed, int):
        raise ValueError(f"Expected int value for `day_diff_allowed` but found "
                         f"{type(day_diff_allowed)}. Fix this in config.py.")
    diff = str(day_diff_allowed) + ' days'

    # Merge project start dates with publication data
    df_merged = pd.merge(left=df_publications, 
                        right=df_projects[['project_id','project_start_date']], 
                        how='left',left_on='coreproject', right_on='project_id',
                        # Ignore duplicate projects listed for multiple programs
                        ).drop_duplicates()

    # Get subset where publication pub date is earlier than project start date
    df_early_pubs = df_merged[
        # Convert dates to datetime for more accurate comparison
        (pd.to_datetime(df_merged['publication_date'],utc=True) 
         # Add in allowable day difference
        + pd.Timedelta(diff))
        < pd.to_datetime(df_merged['project_start_date'],utc=True)
        # Set on a copy to avoid SettingWithCopyWarning
        ].copy()

    # Get publications to keep after early pubs removed
    df_keep_pubs = df_merged[~df_merged.index.isin(df_early_pubs.index)]

    # Drop merged columns and reindex filtered publications
    df_keep_pubs = df_keep_pubs.drop(columns=['project_id','project_start_date']
                                     ).reset_index(drop=True)
    
    # Add column showing difference between pub date and project start date
    df_early_pubs.loc[:,'day_diff'] = (
        pd.to_datetime(df_early_pubs['project_start_date'],utc=True) -
        pd.to_datetime(df_early_pubs['publication_date'],utc=True)
        ).dt.days

    
    if not df_early_pubs.empty:

        # Define and make directory for report
        early_pub_report_path = config.REMOVED_EARLY_PUBLICATIONS
        os.makedirs(os.path.dirname(early_pub_report_path), exist_ok=True)

        # Export report of removed early publications
        df_early_pubs.to_csv(early_pub_report_path, index=False)
        print(f"Early Publications detected:\n"
              f"{len(df_early_pubs)} Publications with publication date more "
              f"than {diff} before the associated project start date were "
              f"removed and saved to {early_pub_report_path}")
        
    # Validate that columns have not changed in returned filtered list
    if df_publications.columns.tolist() != df_keep_pubs.columns.tolist():
        raise ValueError(f"Unexpected change in columns during "
                         f"`remove_publications_before_projects()`.\n"
                         f"Expected: {df_publications.columns.tolist()}\n"
                         f"Actual:   {df_keep_pubs.columns.tolist()}")

    return df_keep_pubs



def get_enum_values(df, cols, output="enums.txt"):
    """Extracts unique values from specified columns in a DataFrame,
    sorts them, and writes them to a YAML-like text file.

    Args:
        df (pd.DataFrame): The DataFrame to process.
        cols (list): A list of column names to extract unique values from.
        output (str, optional): The output file path. Defaults to "enums.txt".
    """

    # Make directory if it doesn't already exist
    os.makedirs(os.path.dirname(output), exist_ok=True)

    # Get data type
    df_type = df['type'][0]

    # Build empty dict
    enums_dict = {}

    for col in cols:
        
        # Skip current iteration if col contains non-string (None) values
        if not pd.api.types.is_string_dtype(df[col]):
            continue

        # Get unique values from semicolon-separated strings
        unique_vals = df[col].dropna().str.split(";").explode().unique()
        # Sort for conistency
        unique_vals = sorted(list(unique_vals))
        # Add "Enum" header to match data model yaml structure
        enums_dict[col] = {"Enum": unique_vals}

    # Save as txt with yaml-like formatting 
    with open(output, "w") as f:
        # Use 2-space indents
        indent = "  "

        # Manually imitate formatting of yaml structure and get data type
        f.write(f"PropDefinitions:\n"
                f"{indent}#properties of {df_type}\n"
                f"{indent}...\n")
        
        # Add columns as headers
        for col, enum_data in enums_dict.items():
            f.write(f"{indent}{col}:\n"
                        f"{indent*2}... \n"
                        f"{indent*2}Type:\n"
                                f"{indent*3}value_type: list\n"
                                f"{indent*3}Enum:\n")
            
            for val in enum_data["Enum"]:
                            # Iterate through unique values as enums
                            f.write(f"{indent*4}- {val}\n")
            # Add trailing ... to end of section
            f.write(f"{indent*2}...\n")

    print(f"Done! Enumerated values for {df_type}s: {cols} saved to {output}")



def generate_md5_hash(file_path):
    """Generates the MD5 hash for a given file."""

    with open(file_path, "rb") as f:
        data = f.read()
        md5_hash = hashlib.md5(data).hexdigest()

    return md5_hash



def generate_md5_hash_file(directory):
    """Generate MD5s for all TSVs in a directory and list in a text file."""

    hashes_file = os.path.join(directory, "_md5.txt")

    with open(hashes_file, "w") as f:
        for root, _, files in os.walk(directory):
            for file in files:
                if file.endswith(".tsv"):
                    file_path = os.path.join(root, file)
                    md5_hash = generate_md5_hash(file_path)
                    f.write(f"{md5_hash}\t{file_path}\n")
    
    print(f"Done! MD5 hashes saved to {hashes_file}.")



def scan_output_tsv(filepath, datatype, column_configs):
    """Scan one packaged TSV for special character problems and return a
    list of finding dicts, one per value (or per file for file_encoding)."""

    datatype_config = column_configs.get(datatype, {})
    exclude_cols = datatype_config.get('exclude_special_char_processing', [])
    html_tag_cols = datatype_config.get('html_tag_cols') or []
    entity_allowed_cols = (DATASET_HTML_ENTITY_EXCLUDE_COLS
                           if datatype in DATASET_TYPES else [])

    findings = []
    record_ids = None

    def add_finding(check, idx, column, detail, snippet):
        findings.append({
            'file': str(filepath),
            'datatype': datatype,
            'record_id': record_ids[idx] if idx is not None else None,
            'check': check,
            'severity': OUTPUT_QA_CHECKS[check][0],
            # Header is Excel row 1, so first data row is row 2
            'excel_row': idx + 2 if idx is not None else None,
            'column': column,
            'detail': detail,
            'snippet': snippet,
        })

    with open(filepath, 'rb') as f:
        raw = f.read()
    if raw.startswith(b'\xef\xbb\xbf'):
        add_finding('file_encoding', None, None, 'Starts with UTF-8 BOM', '')
    try:
        raw.decode('utf-8')
    except UnicodeDecodeError as e:
        add_finding('file_encoding', None, None, f'Not valid UTF-8: {e}', '')

    df = pd.read_csv(filepath, sep='\t', dtype=str, keep_default_na=False,
                     encoding='utf-8-sig', encoding_errors='replace')

    # Dataset source IDs (e.g. phs000001) are easier to search than UUIDs
    id_col = ('dataset_source_id' if 'dataset_source_id' in df.columns
              else datatype_config.get('node_id'))
    record_ids = (df[id_col] if id_col in df.columns
                  else pd.Series('', index=df.index))

    for col in df.columns:
        values = df[col]
        pattern_checks = [('replacement_char', re.compile('\ufffd')),
                          ('garbled_text', MOJIBAKE_PATTERN),
                          ('control_character', CONTROL_CHAR_PATTERN),
                          ('placeholder_char', PLACEHOLDER_PATTERN)]
        if col not in entity_allowed_cols:
            pattern_checks.append(('html_entity', HTML_ENTITY_PATTERN))
        if col in html_tag_cols:
            pattern_checks.append(('html_tag', HTML_TAG_PATTERN))

        for check, pattern in pattern_checks:
            for idx, value in values[values.str.contains(pattern)].items():
                text = (URL_PATTERN.sub(' ', value)
                        if check == 'placeholder_char' else value)
                match = pattern.search(text)
                if match:
                    snippet = text[max(0, match.start() - 30):match.end() + 30]
                    add_finding(check, idx, col, repr(match.group()), snippet)

        if col in exclude_cols:
            continue
        ascii_only = col in ASCII_ONLY_COLS
        for idx, value in values[values.str.contains(NON_ASCII_PATTERN)].items():
            unexpected = sorted({
                char for char in NON_ASCII_PATTERN.findall(value)
                if char != '\ufffd'
                and (ascii_only or not is_preserved_character(char))})
            if unexpected:
                first = value.find(unexpected[0])
                detail = ', '.join(
                    f"{char} U+{ord(char):04X} "
                    f"{unicodedata.name(char, 'UNNAMED')}" for char in unexpected)
                add_finding('unexpected_character', idx, col, detail,
                            value[max(0, first - 30):first + 30])

    return findings



def run_output_character_qa(output_files, column_configs, summary_path,
                            details_path):
    """Scan all packaged TSVs, save a per-file/per-check summary CSV and a
    detail CSV with every finding, and print an overall PASS/FAIL result.
    Reports are only written when there is at least one finding.

    Args:
        output_files (list): (datatype, filepath) tuples of packaged TSVs.
        column_configs (dict): Column configuration from config.py.
        summary_path (str): Filepath for the summary CSV.
        details_path (str): Filepath for the detail CSV.

    Returns:
        tuple: (summary_df, details_df)
    """

    findings = []
    summary_rows = []
    for datatype, filepath in output_files:
        file_findings = scan_output_tsv(filepath, datatype, column_configs)
        findings.extend(file_findings)
        for check, (severity, question) in OUTPUT_QA_CHECKS.items():
            count = sum(f['check'] == check for f in file_findings)
            summary_rows.append({
                'file': str(filepath),
                'datatype': datatype,
                'check': check,
                'question': question,
                'status': severity if count else 'PASS',
                'count': count,
            })

    summary_df = pd.DataFrame(summary_rows, columns=[
        'file', 'datatype', 'check', 'question', 'status', 'count'])
    details_df = pd.DataFrame(findings, columns=[
        'file', 'datatype', 'record_id', 'check', 'severity', 'excel_row',
        'column', 'detail', 'snippet'])

    if findings:
        for path, df in [(summary_path, summary_df),
                         (details_path, details_df)]:
            os.makedirs(os.path.dirname(path), exist_ok=True)
            # BOM so snippets with accents/symbols display correctly in Excel
            df.to_csv(path, index=False, encoding='utf-8-sig')

    print(f"Output character QA across {len(output_files)} file(s):")
    for check, (severity, question) in OUTPUT_QA_CHECKS.items():
        check_rows = summary_df[summary_df['check'] == check]
        total = int(check_rows['count'].sum())
        status = severity if total else 'PASS'
        print(f"    [{status:<4}] {question} ({total})")
        for _, row in check_rows[check_rows['count'] > 0].iterrows():
            print(f"           {row['count']} in {row['file']}")

    failed = (summary_df['status'] == 'FAIL').any()
    print(f"Overall: {'FAIL - review ' + details_path if failed else 'PASS'}")
    if findings:
        print(f"Summary saved to {summary_path}")
        print(f"Details saved to {details_path}")
    else:
        print(f"No findings. Summary and details reports not written.")

    return summary_df, details_df



def package_output_data():
    """Run all data packaging steps for all data types."""

    print(f"\n---\nDATA PACKAGING:\n"
          f"Performing final data packaging steps...\n---\n")

    # Single data model dict defined in config
    column_configs = config.COLUMN_CONFIGS 

    # Special character findings from all curated files loaded below
    special_char_findings = []

    # Load programs data
    if os.path.exists(config.PROGRAMS_INTERMED_PATH):
        programs_exist = True
        df_programs = pd.read_csv(config.PROGRAMS_INTERMED_PATH)
        print(f"Loaded Programs file from {config.PROGRAMS_INTERMED_PATH}")
    else: 
        programs_exist = False
        print(f"No Programs file found.")
    
    # Load grants data
    if os.path.exists(config.GRANTS_INTERMED_PATH):
        grants_exist = True
        df_grants = pd.read_csv(config.GRANTS_INTERMED_PATH)
        print(f"Loaded Grants file from {config.GRANTS_INTERMED_PATH}")
    else: 
        grants_exist = False
        print(f"No Grants file found.")

    # Load projects data
    if os.path.exists(config.PROJECTS_INTERMED_PATH):
        projects_exist = True
        df_projects = pd.read_csv(config.PROJECTS_INTERMED_PATH)
        print(f"Loaded Projects file from {config.PROJECTS_INTERMED_PATH}")
    else: 
        projects_exist = False
        print(f"No Projects file found.")

    # Load publications data
    if os.path.exists(config.PUBLICATIONS_INTERMED_PATH):
        publications_exist = True
        df_publications = pd.read_csv(config.PUBLICATIONS_INTERMED_PATH)
        print(f"Loaded Publications file from {config.PUBLICATIONS_INTERMED_PATH}")
    else: 
        publications_exist = False
        print(f"No Publications file found.")

    # Load dbGaP datasets data. Prefer curated file if it exists.
    if os.path.exists(config.DBGAP_CURATED_INTERMED_PATH):
        dbgap_datasets_exist = True
        dbgap_curated = True
        df_dbgap_datasets = read_csv_with_encoding_fallback(
            config.DBGAP_CURATED_INTERMED_PATH, sep='\t',
            findings=special_char_findings)
        df_dbgap_datasets = backfill_curated_dbgap_accessions(
            df_dbgap_datasets, config.DBGAP_INTERMED_PATH)
        print(f"Loaded dbGaP Datasets (Curated) file from {config.DBGAP_CURATED_INTERMED_PATH}")
    elif os.path.exists(config.DBGAP_INTERMED_PATH):
        dbgap_datasets_exist = True
        dbgap_curated = False
        df_dbgap_datasets = pd.read_csv(config.DBGAP_INTERMED_PATH)
        print(f"Loaded dbGaP Datasets file from {config.DBGAP_INTERMED_PATH}")
    else: 
        dbgap_datasets_exist = False
        dbgap_curated = False
        print(f"No dbGaP Datasets file found.")

    # Load GEO datasets data
    if os.path.exists(config.GEO_INTERMED_PATH):
        geo_datasets_exist = True
        df_geo_datasets = pd.read_csv(config.GEO_INTERMED_PATH,
                                      dtype={'dataset_pmid':str})
        print(f"Loaded GEO Datasets file from {config.GEO_INTERMED_PATH}")
    else: 
        geo_datasets_exist = False
        print(f"No GEO Datasets file found.")

    # Load SRA datasets data. Prefer curated file if it exists.
    if os.path.exists(config.SRA_CURATED_INTERMED_PATH):
        sra_datasets_exist = True
        sra_curated = True
        df_sra_datasets = read_csv_with_encoding_fallback(
            config.SRA_CURATED_INTERMED_PATH, sep='\t', dtype={'dataset_pmid':str},
            findings=special_char_findings)
        print(f"Loaded SRA Datasets (Curated) file from {config.SRA_CURATED_INTERMED_PATH}")
    elif os.path.exists(config.SRA_INTERMED_PATH):
        sra_datasets_exist = True
        sra_curated = False
        df_sra_datasets = pd.read_csv(config.SRA_INTERMED_PATH,
                                      dtype={'dataset_pmid':str})
        print(f"Loaded SRA Datasets file from {config.SRA_INTERMED_PATH}")
    else: 
        sra_datasets_exist = False
        sra_curated = False
        print(f"No SRA Datasets file found.")

    # Load CEDCD datasets data
    if os.path.exists(config.CEDCD_INTERMED_CSV):
        cedcd_datasets_exist = True
        df_cedcd_datasets = pd.read_csv(config.CEDCD_INTERMED_CSV)
        print(f"Loaded CEDCD Datasets file from {config.CEDCD_INTERMED_CSV}")
    else: 
        cedcd_datasets_exist = False
        print(f"No CEDCD Datasets file found.")

    # Load DCEG Cohorts curated data
    if os.path.exists(config.DCEG_CURATED_INTERMED_PATH):
        dceg_datasets_exist = True
        df_dceg_datasets = read_csv_with_encoding_fallback(
            config.DCEG_CURATED_INTERMED_PATH, sep='\t',
            findings=special_char_findings)
        print(f"Loaded DCEG Cohorts file from {config.DCEG_CURATED_INTERMED_PATH}")
    else:
        dceg_datasets_exist = False
        print(f"No DCEG Cohorts file found.")

    # Load NCCR curated data
    if os.path.exists(config.NCCR_CURATED_INTERMED_PATH):
        nccr_datasets_exist = True
        df_nccr_datasets = read_csv_with_encoding_fallback(
            config.NCCR_CURATED_INTERMED_PATH, sep='\t',
            findings=special_char_findings)
        print(f"Loaded NCCR Datasets file from {config.NCCR_CURATED_INTERMED_PATH}")
    else:
        nccr_datasets_exist = False
        print(f"No NCCR Datasets file found.")

    # Load DCTD curated data
    if os.path.exists(config.DCTD_CURATED_INTERMED_PATH):
        dctd_datasets_exist = True
        df_dctd_datasets = read_csv_with_encoding_fallback(
            config.DCTD_CURATED_INTERMED_PATH, sep='\t',
            findings=special_char_findings)
        print(f"Loaded DCTD Datasets file from {config.DCTD_CURATED_INTERMED_PATH}")
    else:
        dctd_datasets_exist = False
        print(f"No DCTD Datasets file found.")

    # Load CTD2 datasets data
    if os.path.exists(config.CTD2_DATASET_CURATED_LOCKED_PATH):
        ctd2_datasets_exist = True
        ctd2_curated = True
        df_ctd2_datasets = read_csv_with_encoding_fallback(
            config.CTD2_DATASET_CURATED_LOCKED_PATH, sep='\t',
            findings=special_char_findings)
        print(f"Loaded CTD2 Datasets (Locked Curation) file from {config.CTD2_DATASET_CURATED_LOCKED_PATH}")
    elif os.path.exists(config.CTD2_DATASET_INTERMED_CSV):
        ctd2_datasets_exist = True
        ctd2_curated = False
        df_ctd2_datasets = pd.read_csv(config.CTD2_DATASET_INTERMED_CSV)
        print(f"Loaded CTD2 Datasets file from {config.CTD2_DATASET_INTERMED_CSV}")
    else:
        ctd2_datasets_exist = False
        ctd2_curated = False
        print(f"No CTD2 Datasets file found.")

    # Load CTD2 file metadata
    if os.path.exists(config.CTD2_FILE_INTERMED_CSV):
        ctd2_filedata_exist = True
        df_ctd2_filedata = pd.read_csv(config.CTD2_FILE_INTERMED_CSV)
        print(f"Loaded CTD2 file metadata from {config.CTD2_FILE_INTERMED_CSV}")
    else:
        ctd2_filedata_exist = False
        print(f"No CTD2 file metadata found.")

    write_special_character_scan(special_char_findings,
                                 config.SPECIAL_CHARACTER_SCAN_PATH)

    # Special handling
    print(f"\n\nApplying special handling steps...")

    print(f"---\nChecking for publications published before projects...")
    # Remove publications published before associated projects
    if publications_exist and projects_exist:
        df_publications = remove_publications_before_projects(
                                df_publications,
                                df_projects,
                                day_diff_allowed=config.PUB_PROJECT_DAY_DIFF)
    
    print(f"---\nAdding dataset UUIDs to curated dbGaP datasets if not present...")
    # Add dbgap dataset uuids for any that do not exist
    # Uses deterministic UUID5 based on source repo + source id
    if dbgap_curated:
        DBGAP_UUID_NAMESPACE = uuid.UUID('12345678-1234-5678-1234-567812345678')
        df_dbgap_datasets['dataset_uuid'] = df_dbgap_datasets.apply(
            lambda row: uuid.uuid5(
                DBGAP_UUID_NAMESPACE,
                '||'.join([str(row['dataset_source_repo']),
                           str(row['dataset_source_id'])]))
            if pd.isna(row['dataset_uuid'])
            or str(row['dataset_uuid']).lower() == 'na'
            or str(row['dataset_uuid']) == ''
            else row['dataset_uuid'],
            axis=1,
        )

    print(f"\n\nFinal packaging steps...")

    # Final packaging
    if programs_exist:
        df_programs_out = package_programs(df_programs, column_configs)
    if grants_exist:
        df_grants_out = package_grants(df_grants, column_configs)
    if projects_exist:
        df_projects_out = package_projects(df_projects, column_configs)
    if publications_exist:
        df_publications_out = package_publications(df_publications, column_configs)
    if dbgap_datasets_exist:
        df_dbgap_datasets_out = package_dbgap_datasets(df_dbgap_datasets, column_configs,
                                                       dbgap_curated)
    if geo_datasets_exist:
        df_geo_datasets_out = package_geo_datasets(df_geo_datasets, column_configs)
    if sra_datasets_exist:
        df_sra_datasets_out = package_sra_datasets(df_sra_datasets, column_configs,
                                                    sra_curated)
    if cedcd_datasets_exist:
        df_cedcd_datasets_out = package_cedcd_datasets(df_cedcd_datasets, column_configs)
    if dceg_datasets_exist:
        df_dceg_datasets_out = package_dceg_datasets(df_dceg_datasets, column_configs)
    if nccr_datasets_exist:
        df_nccr_datasets_out = package_nccr_datasets(df_nccr_datasets, column_configs)
    if dctd_datasets_exist:
        df_dctd_datasets_out = package_dctd_datasets(df_dctd_datasets, column_configs)
    if ctd2_datasets_exist:
        df_ctd2_datasets_out = package_ctd2_datasets(df_ctd2_datasets, column_configs,
                                                     ctd2_curated)
    if ctd2_filedata_exist:
        df_ctd2_filedata_out = package_ctd2_filedata(df_ctd2_filedata, column_configs)

    print(f"\n\n Completing post-packaging steps...")

    # Special post-processing handling
    print(f"---\nGenerating enumerated values for data model...")
    if programs_exist:
        get_enum_values(df_programs_out, 
                        cols = config.ENUM_PROGRAM_COLS, 
                        output = config.ENUM_PROGRAM_PATH)

    print(f"---\nScanning packaged TSVs for special character issues...")
    output_files = [(datatype, path) for exists, datatype, path in [
        (programs_exist, 'program', config.PROGRAMS_OUTPUT_PATH),
        (grants_exist, 'grant', config.GRANTS_OUTPUT_PATH),
        (projects_exist, 'project', config.PROJECTS_OUTPUT_PATH),
        (publications_exist, 'publication', config.PUBLICATIONS_OUTPUT_PATH),
        (dbgap_datasets_exist, 'dbgap_dataset',
         config.DBGAP_OUTPUT_CURATED_CLEANED if dbgap_curated
         else config.DBGAP_OUTPUT_PATH),
        (geo_datasets_exist, 'geo_dataset', config.GEO_OUTPUT_PATH),
        (sra_datasets_exist, 'sra_dataset',
         config.SRA_OUTPUT_CURATED_CLEANED if sra_curated
         else config.SRA_OUTPUT_PATH),
        (cedcd_datasets_exist, 'cedcd_dataset', config.CEDCD_OUTPUT_PATH),
        (dceg_datasets_exist, 'dceg_dataset', config.DCEG_OUTPUT_PATH),
        (nccr_datasets_exist, 'nccr_dataset', config.NCCR_OUTPUT_PATH),
        (dctd_datasets_exist, 'dctd_dataset', config.DCTD_OUTPUT_PATH),
        (ctd2_datasets_exist, 'ctd2_dataset',
         config.CTD2_DATASET_CURATED_LOCKED_PATH if ctd2_curated
         else config.CTD2_DATASET_OUTPUT_PATH),
        (ctd2_filedata_exist, 'file', config.CTD2_FILE_OUTPUT_PATH),
    ] if exists]
    run_output_character_qa(output_files, column_configs,
                            config.OUTPUT_CHARACTER_QA_SUMMARY_PATH,
                            config.OUTPUT_CHARACTER_QA_DETAILS_PATH)
        
    # Generate md5.txt
    print(f"---\nGenerating md5 hashes for file validation...")
    generate_md5_hash_file(config.OUTPUT_GATHERED_DIR)

    print(f"\n\n---\nSuccess! INS Data packaging complete.\n---")

# Run module as a standalone script when called directly
if __name__ == "__main__":

    print(f"Running {os.path.basename(__file__)} as standalone module...")

    # Run main packaging function for all data
    package_output_data()