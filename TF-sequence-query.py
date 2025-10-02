import pandas as pd
import requests
import time
import sys
from tqdm import tqdm

# --- Configuration ---
# Define input and output file paths
INPUT_FILE = "data/pairs_htf_htf.csv"
OUTPUT_FILE = "data/pairs_htf_htf_with_sequences.csv" 

# Define column names for the input CSV file
INPUT_HTF1_COLUMN = "hTF1" 
INPUT_HTF2_COLUMN = "hTF2" 

# Define column names for the output CSV file
OUTPUT_HTF1_COLUMN = "hTF1"
OUTPUT_SEQUENCE_HTF1_COLUMN = "Sequence_hTF1"
OUTPUT_HTF2_COLUMN = "hTF2"
OUTPUT_SEQUENCE_HTF2_COLUMN = "Sequence_hTF2"
OUTPUT_SIMILARITY_COLUMN = "similarity" # To include the similarity column from input

# UniProt API base URL for sequence search
UNIPROT_API_SEARCH_URL = "https://rest.uniprot.org/uniprotkb/search"

# Rate limit for UniProt API requests (in seconds)
# This delay prevents hitting UniProt's rate limits and ensures polite API usage.
# UniProt recommends max 1 request/second for batch, or 0.1-0.5 for individual.
API_REQUEST_DELAY = 0.2  # 200 milliseconds

# --- Function to fetch sequence from UniProt by protein name ---
def fetch_sequence_from_protein_name(protein_name: str) -> str | None:
    """
    Fetches the protein sequence for a given human transcription factor (hTF) name from UniProtKB.
    It queries UniProt for reviewed human entries matching the protein name and extracts
    the primary protein sequence.

    Args:
        protein_name (str): The name of the human transcription factor (e.g., HOXB8).

    Returns:
        str | None: The protein sequence if found, otherwise None.
    """
    # Construct the UniProt query: protein_name AND human (taxonomy_id:9606) AND reviewed
    query_params = {
        "query": f"{protein_name} AND (taxonomy_id:9606) AND (reviewed:true)",
        "fields": "sequence",  # Request only the sequence field
        "format": "json"
    }

    try:
        response = requests.get(UNIPROT_API_SEARCH_URL, params=query_params)
        response.raise_for_status()  # Raise an HTTPError for bad responses (4xx or 5xx)
        data = response.json()

        if data and 'results' in data and len(data['results']) > 0:
            # Return the sequence from the first result found (typically the canonical one)
            sequence = data['results'][0]['sequence']['value']
            return sequence
        else:
            # No reviewed UniProt entry found for the given protein name
            # print(f"  No reviewed human sequence found for '{protein_name}'.", file=sys.stderr)
            return None

    except requests.exceptions.RequestException as e:
        print(f"\nError fetching sequence for protein '{protein_name}': {e}", file=sys.stderr)
        return None
    except KeyError as e:
        print(f"\nError parsing UniProt response for protein '{protein_name}'. Missing key: {e}", file=sys.stderr)
        return None
    except IndexError:
        print(f"\nNo sequence value found in UniProt response for protein '{protein_name}'.", file=sys.stderr)
        return None

# --- Main Script Execution ---
def main():
    """
    Main function to orchestrate the script execution.
    It reads interacting hTF pairs, fetches their sequences from UniProt,
    and saves the combined data to a new CSV file.
    """
    print(f"Starting to process {INPUT_FILE}...")

    try:
        # Read the input CSV file. It has a header and uses commas as separators.
        df_input = pd.read_csv(INPUT_FILE, sep=',')
        print(f"Read {len(df_input)} protein pairs from {INPUT_FILE}.")
    except FileNotFoundError:
        print(f"Error: Input file '{INPUT_FILE}' not found. Please check the path.", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"Error reading input file {INPUT_FILE}: {e}", file=sys.stderr)
        sys.exit(1)

    # 1. Collect all unique hTF names from both columns
    all_htf_names = set()
    all_htf_names.update(df_input[INPUT_HTF1_COLUMN].unique())
    all_htf_names.update(df_input[INPUT_HTF2_COLUMN].unique())
    print(f"Found {len(all_htf_names)} unique hTF names to query for sequences.")

    # 2. Fetch sequences for all unique names and store them in a dictionary (acting as a cache)
    sequence_map = {}
    # tqdm provides a progress bar for the API calls, enhancing user experience.
    for htf_name in tqdm(list(all_htf_names), desc="Fetching sequences from UniProt"):
        if htf_name not in sequence_map:
            sequence = fetch_sequence_from_protein_name(htf_name)
            sequence_map[htf_name] = sequence
            time.sleep(API_REQUEST_DELAY)  # Pause to respect UniProt API rate limits

    print("\nFinished fetching all unique sequences.")

    # 3. Create the output DataFrame by mapping fetched sequences to the original hTFs
    # Initialize with original columns
    df_output = df_input.copy()

    # Map hTF names to their fetched sequences. Names for which no sequence was found
    # (i.e., `sequence_map` contains None) will be filled with "NOT_FOUND".
    df_output[OUTPUT_SEQUENCE_HTF1_COLUMN] = df_output[OUTPUT_HTF1_COLUMN].map(sequence_map).fillna("NOT_FOUND")
    df_output[OUTPUT_SEQUENCE_HTF2_COLUMN] = df_output[OUTPUT_HTF2_COLUMN].map(sequence_map).fillna("NOT_FOUND")

    # Reorder columns for clarity (optional, but good practice)
    df_output = df_output[[
        OUTPUT_HTF1_COLUMN, OUTPUT_SEQUENCE_HTF1_COLUMN,
        OUTPUT_HTF2_COLUMN, OUTPUT_SEQUENCE_HTF2_COLUMN,
        OUTPUT_SIMILARITY_COLUMN # Include the similarity column
    ]]


    # 4. Save the resulting DataFrame to a CSV file
    try:
        df_output.to_csv(OUTPUT_FILE, index=False)
        print(f"Successfully saved full dataset with sequences to {OUTPUT_FILE}.")
    except Exception as e:
        print(f"Error saving output file {OUTPUT_FILE}: {e}", file=sys.stderr)
        sys.exit(1)

    print("Script finished.")

if __name__ == "__main__":
    main()