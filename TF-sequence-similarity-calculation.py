import pandas as pd
import parasail
import multiprocessing
import argparse
from tqdm import tqdm
import sys
import math

# Constants for Alignment
# SCORING_MATRIX: BLOSUM62 is a widely used substitution matrix for protein alignment
SCORING_MATRIX = parasail.blosum62
# GAP_OPEN_PENALTY: Penalty for opening a gap in the alignment
GAP_OPEN_PENALTY = 10
# GAP_EXTEND_PENALTY: Penalty for extending an existing gap
GAP_EXTEND_PENALTY = 1

# Function to calculate similarity for a single pair
def calculate_similarity_for_pair(row: dict) -> tuple:
    """
    Calculates the global alignment (Needleman-Wunsch) percent identity
    for a pair of protein sequences using the parasail library.

    Args:
        row (dict): A dictionary representing a row from the input DataFrame,
                    expected to contain 'hTF1', 'Sequence_hTF1', 'hTF2', 'Sequence_hTF2'.

    Returns:
        tuple: A tuple (hTF1, hTF2, similarity_score).
               similarity_score will be a float (percent identity) or math.nan
               if sequences were not found or are invalid.
    """
    hTF1 = row['hTF1']
    hTF2 = row['hTF2']
    seq1 = row['Sequence_hTF1']
    seq2 = row['Sequence_hTF2']
    similarity = row['similarity']

    # Handle cases where sequences were marked as "NOT_FOUND" in the previous step
    if seq1 == "NOT_FOUND" or seq2 == "NOT_FOUND":
        return hTF1, hTF2, math.nan

    # Ensure sequences are valid strings and not empty
    if not isinstance(seq1, str) or not isinstance(seq2, str) or not seq1 or not seq2:
        return hTF1, hTF2, math.nan # Return NaN for invalid or empty sequences

    try:
        # Perform global alignment using Needleman-Wunsch (nw) algorithm from parasail
        # Use parasail.nw_stats to get alignment results with detailed statistics
        alignment_result = parasail.nw_stats(seq1, seq2, GAP_OPEN_PENALTY, GAP_EXTEND_PENALTY, SCORING_MATRIX)

        # Check if alignment_result is not None and its length is valid before calculating percent identity
        if alignment_result and alignment_result.length is not None and alignment_result.length > 0:
            percent_identity = (alignment_result.matches / alignment_result.length) * 100
        else:
            # If alignment fails or results in zero length, consider similarity as 0% or NaN
            percent_identity = 0.0 # Or math.nan if you prefer to distinguish

    except Exception as e:
        # Catch any unexpected errors during the alignment process
        print(f"Error during alignment for pair {hTF1}-{hTF2}: {e}", file=sys.stderr)
        percent_identity = math.nan # Indicate an error for this pair

    return hTF1, hTF2, seq1, seq2, percent_identity, similarity

# Main Script Execution
def main():
    """
    Main function to orchestrate the sequence similarity calculation.
    It reads hTF pairs and their sequences from an input CSV,
    calculates pairwise percent identity using global alignment (Needleman-Wunsch)
    with parallel processing, and saves the results to an output CSV.
    """
    parser = argparse.ArgumentParser(
        description="Calculate pairwise protein sequence similarity (percent identity) "
                    "for human transcription factor pairs using global alignment (Needleman-Wunsch) "
                    "with parallel processing. Expects an input CSV with sequences."
    )
    parser.add_argument(
        "-i", "--input",
        required=True,
        help="Path to the input CSV file containing hTF pairs and their sequences. "
             "Expected columns: 'hTF1', 'Sequence_hTF1', 'hTF2', 'Sequence_hTF2'."
    )
    parser.add_argument(
        "-o", "--output",
        required=True,
        help="Path to the output CSV file where similarity results will be saved. "
             "Columns will be: 'hTF1', 'hTF2', 'Similarity_PercentIdentity'."
    )
    parser.add_argument(
        "-p", "--num_processes",
        type=int,
        default=1, 
        help="Number of parallel processes to use for similarity calculation. "
             "Defaults to (CPU_count - 1) to leave one core free, or 1 if only one CPU is available."
    )
    args = parser.parse_args()

    print(f"Starting sequence similarity calculation script...")
    print(f"Input file: {args.input}")
    print(f"Output file: {args.output}")
    print(f"Using {args.num_processes} parallel processes.")

    try:
        # Read the input CSV file into a pandas DataFrame
        df_input = pd.read_csv(args.input)
        print(f"Successfully loaded {len(df_input)} protein pairs from {args.input}.")
    except FileNotFoundError:
        print(f"Error: Input file '{args.input}' not found. Please check the path.", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"Error reading input file {args.input}: {e}", file=sys.stderr)
        sys.exit(1)

    # Validate that required columns are present in the input DataFrame
    required_columns = ['hTF1', 'Sequence_hTF1', 'hTF2', 'Sequence_hTF2', 'similarity'] 
    if not all(col in df_input.columns for col in required_columns):
        print(f"Error: Input CSV must contain all of the following columns: {required_columns}", file=sys.stderr)
        sys.exit(1)

    # Convert the DataFrame rows into a list of dictionaries. Each dictionary
    # represents a pair and will be passed to `calculate_similarity_for_pair`
    data_for_pool = df_input.to_dict('records')

    print(f"Beginning similarity calculations for {len(data_for_pool)} pairs...")

    # Use multiprocessing.Pool to parallelize the calculations
    with multiprocessing.Pool(processes=args.num_processes) as pool:
        # `imap` is used with `tqdm` to show a dynamic progress bar as results are generated.
        # `chunksize` can be adjusted for performance, but the default often works well
        results = list(tqdm(pool.imap(calculate_similarity_for_pair, data_for_pool),
                            total=len(data_for_pool),
                            desc="Calculating similarities"))

    print("\nAll similarity calculations completed.")

    # Create a new DataFrame from the collected results
    # The columns are defined here to match the output of `calculate_similarity_for_pair`
    df_output = pd.DataFrame(results, columns=['hTF1', 'hTF2', 'Sequence_hTF1', 'Sequence_hTF2', 'Similarity_PercentIdentity', 'Similarity']) 

    # Save the final results to the specified output CSV file
    try:
        df_output.to_csv(args.output, index=False)
        print(f"Successfully saved similarity results to {args.output}.")
    except Exception as e:
        print(f"Error saving output file {args.output}: {e}", file=sys.stderr)
        sys.exit(1)

    print("Script finished.")

if __name__ == "__main__":
    main()