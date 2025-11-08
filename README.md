# hTF-Analysis: Protein Sequence Similarity Calculations

This repository provides a comprehensive pipeline for analyzing the sequence similarity of human Transcription Factor (hTF) protein pairs. It automates the retrieval of protein sequences from UniProt, calculates pairwise global sequence identity using the ultra-fast `parasail` library, and offers a Jupyter Notebook for visualizing the distributions and relationships of these similarity scores with other metrics (e.g., Jaccard similarity).

The pipeline is optimized for execution on the Boston University Shared Computing Cluster (BU SCC) using the SGE job scheduler, enabling efficient parallel processing for computationally intensive steps on large datasets.


## Workflow

The pipeline consists of the following sequential stages:

1.  **Initial Data Preparation & Sequence Retrieval**
    *   **Description:** This stage fetches canonical protein sequences for a given list of human Transcription Factor (hTF) pairs from the UniProt database. The Python script handles UniProt API rate limits to ensure polite usage and identifies hTFs for which sequences cannot be found, marking them appropriately.
    *   **Scripts:** `TF-sequence-query.py` (Python script), `run-sequence-api-requests.sh` (SGE job submission script).
    *   **Input:** A CSV file (e.g., `data/pairs_htf_htf.csv`) listing hTF pairs with columns like `hTF1` and `hTF2`.
    *   **Output:** A CSV file (e.g., `data/pairs_htf_htf_with_sequences.csv`) containing the original hTF IDs along with their fetched protein sequences.

2.  **Protein Sequence Similarity Calculation**
    *   **Description:** After sequence retrieval, this stage calculates the global alignment (Needleman-Wunsch) percent identity for all protein pairs. It leverages the `parasail` library for ultra-fast alignment and utilizes multiprocessing for efficient parallel execution across multiple CPU cores on the cluster. BLOSUM62 is used as the scoring matrix.
    *   **Scripts:** `TF-sequence-similarity-calculation.py` (Python script), `run-sequence-similarity-calculation.sh` (SGE job submission script).
    *   **Input:** The CSV file generated in the previous step (e.g., `data/pairs_htf_htf_with_sequences.csv`).
    *   **Output:** A new CSV file (e.g., `data/pairs_htf_htf_with_sequences_with_sequence-similarity.csv`) which includes an additional column `Similarity_PercentIdentity` for each pair.

3.  **Analysis and Visualization of Results**
    *   **Description:** This final stage involves processing and visualizing the calculated sequence similarity results. The provided Jupyter Notebook can be used to generate various plots (e.g., scatter plots, histograms) and summary statistics, allowing for insights into the distribution of sequence similarity and its relationship with other metrics, such as Jaccard similarity.
    *   **Script:** `analyse-hTF-pair-similarity.ipynb` (Jupyter Notebook).
    *   **Input:** The CSV file generated from the similarity calculation step (e.g., `data/pairs_htf_htf_with_sequences_with_sequence-similarity_adjusted-headers.csv`). *Note: Ensure the column names in your input file match what the notebook expects (e.g., `Sequence_PercentIdentity`, `Jaccard_Similarity`).*
    *   **Output:** Various plots (saved figures or displayed in the notebook) and printed summary tables for detailed analysis.

## Directory Structure

The expected directory structure and key files within this repository are as follows:

.
├── data/                                   # Directory for all input and intermediate/final data files
│   ├── pairs_htf_htf.csv                   # Initial list of hTF pairs for sequence retrieval
│   ├── pairs_htf_htf_with_sequences.csv    # Output from the sequence retrieval step
│   └── pairs_htf_htf_with_sequences_with_sequence-similarity.csv # Final output from similarity calculation
├── TF-sequence-query.py                    # Python script for querying UniProt for protein sequences
├── TF-sequence-similarity-calculation.py   # Python script for calculating pairwise sequence similarity
├── analyse-hTF-pair-similarity.ipynb       # Jupyter Notebook for results analysis and visualization
├── parasail-python_env.yml                 # Conda environment definition for parasail calculations
├── run-sequence-api-requests.sh            # SGE job submission script for sequence retrieval
└── run-sequence-similarity-calculation.sh  # SGE job submission script for similarity calculation


## Setup and Prerequisites

### Conda Environments

You will need to set up at least two Conda environments to run the entire pipeline effectively:

1.  **`jupyter_env` (for sequence retrieval and general data handling):**
    This environment is used by `TF-sequence-query.py` and for running the `analyse-hTF-pair-similarity.ipynb` notebook. It should contain `python`, `pandas`, `requests`, `tqdm`, `seaborn`, and `matplotlib`. If you have an existing `jupyter_env.yml` (as provided in previous requests), you can use that. Otherwise, a minimal setup would look like this:

    ```bash
    # Create the environment (adjust Python version if needed)
    conda create -n jupyter_env python=3.12 pandas requests tqdm jupyterlab ipykernel seaborn matplotlib -y

    # Activate it
    conda activate jupyter_env

    # Install ipykernel for Jupyter Notebook integration
    python -m ipykernel install --user --name=jupyter_env --display-name "Python (jupyter_env)"
    ```

2.  **`parasail_env` (for sequence similarity calculation):**
    This environment is specifically defined by the `parasail-python_env.yml` file and is used by `TF-sequence-similarity-calculation.py`.

    ```bash
    conda env create -f parasail-python_env.yml
    conda activate parasail_env
    ```

### BU SCC Specifics

*   The provided `.sh` scripts (`run-sequence-api-requests.sh`, `run-sequence-similarity-calculation.sh`) are designed for submission to the Boston University Shared Computing Cluster (BU SCC) using the SGE job scheduler.
*   Remember to update the `#$ -P myproject` directive (e.g., `#$ -P cancergrp`) in your `.sh` scripts to your specific SCC project.
*   Resource requests such as hard time limits (`-l h_rt`), memory per core (`-l mem_per_core`), and the number of cores for parallel execution (`-pe omp`) should be adjusted based on your dataset size and the cluster's available resources.
*   The `TF-sequence-similarity-calculation.py` script automatically utilizes the `$NSLOTS` environment variable (provided by SGE for parallel environments) to configure the number of parallel processes, ensuring optimal resource usage.

---

## Usage

Follow these steps to execute the hTF sequence similarity analysis pipeline on the BU SCC:

### 1. Initial Data Preparation & Sequence Retrieval

Ensure your initial input file (`data/pairs_htf_htf.csv`) is present in the `data/` directory relative to where you submit the job.

Submit the sequence retrieval script to the SCC:

```bash
qsub run-sequence-api-requests.sh

    Input File Example: data/pairs_htf_htf.csv (expected to contain hTF1 and hTF2 columns).
    Output File Example: data/pairs_htf_htf_with_sequences.csv.
```

### 2. Protein Sequence Similarity Calculation

Once the sequence retrieval job completes and data/pairs_htf_htf_with_sequences.csv is generated, proceed to calculate the pairwise percent identity.

Important: The run-sequence-similarity-calculation.sh script currently uses hardcoded absolute paths for its input and output files (/projectnb/cancergrp/Philipp/data/...). Ensure that the output from Step 1 is located at the specified input path, or adjust the paths within the run-sequence-similarity-calculation.sh script to match your setup.

Submit the similarity calculation script to the SCC:

```bash
qsub run-sequence-similarity-calculation.sh

    Input File Example: /projectnb/cancergrp/Philipp/data/pairs_htf_htf_with_sequences.csv (this path must contain the output from Step 1).
    Output File Example: /projectnb/cancergrp/Philipp/data/pairs_htf_htf_with_sequences_with_sequence-similarity.csv.
```

### 3. Analysis and Visualization of Results

After the similarity calculations are complete, you can analyze and visualize the results using the provided Jupyter Notebook. This step is typically run interactively.

  Activate your jupyter_env:

   ```bash
    conda activate jupyter_env
   ```
    
Open and run all cells in the notebook: You can do this by launching Jupyter Lab or Jupyter Notebook and navigating to the file.

```bash
    jupyter lab analyse-hTF-pair-similarity.ipynb
    # Or, on SCC, you might use an interactive job or a dedicated Jupyter server setup

    Input File Example: data/pairs_htf_htf_with_sequences_with_sequence-similarity_adjusted-headers.csv (ensure your file's headers, especially for Sequence_PercentIdentity and Jaccard_Similarity, match what the notebook expects).
    Output: The notebook will display various plots and printed summary statistics, providing insights into the sequence similarity data.
```

## Authorship

This pipeline and its associated scripts were solely developed by Philipp Trollmann during his second PhD rotation in Dr. Juan Fuxman Bass's lab at Boston University.

---
For further details on specific script functionalities, please refer to the extensive annotations within each script file, or feel free to reach out.

If you use or adapt this pipeline for your research, please cite this repository and attribute appropriately.
