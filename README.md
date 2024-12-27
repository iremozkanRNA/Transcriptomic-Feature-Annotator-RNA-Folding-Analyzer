# Transcriptomic-Feature-Annotator-RNA-Folding-Analyzer
This script processes genomic and transcriptomic data, extracts upstream sequences, predicts RNA folding patterns, identifies U-track-rich regions, annotates genomic features (e.g., UTRs, loci), and analyzes transcriptomic distributions across different cases. The pipeline is designed to work with FASTA files, CSV files containing peak coordinates, and gene annotation files.

# Features
1. Library Loading: Automatically loads required libraries like `pandas` and `Biopython`.
2. 	Input File Handling:
    Reads a FASTA file containing genomic sequences.
    Processes CSV files with peak coordinates.
   	Reads gene annotation files for genome feature annotation.
3. Upstream Sequence Extraction:
    Extracts 50 nucleotides upstream (or downstream for complementary strands) of transcription termination sites (TTS).
4. RNA Folding Prediction:
    Generates RNA folding patterns and calculates free energy using the external `RNAfold` tool.
5. U-Track Detection:
    Identifies U-track-rich regions in RNA folding patterns.
6. Genome Annotation:
    Annotates peaks with genomic features such as 5’UTR, 3’UTR, genic regions, or orphan regions.
7. Transcriptomic Distribution Analysis:
    Categorizes peaks into cases based on signal conditions and calculates their distribution across genomic features.
# PREREQUISITES
# Software Requirements
R: Version 4.x or higher.
RStudio: Recommended for an integrated development environment.
# Required R Packages
Biostrings: For handling biological sequences.
tidyverse: For data manipulation and visualization.
stringr: For string manipulation.
dplyr: For advanced data manipulation.
# External Tools
RNAfold (from the ViennaRNA package): Required for RNA secondary structure prediction.
# Input Files
FASTA File: Contains genomic sequences (e.g., `genome.fasta`).
CSV Files: Contain peak coordinates (e.g., `TopStrandSigPeaks_EnhancedCompare.csv`).
Gene Annotation File: Tab-delimited file with gene information (e.g., `genome.genes`).

# Usage
# Running the Script
Open the script in RStudio or your preferred R environment.
Ensure all required packages are installed using the following commands:
> if (!requireNamespace("BiocManager", quietly = TRUE))
	install.packages("BiocManager")
> BiocManager::install("Biostrings")
> install.packages("tidyverse")
> install.packages("stringr")
> install.packages("dplyr")

Run the script with the appropriate command-line arguments:
> Rscript script.R <fastadirectory> <inputdir> <outputfilename>
Example
> Rscript script.R genome.fasta input_data/ output_results.csv

# Outputs
1. RNAfold Input File (`RNAfoldinput.csv`):
    Contains sequences for RNA folding prediction.
2. RNAfold Output File (`RNAfoldoutput.csv`):
    Contains RNA folding patterns and free energy values.
3. Final Output File (`output_results.csv`):
    Annotated peaks with genomic features, RNA folding data, U-track analysis, and case classifications.
4. Enriched Output File (`Termseq-EnrichedOutput.csv`):
    Summary of transcriptomic distribution across cases.

# Workflow Breakdown
Step 1: Load Libraries and Read Input Files
	•	Loads necessary libraries and reads input FASTA and CSV files.
Step 2: Extract Upstream Sequences
	•	Extracts sequences upstream of TTS based on peak coordinates.
Step 3: Predict RNA Folding Patterns
	•	Uses `RNAfold` to predict RNA structures and calculate free energy values.
Step 4: Identify U-Track-Rich Regions
	•	Analyzes RNA folding patterns for U-rich regions.
Step 5: Annotate Genomic Features
	•	Annotates peaks with genomic features using gene annotation data.
Step 6: Analyze Transcriptomic Distribution
	•	Categorizes peaks into cases based on signal conditions and analyzes their distribution across genomic features.

# Notes
Ensure that `RNAfold` is installed and accessible from the command line.
Modify column names in the script if your input files have different headers than expected.
The script is optimized for typical datasets but may require adjustments for very large datasets.

# Troubleshooting
1. Missing Libraries:
	•	Install missing R packages using `install.packages()` or `BiocManager::install()` as needed.
2. RNAfold Not Found:
	•	Verify that the ViennaRNA package is installed and added to your system’s PATH.
3. Incorrect Input Format:
	•	Ensure input files are correctly formatted according to expectations (FASTA for sequences, CSV for peaks).
4. Empty Outputs:
	•	Check input files for valid data matching expected formats.

# Author & License
This script was written in R for transcriptomic data analysis and genomic feature annotation.

Author: Irem Ozkan

Date: December 27, 2024
