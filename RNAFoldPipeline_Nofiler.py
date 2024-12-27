#!/usr/bin/env python3
import sys
import os
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq

# ############################################################
print("Loading the libraries")
# ############################################################

# Ensure required libraries are installed
try:
    import numpy as np
except ImportError:
    print("Please install numpy using pip install numpy")

# ############################################################
print("Reading the input files")
# ############################################################

# Parse command-line arguments
args = sys.argv[1:]
if len(args) < 3:
    print("Usage: script.py <fastadirectory> <inputdir> <outputfilename>")
    sys.exit(1)

fastadirectory = args[0]
inputdir = args[1]
outputfilename = args[2]

# Read FASTA file
fasta_sequences = list(SeqIO.parse(fastadirectory, "fasta"))
fasta_df = pd.DataFrame({
    "seqname": [record.id for record in fasta_sequences],
    "sequence": [str(record.seq) for record in fasta_sequences]
})

# Read all CSV files in the input directory
csv_files = [file for file in os.listdir(inputdir) if file.endswith(".csv")]
csv_data = {file: pd.read_csv(os.path.join(inputdir, file)) for file in csv_files}

# Extract specific CSVs
top = csv_data.get("TopStrandSigPeaks_EnhancedCompare.csv")
comp = csv_data.get("CompStrandSigPeaks_EnhancedCompare.csv")

# Add new columns to top and comp DataFrames
top["Sequence..50nt.upstream"] = ""
top["Utrack"] = ""
top.insert(top.columns.get_loc("PeakCoord") + 1, "strand", "+")
comp["Sequence..50nt.upstream"] = ""
comp["Utrack"] = ""
comp.insert(comp.columns.get_loc("PeakCoord") + 1, "strand", "-")

# ############################################################
print("Calling 50 nt Upstream of every TTS")
# ############################################################

def extract_upstream_sequence(row, sequence_column, is_top=True):
    """Extract upstream sequence based on strand and coordinates."""
    peak_coord = row["PeakCoord"]
    sequence = row[sequence_column]
    if is_top:
        return sequence[max(0, peak_coord - 50):peak_coord]
    else:
        return sequence[peak_coord:peak_coord + 50]

top["Sequence..50nt.upstream"] = top.apply(
    lambda row: extract_upstream_sequence(row, "sequence", is_top=True), axis=1)
comp["Sequence..50nt.upstream"] = comp.apply(
    lambda row: extract_upstream_sequence(row, "sequence", is_top=False), axis=1)

# Reverse complement for complementary strand sequences
comp["Sequence..50nt.upstream"] = comp["Sequence..50nt.upstream"].apply(
    lambda seq: str(Seq(seq).reverse_complement()))

# ############################################################
print("Running RNA fold")
# ############################################################

# Combine top and comp DataFrames for RNAfold input
rnafold_input_df = pd.concat([top, comp]).sort_values(by="PeakCoord")
rnafold_input_df["RNAfolding"] = ""
rnafold_input_df["FreefoldingEnergy"] = ""

# Save RNAfold input sequences to a file
rnafold_input_file = os.path.join(inputdir, "RNAfoldinput.csv")
rnafold_input_df[["Sequence..50nt.upstream"]].to_csv(rnafold_input_file, index=False)

# Run RNAfold (assuming RNAfold is installed and accessible from the command line)
os.system(f"RNAfold < {rnafold_input_file} > RNAfoldoutput.csv")

# Read RNAfold output and parse folding patterns and free energy values
rnafold_output_file = os.path.join(inputdir, "RNAfoldoutput.csv")
rnafold_output_df = pd.read_csv(rnafold_output_file)

# Parse RNA folding patterns and free energy values into the DataFrame
for index, row in rnafold_input_df.iterrows():
    folding_pattern_row = rnafold_output_df.iloc[index]
    rnafold_input_df.at[index, "RNAfolding"] = folding_pattern_row["FoldingPattern"]
    rnafold_input_df.at[index, "FreefoldingEnergy"] = folding_pattern_row["FreeEnergy"]

# ############################################################
print("Calculating the U-track rich regions")
# ############################################################

def calculate_utrack(sequence):
    """Determine if a sequence contains U-track-rich regions."""
    for i in range(len(sequence) - 6):
        if sequence[i:i + 7].count("U") >= 5:
            return "yes"
    return "no"

rnafold_input_df["Utrack"] = rnafold_input_df["FoldingPattern"].apply(calculate_utrack)

# ############################################################
print("Annotating the genome")
# ############################################################

genes_file = [file for file in os.listdir(inputdir) if file.endswith(".genes")][0]
genes_df = pd.read_csv(os.path.join(inputdir, genes_file), sep="\t", header=None)
genes_df.columns = ["genome", "from", "to", "strand", "name", "old_name",
                    "length", "protein_name", "symbol"]

genes_top_strand = genes_df[genes_df["strand"] == "+"]
genes_comp_strand = genes_df[genes_df["strand"] == "-"]

def annotate_genome(peaks_df, genes_strand):
    """Annotate genome regions (TSSUTR, Locus, termUTR) based on peak coordinates."""
    peaks_df["TSSUTR"] = ""
    peaks_df["Locus"] = ""
    peaks_df["termUTR"] = ""

    for i, peak in peaks_df.iterrows():
        peak_coord = peak["PeakCoord"]
        for j, gene in genes_strand.iterrows():
            if gene["to"] < peak_coord <= gene.get("from", float("inf")) + 50:
                peaks_df.at[i, "termUTR"] = gene["old_name"]
            elif gene.get("from", float("-inf")) + 50 <= peak_coord <= gene.get("to", float("inf")):
                peaks_df.at[i, "Locus"] = gene["old_name"]
            elif gene.get("to", float("-inf")) + 300 < peak_coord < gene.get("from", float("inf")) - 300:
                peaks_df.at[i, "TSSUTR"] = gene["old_name"]

annotate_genome(top, genes_top_strand)
annotate_genome(comp, genes_comp_strand)

# ############################################################
print("Analyzing transcriptomic distribution")
# ############################################################

# Define cases based on signal conditions
def classify_cases(df):
    """Classify data into cases based on signal conditions."""
    case1 = df[(df['Sig_Control'] == 1) & (df['Sig_Cond1'] == 1) & (df['Sig_Cond2'] == 1)]
    case2 = df[(df['Sig_Control'] == 1) & (df['Sig_Cond1'] == 0) & (df['Sig_Cond2'] == 0)]
    case3 = df[((df['Sig_Control'] == 0) & (df['Sig_Cond1'] == 0) & (df['Sig_Cond2'] == 1)) |
               ((df['Sig_Control'] == 0) & (df['Sig_Cond1'] == 1) & (df['Sig_Cond2'] == 1))]
    return case1, case2, case3

# Classify top and comp strands into cases
case1_top, case2_top, case3_top = classify_cases(top)
case1_comp, case2_comp, case3_comp = classify_cases(comp)

# Function to calculate transcriptomic distribution
def calculate_distribution(case_df, strand):
    """Calculate the distribution of transcriptomic features for a given case."""
    genic = case_df[case_df["Locus"] != ""]
    TSSUTR = case_df[case_df["TSSUTR"] != ""]
    orphan = case_df[(case_df["Locus"] == "") & (case_df["TSSUTR"] == "") & (case_df["termUTR"] == "")]
    termUTR = case_df[case_df["termUTR"] != ""]
    
    total_number = len(case_df)
    genic_number = len(genic)
    TSSUTR_number = len(TSSUTR)
    orphan_number = len(orphan)
    termUTR_number = len(termUTR)
    
    return {
        "strand": strand,
        "total_number": total_number,
        "genic_number": genic_number,
        "TSSUTR_number": TSSUTR_number,
        "orphan_number": orphan_number,
        "termUTR_number": termUTR_number,
        "genic_percentage": genic_number / total_number * 100 if total_number > 0 else 0,
        "TSSUTR_percentage": TSSUTR_number / total_number * 100 if total_number > 0 else 0,
        "orphan_percentage": orphan_number / total_number * 100 if total_number > 0 else 0,
        "termUTR_percentage": termUTR_number / total_number * 100 if total_number > 0 else 0
    }

# Calculate distributions for each case and strand
distributions = []
for i, (case_top, case_comp) in enumerate(zip([case1_top, case2_top, case3_top], [case1_comp, case2_comp, case3_comp]), start=1):
    distributions.append(calculate_distribution(case_top, f"Case {i} Top"))
    distributions.append(calculate_distribution(case_comp, f"Case {i} Comp"))

# Convert distributions to a DataFrame for easy viewing and export
distribution_df = pd.DataFrame(distributions)

# Save the distribution data to a CSV file
distribution_output_file = os.path.join(inputdir, "Termseq-EnrichedOutput.csv")
distribution_df.to_csv(distribution_output_file, index=False)

print("Transcriptomic distribution analysis complete.")
