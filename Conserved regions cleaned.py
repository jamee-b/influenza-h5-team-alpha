import subprocess
from Bio import SeqIO, AlignIO
import os
import numpy as np
import matplotlib.pyplot as plt

# Define file paths
fasta_file_path = "BVBRC_genome_sequence_2024_H5N1_Texas.fasta"
output_dir = "alignment_results"
os.makedirs(output_dir, exist_ok=True)

# Read and group sequences by segment
segment_sequences = {}
for record in SeqIO.parse(fasta_file_path, "fasta"):
    header = record.description
    try:
        segment_number = header.split("segment ")[1].split(" ")[0].strip()
    except IndexError:
        continue  # Skip if segment number isn't found
    segment_sequences.setdefault(segment_number, []).append(record)

# Function to run Clustal Omega via subprocess
def run_clustalo(input_fasta, output_fasta):
    """Runs Clustal Omega for sequence alignment."""
    clustal_path = r"C:\Users\brandlja\Documents\Clustal Omega\clustal-omega-1.2.2-win64\clustalo.exe"  # Update this with the actual path
    try:
        subprocess.run([clustal_path, "-i", input_fasta, "-o", output_fasta, "--force", "-v"], check=True)
    except subprocess.CalledProcessError as e:
        print(f"Error running Clustal Omega: {e}")

# Perform MSA and store results
aligned_files = []
for segment, records in segment_sequences.items():
    # Only process segments 1, 4, and 6
    if segment in ['1', '4', '6']:
        segment_fasta = os.path.join(output_dir, f"segment_{segment}.fasta")
        aligned_file = os.path.join(output_dir, f"segment_{segment}_aligned.fasta")

        SeqIO.write(records, segment_fasta, "fasta")
        run_clustalo(segment_fasta, aligned_file)
        aligned_files.append(aligned_file)

# Function to calculate conservation scores
def calculate_conservation(alignment_file):
    alignment = AlignIO.read(alignment_file, "fasta")
    seq_count = len(alignment)
    seq_length = alignment.get_alignment_length()
    
    conservation_scores = []
    for i in range(seq_length):
        column = [seq[i] for seq in alignment]
        unique_bases = set(column)
        conservation_scores.append(1 - (len(unique_bases) / seq_count))  # Higher means conserved
    
    return np.array(conservation_scores)

# Plot Conservation Scores for segments 1, 4, and 6
plt.figure(figsize=(12, 6))
for aligned_file in aligned_files:
    segment = aligned_file.split("_")[-2]
    conservation_scores = calculate_conservation(aligned_file)
    plt.plot(conservation_scores, label=f"Segment {segment}")

plt.title("Conserved and Variable Regions Across Segments 1, 4, and 6")
plt.xlabel("Position in Sequence")
plt.ylabel("Conservation Score (1 = Highly Conserved)")
plt.legend()
plt.show()
