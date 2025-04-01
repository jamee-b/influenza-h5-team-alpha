import os
import matplotlib.pyplot as plt
import numpy as np
from itertools import product
from Bio import SeqIO
from Bio.Align import PairwiseAligner

# Debug: Print the current working directory
print("Current Working Directory:", os.getcwd())

# Input and output files
fasta_file = "BVBRC_genome_sequence_2024_H5N1_Texas2.fasta"
output_file = "group_pairwise_alignment4.txt"

# Parse FASTA file and store sequences in a dictionary
sequences = {record.description: str(record.seq) for record in SeqIO.parse(fasta_file, "fasta")}

# Extract HA sequences based on organism categories
bird_HA = [seq for id, seq in sequences.items() if any(x in id.lower() for x in ["chicken", "common grackle", "blackbird", "grackle"]) and "(HA)" in id]
mammal_HA = [seq for id, seq in sequences.items() if any(x in id.lower() for x in ["dairy cow", "cat", "feline", "cattle", "bovine"]) and "(HA)" in id]
human_HA = [seq for id, seq in sequences.items() if "viet nam" in id.lower() and "(HA)" in id]

# Debugging: Check if sequences are found
print(f"Found {len(bird_HA)} Bird HA sequences.")
print(f"Found {len(mammal_HA)} Mammal HA sequences.")
print(f"Found {len(human_HA)} Human HA sequences.")

# Exit if any group is empty
if not bird_HA or not mammal_HA or not human_HA:
    print("No matching sequences found. Check your FASTA headers.")
    print("Available Headers:", list(sequences.keys())[:10])  # Show some headers for debugging
    exit()

# Initialize aligner with adjusted scoring
aligner = PairwiseAligner()
aligner.mode = "global"
aligner.match_score = 2
aligner.mismatch_score = -1
aligner.open_gap_score = -3
aligner.extend_gap_score = -1

# Function to compute identity percentage
def compute_identity(seq1, seq2):
    alignments = aligner.align(seq1, seq2)
    best_alignment = alignments[0]  # We are using the first (best) alignment

    # Extract aligned sequences (aligned indices)
    aligned_seq1, aligned_seq2 = best_alignment[0], best_alignment[1]

    # Remove gaps (represented by '-')
    aligned_seq1_str = [base for base in aligned_seq1 if base != '-']
    aligned_seq2_str = [base for base in aligned_seq2 if base != '-']

    # Calculate the number of matches between the two aligned sequences
    matches = sum(1 for a, b in zip(aligned_seq1_str, aligned_seq2_str) if a == b)
    total = len(aligned_seq1_str)  # This is the total number of bases excluding gaps

    # Calculate identity percentage
    identity = (matches / total) * 100
    return identity

# Function to generate a dot plot with custom labels for axes
def dot_plot(seq1, seq2, ax, label1, label2):
    matrix = np.zeros((len(seq1), len(seq2)))

    for i in range(len(seq1)):
        for j in range(len(seq2)):
            if seq1[i] == seq2[j]:
                matrix[i, j] = 1  # Mark matching positions

    ax.imshow(matrix, cmap="Greys", aspect="auto")
    ax.set_title(f"Dot Plot: {label1} vs {label2}")
    ax.set_xlabel(label2)
    ax.set_ylabel(label1)

    # Set x and y ticks to correspond to sequence names
    ax.set_xticks([i for i in range(len(seq2))])
    ax.set_yticks([i for i in range(len(seq1))])
    
    # Label the ticks with "Bird", "Mammal", "Human"
    ax.set_xticklabels([label2] * len(seq2))  # All x-axis ticks labeled as "Mammal" or "Human"
    ax.set_yticklabels([label1] * len(seq1))  # All y-axis ticks labeled as "Bird" or "Human"

# Open output file for writing alignments and analysis
with open(output_file, "w") as out_f:
    # Create the figure with subplots for three comparisons
    fig, axes = plt.subplots(1, 3, figsize=(18, 6))
    
    # Perform pairwise alignments and write to file
    for idx, (seq1, seq2, comparison_name, label1, label2) in enumerate(
        [(bird_HA[0], mammal_HA[0], "Bird vs Mammal", "Bird", "Mammal"),
         (bird_HA[0], human_HA[0], "Bird vs Human", "Bird", "Human"),
         (mammal_HA[0], human_HA[0], "Mammal vs Human", "Mammal", "Human")]):
        
        alignments = aligner.align(seq1, seq2)
        score = alignments[0].score  # Extract the score
        identity = compute_identity(seq1, seq2)  # Calculate identity percentage
        out_f.write(f"{comparison_name} Alignment\nScore={score}\nIdentity={identity:.2f}%\n{alignments[0]}\n\n")
        
        # Generate dot plot for this alignment in the appropriate subplot
        dot_plot(seq1, seq2, axes[idx], label1, label2)

    # Show the plot
    plt.tight_layout()
    plt.show()

print(f"Alignments and analyses have been written to {output_file}")


