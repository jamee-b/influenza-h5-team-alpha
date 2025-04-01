import os
from itertools import product
from Bio import SeqIO
from Bio.Align import PairwiseAligner

# Debug: Print the current working directory
print("Current Working Directory:", os.getcwd())

# Input and output files
fasta_file = "BVBRC_genome_sequence_2024_H5N1_Texas2.fasta"
output_file = "group_pairwise_alignment2.txt"

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



# Open output file for writing alignments
with open(output_file, "w") as out_f:
    # Perform pairwise alignments and write to file
    for seq1, seq2 in product(bird_HA, mammal_HA):
        alignments = aligner.align(seq1, seq2)
        identity = compute_identity(seq1, seq2)  # Calculate identity percentage
        score = alignments[0].score  # Extract the score
        out_f.write(f"Bird vs Mammal Alignment\nScore={score}\nIdentity={identity:.2f}%\n{alignments[0]}\n\n")

    for seq1, seq3 in product(bird_HA, human_HA):
        alignments = aligner.align(seq1, seq3)
        identity = compute_identity(seq1, seq3)  # Calculate identity percentage
        score = alignments[0].score
        out_f.write(f"Bird vs Human Alignment\nScore={score}\nIdentity={identity:.2f}%\n{alignments[0]}\n\n")

    for seq2, seq3 in product(mammal_HA, human_HA):
        alignments = aligner.align(seq2, seq3)
        identity = compute_identity(seq2, seq3)  # Calculate identity percentage
        score = alignments[0].score
        out_f.write(f"Mammal vs Human Alignment\nScore={score}\nIdentity={identity:.2f}%\n{alignments[0]}\n\n")

print(f"Alignments have been written to {output_file}")

