from itertools import product
from Bio import SeqIO
from Bio.Align import PairwiseAligner

fasta_file = "BVBRC_genome_sequence_2024_H5N1_Texas.fasta"
output_file = "chicken_dairycow_pairwise_alignment.txt"

# Parse FASTA file and store sequences in a dictionary
sequences = {record.description: str(record.seq) for record in SeqIO.parse(fasta_file, "fasta")}

# Extract HA sequences
chicken_HA = [seq for id, seq in sequences.items() if "chicken" in id.lower() and "(HA)" in id]
dairy_cow_HA = [seq for id, seq in sequences.items() if any(x in id.lower() for x in ["dairy cow"]) and "(HA)" in id]

# Debug: Check if we have sequences
print(f"Found {len(chicken_HA)} Chicken HA sequences.")
print(f"Found {len(dairy_cow_HA)} Dairy Cow HA sequences.")

if not chicken_HA or not dairy_cow_HA:
    print("No matching sequences found. Check your FASTA headers.")
else:
    # Initialize aligner with adjusted scoring
    aligner = PairwiseAligner()
    aligner.mode = "global"
    aligner.match_score = 2  # Reward for matching bases
    aligner.mismatch_score = -1  # Penalty for mismatches
    aligner.open_gap_score = -3  # Open gap penalty
    aligner.extend_gap_score = -1  # Gap extension penalty

    # Compare each chicken HA with each dairy cow HA
    for seq1, seq2 in product(chicken_HA, dairy_cow_HA):
        alignments = aligner.align(seq1, seq2)
        
        # Check if there are alignments and print only the first one
        if alignments:
            print(alignments[0])  # Print the first alignment only
        else:
            print("No alignment found for this pair.")

    print(f"Alignments have been written to {output_file}")
