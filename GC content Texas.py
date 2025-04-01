from Bio import SeqIO
import pandas as pd
import os
import re  # Import regular expressions for extracting species

def calculate_gc_content(sequence):
    """Calculate GC content percentage."""
    gc_count = sequence.count("g") + sequence.count("c")
    return (gc_count / len(sequence)) * 100 if len(sequence) > 0 else 0

def extract_species_and_segment(record):
    """Extract species name and segment from the FASTA header."""
    species_match = re.search(r"\[(.*?)\]", record.description)  # Find text inside brackets for species
    segment_match = re.search(r"segment (\d+)", record.description)  # Extract segment number
    
    species = species_match.group(1) if species_match else "Unknown"
    segment = segment_match.group(1) if segment_match else "Unknown"  # Default to "Unknown" if no segment is found
    return species, segment

def analyze_fasta(file_path, file_name):
    """Extract GC content, sequence length, and species for each record in a FASTA file."""
    results = []
    for record in SeqIO.parse(file_path, "fasta"):
        gc_content = calculate_gc_content(str(record.seq))
        species, segment = extract_species_and_segment(record)  # Extract species name from the header
        results.append([record.id, len(record.seq), gc_content, species, segment, file_name])
    return results

# File path
directory = "C:/Users/brandlja/Documents/Capstone/"
file_name = "BVBRC_genome_sequence_2024_H5N1_Texas.fasta"
file_path = os.path.join(directory, file_name)

# Analyze the FASTA file
all_results = analyze_fasta(file_path, file_name)

# Convert results to DataFrame
final_df = pd.DataFrame(all_results, columns=["Sequence_ID", "Length", "GC_Content", "Species", "Segment", "File"])

# Save to CSV
output_file = "gc_analysis_tx.csv"
final_df.to_csv(output_file, index=False)
print(f"Analysis complete. Results saved to {output_file}.")
print(final_df.head())




