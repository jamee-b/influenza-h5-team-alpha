import subprocess
from Bio import SeqIO
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import re  # Import regex for pattern matching

# Define the path to your FASTA file
fasta_file_path = 'BVBRC_genome_sequence_2024_H5N1_Texas.fasta'  

# Initialize an empty list to store sequence information
sequence_data = []

# Read sequences from the FASTA file
for record in SeqIO.parse(fasta_file_path, "fasta"):
    header = record.description  # Header contains strain and segment info
    sequence = str(record.seq)
    sequence_length = len(sequence)
    
    # Extract species and segment from the header
    try:
        species_segment = header.split('|')[1]  # Extract species and segment
        segment = species_segment.split('segment ')[1].split(' ')[0]  # Extract segment number
    except IndexError:
        species, segment = "Unknown", "Unknown"  # Handle unexpected formats
    
    # Append the information to the list
    sequence_data.append([header, segment, sequence_length])  # Store full header initially

# Convert to a pandas DataFrame
df = pd.DataFrame(sequence_data, columns=['species', 'segment', 'sequence_length'])

# Function to extract species name
def extract_species_name(species):
    match = re.search(r'Influenza A virus A/(.*?)/Texas', species) or re.search(r'Influenza A virus A/(.*?)/TX', species)
    return match.group(1) if match else species  # Extract matched group or return original

# Apply function to create a new cleaned 'species' column
df['species'] = df['species'].apply(lambda x: extract_species_name(str(x)))

# Clean up segment formatting
df['segment'] = df['segment'].str.strip()

# Map species to desired groupings
def map_species_to_group(species):
    species_map = {
        'cow': 'cow', 'bovine': 'cow', 'dairy cow': 'cow', 'cattle': 'cow',
        'cat': 'cat', 'feline': 'cat',
        'chicken': 'bird', 'common grackle': 'bird', 'grackle': 'bird', 'blackbird': 'bird'
    }
    # Return the mapped group if exists, else return the original species
    return species_map.get(species.lower(), species)

# Apply the mapping to the 'species' column
df['species'] = df['species'].apply(map_species_to_group)

# Group by 'species' and 'segment' to analyze sequence length variations
grouped = df.groupby(['species', 'segment']).agg({'sequence_length': ['mean', 'std', 'min', 'max']}).reset_index()

# Flatten column names after aggregation
grouped.columns = ['species', 'segment', 'mean_length', 'std_length', 'min_length', 'max_length']

# Plot sequence length variations
plt.figure(figsize=(14, 8))
sns.boxplot(x='species', y='sequence_length', hue='segment', data=df, palette='Set2')

plt.title('Variation in Sequence Lengths Between Species and Segments', fontsize=16)
plt.xlabel('Species', fontsize=12)
plt.ylabel('Sequence Length', fontsize=12)
plt.xticks(rotation=45, ha='right')  # Rotate x-axis labels for better readability
plt.legend(title='Segment', bbox_to_anchor=(1.05, 1), loc='upper left')
plt.tight_layout()

# Show the plot
plt.show()

# Print statistical summary
print(grouped)

