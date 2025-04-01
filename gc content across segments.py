import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import re

# Load the data from CSV file
df = pd.read_csv('gc_analysis_tx.csv')

def extract_segment(segment):
    match = re.search("Segment", segment)
    if match:
        return match.group(1)  
    return segment  # If no match is found, return the original string


# Create the boxplot using the new 'Species_Short' column
plt.figure(figsize=(10, 6))
sns.boxplot(x='Segment', y='GC_Content', data=df)
plt.xticks(rotation=90)  # Rotate the species labels for better readability
plt.title('GC Content Across Different Segments')
plt.xlabel('Segment')
plt.ylabel('GC Content (%)')
plt.tight_layout()
plt.show()
