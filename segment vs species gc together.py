import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import re

# Load data
df = pd.read_csv("no_bovine_milk_gc_content.csv")  # Ensure the CSV has columns: Segment, GC_Content, Species

def extract_species_name(species):
    match = re.search(r'Influenza A virus A/(.*?)/Texas', species) or re.search(r'Influenza A virus A/(.*?)/TX', species)
    if match:
        return match.group(1)  # Extracting the part between 'Influenza A virus A/' and '/Texas' or '/TX'
    return species  # If no match is found, return the original string

# Apply the function to create a new 'Species_Short' column
df['Species_Short'] = df['Species'].apply(extract_species_name)

# Convert segment numbers to string for consistent labeling
df["Segment"] = df["Segment"].astype(str)

# Define y-axis limits
y_min = max(40, df["GC_Content"].min())  # Ensuring minimum of 40
y_max = min(50, df["GC_Content"].max())  # Ensuring maximum of 50

# Set up the figure for multiple subplots
fig, axes = plt.subplots(4, 2, figsize=(16, 12))  # 4 rows, 2 columns

# Flatten axes array for easy iteration
axes = axes.flatten()

# Loop through segments 1 to 8 and create a separate plot for each
for i, segment in enumerate(range(1, 9)):
    ax = axes[i]
    seg_df = df[df["Segment"] == str(segment)]
    
    if not seg_df.empty:
        sns.boxplot(data=seg_df, x="Species_Short", y="GC_Content", ax=ax)
        ax.set_title(f"Segment {segment} GC Content vs Species")
        ax.set_xlabel("Species")
        ax.set_ylabel("GC Content")
        ax.tick_params(axis='x', rotation=45)  # Rotate x-axis labels for readability
        ax.set_ylim(y_min, y_max)  # Enforcing y-axis limits

plt.tight_layout()
plt.show()
