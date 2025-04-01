import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import re

df = pd.read_csv("no_bovine_milk_gc_content.csv")  # Ensure CSV has columns: Segment, GC_Content, Species

# Function to extract species name
def extract_species_name(species):
    match = re.search(r'Influenza A virus A/(.*?)/Texas', species) or re.search(r'Influenza A virus A/(.*?)/TX', species)
    if match:
        return match.group(1)  # Extract species name
    return species  # Return original if no match

# Apply function to "Species" column
df['Species_Short'] = df['Species'].apply(extract_species_name)

# Convert segment numbers to string for consistency
df["Segment"] = df["Segment"].astype(str)

# Debugging: Print min and max GC_Content values
print("GC_Content min:", df["GC_Content"].min())
print("GC_Content max:", df["GC_Content"].max())


# Loop through segments 1 to 8 and create separate plots
for segment in range(1, 9):
    seg_df = df[df["Segment"] == str(segment)]
    
    if not seg_df.empty:
        plt.figure(figsize=(8, 6))
        sns.boxplot(data=seg_df, x="Species_Short", y="GC_Content")

        # Set title and labels
        plt.title(f"Segment {segment} GC Content vs Species")
        plt.xlabel("Species")
        plt.ylabel("GC Content")

        # Rotate x-axis labels for readability
        plt.xticks(rotation=45)

        # Show the plot
        plt.show()

