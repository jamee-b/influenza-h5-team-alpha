import re

# Function to read alignment file and extract scores
def extract_alignment_scores(file_path):
    scores = []
    score_pattern = re.compile(r"Score=([\d\.\-]+)")  # Regex to match scores

    with open(file_path, "r") as file:
        for line in file:
            match = score_pattern.search(line)
            if match:
                scores.append(float(match.group(1)))  # Extract and convert score to float

    return scores

alignment_file = "group_pairwise_alignment.txt"
scores = extract_alignment_scores(alignment_file)

# Ensure there are scores before computing statistics
if scores:
    print(f"Total Alignments: {len(scores)}")
    print(f"Highest Score: {max(scores)}")
    print(f"Lowest Score: {min(scores)}")
    print(f"Average Score: {sum(scores)/len(scores):.2f}")
else:
    print("No alignment scores found.")
