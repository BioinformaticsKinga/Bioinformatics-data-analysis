import sys
import numpy as np
import matplotlib.pyplot as plt

# Initialize variables
sequences = {}
total_length = 0
total_numer_bp = 0

# Open and read the FASTA file
with open("trinity_izoformy_epip.fasta", "r") as fasta_file:
    current_id = ""
    current_seq = ""
    
    for line in fasta_file:
        line = line.strip()
        
        # If the line is a header (starts with '>'), process the previous sequence
        if line.startswith(">"):
            if current_id:
                sequences[current_id] = len(current_seq)
                total_length += len(current_seq)
                total_numer_bp += 1
            
            # Store new sequence ID and reset sequence
            current_id = line[1:]
            current_seq = ""
        else:
            current_seq += line

    # Add the last sequence to the dictionary
    if current_id:
        sequences[current_id] = len(current_seq)
        total_length += len(current_seq)
        total_numer_bp += 1

# Perform calculations
sequence_lengths = list(sequences.values())
sequence_lengths.sort(reverse=True)

total_bp = total_length
minimum_length = min(sequence_lengths)
maximum_length = max(sequence_lengths)
average_length = total_length / total_numer_bp

# Calculate N50
cumulative_length = 0
n50_length = None
for length in sequence_lengths:
    cumulative_length += length
    if cumulative_length >= total_length / 2:
        n50_length = length
        break

# Prepare data for the table
table_data = [
    ["Total Number", total_numer_bp],
    ["Total Length (bp)", total_bp],
    ["N50 Length (bp)", n50_length],
    ["Minimum Length (bp)", minimum_length],
    ["Maximum Length (bp)", maximum_length],
    ["Average Length (bp)", round(average_length, 2)]
]

# Create a table using Matplotlib
fig, ax = plt.subplots(figsize=(8, 6))
ax.axis('off')
table = ax.table(cellText=table_data, cellLoc='center', loc='center')
table.auto_set_font_size(False)
table.set_fontsize(12)
table.scale(1, 1.5)
plt.tight_layout()

# Save the table as an image
plt.savefig("dlugosci.png", bbox_inches='tight')
plt.show()
