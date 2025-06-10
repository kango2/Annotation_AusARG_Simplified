import sys
import re

# Function to extract text after the specified pattern
def extract_text(pattern, string):
    match = re.search(pattern, string)
    if match:
        return match.group(1).split()[0]  # Extract text up to the first whitespace
    else:
        return None

# Check if the correct number of arguments is provided
if len(sys.argv) != 2:
    print("Usage: python script.py input_file.tsv")
    sys.exit(1)

input_file = sys.argv[1]

# Read the input TSV file
with open(input_file, 'r') as file:
    lines = file.readlines()

# Initialize variables
current_group = []
groups = []
comment_lines = []

# Process each line in the file
for line in lines:
    # Check if line starts with "#"
    if line.startswith("#"):
        comment_lines.append(line.strip())
        continue
    
    columns = line.strip().split('\t')
    if columns[2] == 'cDNA_match':
        # If column 3 is A, start a new group
        if current_group:
            groups.append(current_group)
            current_group = []
    # Append the line to the current group
    current_group.append(columns)

# Append the last group
if current_group:
    groups.append(current_group)

# Print comment lines first
for comment_line in comment_lines:
    print(comment_line)

# Process each group
for group in groups:
    target_variable = None
    # Find the first instance of column 3 being B and extract the target variable
    for entry in group:
        if entry[2] == 'cDNA_match_part':
            target_variable = extract_text(r'Target=(.*?)(?:;|$)', entry[8])
            break
    # Append the target variable as a new column for each line in the group
    for entry in group:
        if entry[2] == 'cDNA_match':
            print('\t'.join(entry[:9] + [target_variable]))
        else:
            print('\t'.join(entry[:9] + [target_variable]))
