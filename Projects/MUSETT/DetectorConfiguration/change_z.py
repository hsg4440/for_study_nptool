import numpy as np


dz = np.array([0.25,0.61,0.36,0.47])*0.760

cnt = 0
def process_line(line, dz_value):
    print(line, dz_value)
    # Check if the line contains the '=' character
    if '=' in line:
        # Split the line into key and values
        key, values = line.split('=')
        # Split the values and convert to float
        nums = values.split()
        # Adjust the last number (before 'mm') by the given dz_value
        nums[-2] = f"{float(nums[-2]) + dz_value:.2f}"
        # Reconstruct the line
        new_line = f" {key}= {nums[0]} {nums[1]} {nums[2]} mm"
        return new_line
    else:
        # If the line does not contain '=', return it as is
        return line

# Read the file
with open('MUSETT_NPS.detector', 'r') as file:
    lines = file.readlines()

dz_index =-1
# Process each line
processed_lines = []
for line in lines:
    print(dz_index)
    stripped_line = line.strip()
    if stripped_line == "MUSETT":
        dz_index += 1
    if dz_index >= 0 and dz_index < len(dz):
        processed_lines.append(process_line(stripped_line, dz[dz_index]))
    else:
        processed_lines.append(stripped_line)

# Write the updated lines back to a file
with open('updated_data.txt', 'w') as file:
    for line in processed_lines:
        file.write(line + '\n')
