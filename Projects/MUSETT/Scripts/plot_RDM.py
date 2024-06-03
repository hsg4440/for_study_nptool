import matplotlib.pyplot as plt
import numpy as np

plt.rcParams['font.size'] = 16

# Read data from the file
data_file = "/Users/lh270370/Software/nptool/Outputs/Analysis/Isolde/TEST/RdmCoinc_alpha_efficiency_ALL.txt"  # Replace with the path to your file

# Initialize lists to hold the data
detector_pairs = []
observed_counts = []
observed_errors = []
expected_counts = []
expected_errors = []

# Read and parse the data from the file
with open(data_file, 'r') as file:
    for line in file:
        if 'Detector Pair' in line:
            pair_info = line.split('(')[-1].split(')')[0].replace('Detectors: ', '').strip()
            pair = pair_info.split(', ')
            detector_pairs.append((pair[0], pair[1]))
        elif 'Observed Counts' in line:
            observed, error = line.split(':')[1].split(' +/- ')
            observed_counts.append(float(observed.strip()))
            observed_errors.append(float(error.strip()))
        elif 'Expected Counts' in line:
            expected, error = line.split(':')[1].split(' +/- ')
            expected_counts.append(float(expected.strip()))
            expected_errors.append(float(error.strip()))

# Convert lists to numpy arrays for plotting
x = np.arange(len(detector_pairs))  # the label locations

# Plotting
plt.figure(figsize=(10, 6))

# Observed Counts Bar Graph
plt.bar(x, observed_counts, yerr=observed_errors, capsize=5, label='Observed Counts', color='skyblue', alpha=1)

# Expected Counts Plot on top of the bars
plt.errorbar(x, expected_counts, yerr=expected_errors, fmt='^', capsize=5, label='Expected Counts', color='salmon', linestyle='None')

print(np.array(observed_counts)/np.array(expected_counts))
# Labels and Title
plt.ylabel('Counts')
plt.title('Observed vs Expected Counts with Uncertainties')
plt.xticks(x, [f'({a}, {b})' for a, b in detector_pairs])
plt.legend()
plt.grid(True, linestyle='--', alpha=0.5)

plt.tight_layout()
plt.show()
