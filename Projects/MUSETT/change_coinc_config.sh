#!/bin/bash

# Define the path to the configuration file
CONFIG_FILE="configs/merge_config.txt"

# Define the Python file to store the results
PYTHON_FILE="store_results3.py"

# Initialize the Python file
echo "results = []" > $PYTHON_FILE

# Define an array of values for coinc_time
COINC_TIMES=(1 2 5 10 20 50 100 200 500 1000 2000 5000)

# Loop over each value in the COINC_TIMES array
for coinc_time in "${COINC_TIMES[@]}"; do
    # Backup the original configuration file (optional)
    cp $CONFIG_FILE "${CONFIG_FILE}.bak"

    # Use 'sed' to replace the coinc_time value in the configuration file
    sed -i '' "s/coinc_time=.*/coinc_time= $coinc_time/" $CONFIG_FILE

    # Run the routine and capture the output
    OUTPUT=$(./routine.sh)

    # Extract the required numbers using awk
    measured=$(echo "$OUTPUT" | awk '/RDM_CNT/{print $3}')  # Assuming the correct format is `RDM_CNT = value`
    expected=$(echo "$OUTPUT" | awk '/calcul/{print $3}')   # Assuming the correct format is `calcul  = value`

    # Append the results to the Python file
    echo "results.append({'coinc_time': $coinc_time, 'measured': $measured, 'expected': $expected})" >> $PYTHON_FILE

    # Restore the original configuration file
    mv "${CONFIG_FILE}.bak" $CONFIG_FILE
done

# Move the results Python file to the specified directory (replace '/desired/path/' with your actual path)
PATH_BENCH_COINC=/Users/lh270370/Software/nptool/Outputs/Simulation/Isolde/bench_coinc/
mv $PYTHON_FILE $PATH_BENCH_COINC$PYTHON_FILE

# Assuming you have a Python script named `plot_results.py` in the same directory as $PYTHON_FILE
# that can plot these results
#cd PATH_BENCH_COINC
#python plot_results3.py

echo "Script execution completed."
