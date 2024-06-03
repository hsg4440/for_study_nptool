#!/bin/bash

# -F is for name of the output file
# -O is for name of the output directory



SIM_DIR="/Users/lh270370/Software/nptool/Outputs/Simulation/Isolde"
OUTPUT_DIR="default"
OUTPUT_FILE="output"

MACRO="macros/run1e4.mac"
NSIMU=3

DETECTOR_CONFIG="DetectorConfiguration/MUSETT_NPS.detector"
EVENT_SOURCE="Beam/alpha7_5.source"

# Process command-line options
while getopts "F:D:B:E:N:" opt; do
  case $opt in
    F) OUTPUT_FILE="$OPTARG"
    ;;
    D) OUTPUT_DIR="$OPTARG"
    ;;
    B) MACRO="macros/${OPTARG}.mac"  # Set MACRO based on the provided argument
    ;;
    E) EVENT_SOURCE="Beam/$OPTARG"  # Set MACRO based on the provided argument
    ;;
    N) NSIMU=$OPTARG  # Set MACRO based on the provided argument
    ;;
    \?) echo "Invalid option -$OPTARG" >&2
    ;;
  esac
done

# Ensure the output directory exists
if [ ! -d "$SIM_DIR/$OUTPUT_DIR" ]; then
    mkdir -p "$SIM_DIR/$OUTPUT_DIR"
fi

# Array to store the paths of all the output files
output_files=()
for i in $(seq 1 $NSIMU)
do
  FILE="${OUTPUT_FILE}_${i}"  # Fixed concatenation for output file name
  output_files+=("$SIM_DIR/$OUTPUT_DIR/$FILE.root")
done

#Check if output file does not exist already
final_output_file="$SIM_DIR/$OUTPUT_DIR/${OUTPUT_FILE}.root"
if [ -f "$final_output_file" ]; then
    # File exists, prompt the user
    echo "$final_output_file exists."
    read -p "Do you want to remove all existing output files in the directory and proceed? (y/n): " -r
    echo    # Move to a new line
    if [[ $REPLY =~ ^[Yy]$ ]]; then
        # User chose to overwrite, so remove the files
        echo "Removing existing output files..."
        for file in "${output_files[@]}"; do
            rm -f "$file"
        done
        rm -f "$final_output_file"
        echo "Removed existing files. Proceeding with the simulation..."
    else
        # User chose not to overwrite, exit the script
        echo "Not overwriting the files. Exiting."
        exit 1
    fi
fi


# Run several simulations
for i in $(seq 1 $NSIMU)
do
    let "rdm_seed=i*61"

    FILE="${OUTPUT_FILE}_${i}"  # Fixed concatenation for output file name

    if [ $i -eq $NSIMU ]; then
      npsimulation -D "$DETECTOR_CONFIG" -E "$EVENT_SOURCE" -O "Isolde/$OUTPUT_DIR/$FILE" -B "$MACRO" --random-seed $rdm_seed
    else
      nohup npsimulation -D "$DETECTOR_CONFIG" -E "$EVENT_SOURCE" -O "Isolde/$OUTPUT_DIR/$FILE" -B "$MACRO" --random-seed $rdm_seed &
      pids+=($!)
    fi
done

#echo "${output_files[@]}"
for pid in ${pids[@]}; do
    wait $pid
done

hadd "$SIM_DIR/$OUTPUT_DIR/${OUTPUT_FILE}.root" "${output_files[@]}"
rm "${output_files[@]}"
