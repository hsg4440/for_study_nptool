#!/opt/homebrew/opt/python@3.11/libexec/bin/python
import argparse
import os
import subprocess
import sys
import time
import shutil

def parse_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument("-IF", dest="input_file", help="Name of the file to analyse", default="default_file")
    parser.add_argument("-OF", dest="output_file", help="Name of the output file ", default="default_file")
    parser.add_argument("-OD", dest="output_dir", help="Name of the output directory ", default="default_dir")
    parser.add_argument("-C", dest="calib_file", help="Name of the calibration file", default="CalibrationSimu")
    return parser.parse_args()

def ensure_directory_exists(directory):
    os.makedirs(directory, exist_ok=True)

def check_and_handle_existing_file(file_path):
    if os.path.isfile(file_path):
        print(f"{file_path} exists.")
        user_input = input("Do you want to remove the existing output file and proceed? (y/n): ")
        if user_input.lower() == 'y':
            print("Removing existing output file...")
            os.remove(file_path)
            print("Removed existing file. Proceeding with the analysis...")
        else:
            print("Not overwriting the files. Exiting.")
            sys.exit(1)

def run_command(command):
    print(f'Executing command: {command}')
    process = subprocess.run(command, shell=True, text=True)
    if process.returncode != 0:
        print(f"Error executing command: {command}")
        sys.exit(1)

def replace_line_in_file(file_path, identifier, new_line):
    with open(file_path, 'r') as file:
        lines = file.readlines()

    with open(file_path, 'w') as file:
        write_next = False
        for line in lines:
            if write_next:
                file.write(new_line + '\n')
                write_next = False
            else:
                file.write(line)
            if identifier in line:
                write_next = True

def main():
    args = parse_arguments()

    input_file = args.input_file
    output_file = args.output_file
    output_dir  = args.output_dir
    calib = args.calib_file

    detector_config = "DetectorConfiguration/MUSETT_NPS.detector"
    calibration="calibration/" + calib + ".txt"
    run_to_treat ="RunToTreat.txt"

    time.sleep(1)

    input_dir = "RootR"
    ana_dir = "/Users/lh270370/Software/nptool/Outputs/Analysis/Data"
    raw_output_dir = "/Users/lh270370/Software/nptool/Outputs/Analysis"
    total_dir = os.path.join(ana_dir, output_dir)
    final_output_file = os.path.join(output_dir, f"{output_file}.root")
    final_output_file_path = os.path.join(total_dir, f"{output_file}.root")
    ensure_directory_exists(total_dir)
    check_and_handle_existing_file(final_output_file_path)

    new_path = f"{ana_dir}/{output_dir}/{output_file}.root"

    replace_line_in_file('RunToTreat.txt', 'TTreeName', '\tRD')
    replace_line_in_file('RunToTreat.txt', 'RootFileName', '\t' + input_dir+'/'+input_file)

    if shutil.which("npanalysis"):

        bash = "~/.bashrc"
        #run_command(command)
        run_command(f"source {bash} && npanalysis -D {detector_config} -O {output_file} -C {calibration} -R {run_to_treat}")
        shutil.move(raw_output_dir+f"/{output_file}.root", final_output_file_path)
        #print(f"Moved output file to: {final_output_file_path}")
    else:
        print("Error: npanalysis command not found. Make sure it's installed and in your PATH.")
        sys.exit(1)

if __name__ == "__main__":
    main()
