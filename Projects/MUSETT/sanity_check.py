#!/opt/homebrew/opt/python@3.11/libexec/bin/python

import argparse
import os
import subprocess

def parse_arguments():
    parser = argparse.ArgumentParser(description='Update SanityCheckV2.cxx with new configuration.')
    parser.add_argument("-F", "--phys_file", required=True, help="Name of the physics input file")
    parser.add_argument("-D", "--phys_dir", required=True, help="Name of the physcs input directory")
    return parser.parse_args()

def parse_config(config_path):
    config_vars = {}
    with open(config_path, 'r') as config_file:
        for line in config_file:
            if '=' in line:
                var, value = line.split('=')
                config_vars[var.strip()] = value.strip()
    return config_vars

def modify_cxx_file(cxx_path, config_vars, phys_dir, phys_file):
    input_file_name = f"{phys_dir}/{phys_file}.root"
    output_file_name = f"{phys_dir}/histoSanity_{phys_file}.root"

    with open(cxx_path, 'r') as file:
        lines = file.readlines()

    with open(cxx_path, 'w') as file:
        for line in lines:
            if 'config.inputFileName =' in line:
                file.write(f'    config.inputFileName = "{input_file_name}";\n')
            elif 'config.outputFileName =' in line:
                file.write(f'    config.outputFileName = "{output_file_name}";\n')
            elif 'config.tau =' in line:
                file.write(f'    config.tau = {config_vars["coinc_time"]};\n')
            elif 'config.T =' in line:
                file.write(f'    config.T = {config_vars["t_max"]};\n')
            else:
                file.write(line)

def run_root_script(cxx_path):
    # Path to thisroot.sh (adjust as necessary for your ROOT installation)
    bash = "~/.bashrc"
    command = f"source {bash} && root -l -b -q {cxx_path}"

    try:
        # Execute the command in bash to allow 'source' to work
        subprocess.run(command, shell=True, executable="/bin/bash", check=True)
        print("ROOT script executed successfully.")
    except subprocess.CalledProcessError as e:
        print(f"Failed to execute ROOT script: {e}")


def main():
    args = parse_arguments()

    ana_dir = "/Users/lh270370/Software/nptool/Outputs/Analysis/Simu/Isolde/"
    config_path = "configs/merge_config.txt"
    cxx_path = "Scripts/SanityCheckV2.cxx"

    config_vars = parse_config(config_path)
    modify_cxx_file(cxx_path, config_vars, ana_dir+args.phys_dir, args.phys_file)

    run_root_script(cxx_path)

if __name__ == "__main__":
    main()
