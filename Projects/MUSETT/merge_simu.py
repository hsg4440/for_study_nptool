#!/opt/homebrew/opt/python@3.11/libexec/bin/python

import argparse
import os
import subprocess


def parse_arguments():
    parser = argparse.ArgumentParser(description='Modify IsoldeMerging2.hh based on merge_config.txt and command-line arguments.')
    parser.add_argument("-F", "--sim_output_file", required=True, help="Name of the simulation output file")
    parser.add_argument("-D", "--sim_output_dir", required=True, help="Name of the simulation output directory")
    parser.add_argument("-R", "--change_RunToTreat", required=True, help="1 To change RunToTreat")
    return parser.parse_args()

def parse_config(config_path):
    config_vars = {}
    with open(config_path, 'r') as config_file:
        for line in config_file:
            if '=' in line:
                var, value = line.split('=')
                config_vars[var.strip()] = value.strip()
    return config_vars

def update_config(header_path, config_vars, sim_output_dir, sim_output_file):
    temp_file_path = header_path + '.tmp'
    t_max = config_vars['t_max']
    T_ADC = config_vars['T_ADC']
    coinc_time = config_vars['coinc_time']
    with open(header_path, 'r') as merging_file, open(temp_file_path, 'w') as temp_file:
        for line in merging_file:
            modified = False

            if 'FileName' in line:
                temp_file.write(f'std::string FileName = "{sim_output_dir}/{sim_output_file}";\n')
                modified = True

            if 'double t_max =' in line :
                temp_file.write(f'double t_max = {t_max};\n')
                modified = True;
            elif 'double coinc_time = ' in line :
                temp_file.write(f'double coinc_time = {coinc_time};\n')
                modified = True;
            elif 'double T_ADC = ' in line :
                temp_file.write(f'double T_ADC = {T_ADC};\n')
                modified = True;
            if not modified:
                temp_file.write(line)
    os.replace(temp_file_path, header_path)

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


def run_root_script(cxx_path):
    # Path to thisroot.sh (adjust as necessary for your ROOT installation)
    bash = "~/.bashrc"
    command = f"cd Scripts && source {bash} && root -l  -q {cxx_path}++"

    try:
        # Execute the command in bash to allow 'source' to work
        subprocess.run(command, shell=True, executable="/bin/bash", check=True)
        print("ROOT script executed successfully.")
        subprocess.run("rm temp.root", shell = True)
    except subprocess.CalledProcessError as e :
        print(f"Failed to execute ROOT script: {e}")


def main():
    args = parse_arguments()



    sim_dir = "/Users/lh270370/Software/nptool/Outputs/Simulation/Isolde"
    config_path ="configs/MergeConfig.txt"
    header_path="Scripts/IsoldeMerging3.hh"
    merging_cc = "IsoldeMerging3.cc"

    config_vars = parse_config(config_path)
    update_config(header_path, config_vars, sim_dir+'/'+args.sim_output_dir, args.sim_output_file)
    run_root_script(merging_cc)


    sim_dir = sim_dir+'/'+args.sim_output_dir
    sim_file = args.sim_output_file
    if args.change_RunToTreat :
        new_path = f"{sim_dir}/{sim_file}_MUSETT_tau" + str(config_vars['coinc_time'])
        new_path += "_T" +str(config_vars['t_max']) + "_TADC" + str(config_vars['T_ADC']) + ".root"

        replace_line_in_file('RunToTreat.txt', 'RootFileName', f'\t{new_path}')
        replace_line_in_file('RunToTreat.txt', 'TTreeName', '\tSimulatedTree')

if __name__ == "__main__":
    main()
