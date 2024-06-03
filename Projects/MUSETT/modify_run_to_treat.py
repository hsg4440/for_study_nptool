#!/opt/homebrew/opt/python@3.11/libexec/bin/python
import argparse
import os
import subprocess


def parse_arguments():
    parser = argparse.ArgumentParser(description='Modify IsoldeMerging2.hh based on merge_config.txt and command-line arguments.')
    parser.add_argument("-F", "--sim_output_file", required=True, help="Name of the simulation output file")
    parser.add_argument("-D", "--sim_output_dir", required=True, help="Name of the simulation output directory")
    return parser.parse_args()

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



    sim_dir = "/Users/lh270370/Software/nptool/Outputs/Simulation/Isolde"
    config_path ="configs/merge_config.txt"
    header_path ="Scripts/IsoldeMerging3.hh"
    cxx_path = "IsoldeMerging3.cc"

    sim_dir = sim_dir+'/'+args.sim_output_dir
    sim_file = args.sim_output_file
    new_path = f"{sim_dir}/{sim_file}.root"
    replace_line_in_file('RunToTreat.txt', 'RootFileName', f'\t{new_path}')
    replace_line_in_file('RunToTreat.txt', 'TTreeName', '\tSimulatedTree')


if __name__ == "__main__":
    main()
