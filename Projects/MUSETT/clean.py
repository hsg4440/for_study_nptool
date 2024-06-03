#!/opt/homebrew/opt/python@3.11/libexec/bin/python

import argparse
import os
import subprocess
import shutil
import glob


def parse_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument("-SF", dest="sim_file", help="Name of the simulation output file", default="default_file")
    parser.add_argument("-SD", dest="sim_dir", help="Name of the simulation output directory", default="default_dir")
    parser.add_argument("-AF", dest="ana_file", help="Name of the analysis output directory")
    parser.add_argument("-AD", dest="ana_dir", help="Name of the analysis output file")
    return parser.parse_args()


def clean_sim(directory,file):
    sim_path = "/Users/lh270370/Software/nptool/Outputs/Simulation/Isolde"
    new_directory = sim_path + "/" + directory + "/" + file
    #print("SIM : ", new_directory)
    file_pattern = sim_path + "/" + directory + "/" + "*" + file + "*"

    if not os.path.exists(new_directory):
        os.makedirs(new_directory)

    for file in glob.glob(file_pattern):
        # Avoid trying to move the directory into itself
        if not os.path.isdir(file):
            shutil.move(file, os.path.join(new_directory, os.path.basename(file)))

def clean_ana(directory,file):
    sim_path = "/Users/lh270370/Software/nptool/Outputs/Analysis/Isolde"
    new_directory = sim_path + "/" + directory + "/" + file
    #print("ANA : ", new_directory)
    file_pattern = sim_path + "/" + directory + "/" + "*" + file + "*"

    if not os.path.exists(new_directory):
        os.makedirs(new_directory)

    for file in glob.glob(file_pattern):
        # Avoid trying to move the directory into itself
        if not os.path.isdir(file):
            shutil.move(file, os.path.join(new_directory, os.path.basename(file)))







def main():
    args = parse_arguments()

    #clean_sim(args.sim_dir,args.sim_file)
    clean_ana(args.ana_dir,args.ana_file)


if __name__ == "__main__":
    main()
