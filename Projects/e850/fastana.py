#!/usr/bin/env python3

import subprocess
import os
import time
from time import sleep
from multiprocessing.pool import Pool
import sys
import subprocess as sp
import subprocess as sp2

start = int(sys.argv[1])
#end = int(sys.argv[2])
Nruns = len(sys.argv)-1

def task(irun, fileN):
  run = int(sys.argv[irun+1])
  infile = "/home/e850/workdir/Pierre/Analysis/RootA/r"+str('%04d'% run)+"_"+str('%03d'% fileN)+"a.root"
  outfile = "r"+str('%03d'% run)+"_"+str('%02d'% fileN)+".root"
  bashcommand= ("npanalysis -T "+infile+ " AD -C Calibration.txt -D pista.detector -O "+outfile)
  logfile = open("log/npanalysis_run"+str(run)+"_f"+str(fileN)+".log",'w')
  print(bashcommand)
  time.sleep(int(fileN+2))
  #subprocess.run(bashcommand.split(), stdout=logfile)
  os.system(bashcommand)
  logfile.close()

def hadd():
  for i in range(Nruns):
    run = int(sys.argv[i+1])
    bashcommand= ("hadd -f Outputs/Analysis/run_"+str('%03d'% run)+".root Outputs/Analysis/r"+str('%03d'% run)+"*")
    print(bashcommand)
    os.system(bashcommand)

if __name__=='__main__':
  NbOfFiles=[0 for i in range(Nruns)]
  for i in range(Nruns):
    run = int(sys.argv[i+1])
    NbOfFiles[i]=sp2.getoutput("ls /home/e850/workdir/Pierre/Analysis/RootA/r"+str('%04d'% run)+"*.root -print | wc -l ")
    print(str(i)+" run "+ str(run) +" "+ str(NbOfFiles[i])+" files")
  with Pool(8) as pool:
    items=[(r,f) for r in range(int(Nruns)) for f in range(int(NbOfFiles[r]))]
    for result in pool.starmap(task, items):
      print("Got result: {result}", flush=True)
  hadd()
