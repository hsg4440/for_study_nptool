#!/usr/bin/env python3

import subprocess
import os
import time
from time import sleep
from multiprocessing.pool import Pool
import sys
import subprocess as sp
import subprocess as sp2

Nruns = len(sys.argv)-1

def task(irun, fileN):
  run = int(sys.argv[irun+1])
  with open('.RunToTreat.txt', 'w') as f:
    f.write(str(run)+' '+str(fileN))
  bashcommand= ("root -l -q 'fastFillHistoPISTA.cc(" + str(run) + "," + str(fileN) + ")'")
  os.system(bashcommand)
  f.close()

def hadd():
  for i in range(Nruns):
    run = int(sys.argv[i+1])
    bashcommand= ("hadd -f hist/HistoPISTA_run"+str(run)+".root hist/HistoPISTA_r"+str(run)+"*")
    print(bashcommand)
    os.system(bashcommand)
    bashcommand= ("rm hist/HistoPISTA_r"+str(run)+"*")
    os.system(bashcommand)

if __name__=='__main__':
  NbOfFiles=[0 for i in range(Nruns)]
  for i in range(Nruns):
    run = int(sys.argv[i+1])
    NbOfFiles[i]=sp2.getoutput("ls /home/e850/workdir/Pierre/np_e850/Outputs/Analysis/r"+str('%03d'% run)+"*.root -print | wc -l ")
    print(" run "+ str(run) +" "+ str(NbOfFiles[i])+" files")
  with Pool(15) as pool:
    items=[(r,f) for r in range(int(Nruns)) for f in range(int(NbOfFiles[r]))]
    for result in pool.starmap(task, items):
      print("Got result: {result}", flush=True)
  hadd()
