#!/usr/bin/python2.7
import os, time, sys, multiprocessing, glob

runs = open("/gold/nestar/data/GANIL/e805/MFMRuns/RunRemerged/runlistlist.txt").readlines()
while len(runs) > 0:
  core_tot = multiprocessing.cpu_count()
  core_used = os.getloadavg()[0]
  core_ratio = float(core_used)/float(core_tot)
  mem_str = os.popen("free -m | awk \'/Mem/{print $3}\'").read()
  mem = int(mem_str[:-1])
  irun = ""
  if (mem > 20000) or (core_ratio > 0.75):
    core_tot = multiprocessing.cpu_count()
    core_used = os.getloadavg()[0]
    core_ratio = float(core_used)/float(core_tot)
    mem_str = os.popen("free -m | awk \'/Mem/{print $3}\'").read()
    mem = int(mem_str[:-1])
    sys.stdout.write("\rCurrent Memory occupation : {0}  Current Core occupation : {1}".format(mem,core_ratio))
    sys.stdout.flush()
    time.sleep(1)
  else:
    runs_list = runs.pop(0).strip()
    command = "nohup npmfm/npmfm/install/bin/MFMUnpack -D npmfm/npmfm/examples/E805DataRun/MUGAST.detector -ID /gold/nestar/data/GANIL/e805/MFMRuns/RunRemerged/ -OD /gold/nestar/data/GANIL/e805/NPRootR/RunRemerged/ -exp zdd -R /gold/nestar/data/GANIL/e805/MFMRuns/RunRemerged/"+runs_list+"  > /gold/nestar/data/GANIL/e805/NPRootR/nohup_logs/logs_conversion_"+runs_list+".log 2>&1 </dev/null &"
    print "Starting conversion of runlist "+runs_list
    os.system(command)
    time.sleep(20)


