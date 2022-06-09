#!/bin/bash

rfile='./Reaction/47Kdp_Sim_0143MeV.reaction'
#rfile='./Reaction/47Kdp_Sim_0968MeV.reaction'
#rfile='./Reaction/47Kdp_Sim_3605MeV.reaction'

#rfile='./Reaction/47Kdp_Sim_3605MeV_Flat.reaction'
#rfile='./Reaction/47Kdp_Sim.reaction'
#rfile='./Reaction/IsotropicProtons.reaction'
#rfile='./Reaction/47Kdp_Sim_Flat3500MeV.reaction'

cd ~/Programs/nptool/Projects/e793s;
cmake ./;
make -j6;

directory=' /home/charlottepaxman/Programs/nptool/Outputs/Simulation/'
dotroot='.root'
dash='-'

for x in 1 2 3 4 5
do
        outname=$1$dash$x
	npsimulation -D ./Detector/mugast_08Nov.detector -E $rfile -B ./runsimulation.mac -O $outname;

	filename=$directory$1$dash$x$dotroot

done

sim='Sim_'
outfile=$sim$1

echo "TTreeName" > RunToTreat_AutoGenerated.txt
echo " SimulatedTree" >> RunToTreat_AutoGenerated.txt
echo "RootFileName" >> RunToTreat_AutoGenerated.txt

for x in 1 2 3 4 5
do
	filename=$directory$1$dash$x$dotroot
	echo $filename >> RunToTreat_AutoGenerated.txt
done

npanalysis --definition Sim -R RunToTreat_AutoGenerated.txt -E $rfile -D Detector/mugast_08Nov.detector -O $outfile;
