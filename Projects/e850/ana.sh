#!/bin/bash

#function ana(){
	run_number=$1
	echo "Run Number : $run_number"
	run_name=$(printf 'run_%03d' $run_number)
	echo "Run Name : $run_name"
	
	rm RunToTreat.txt
	printf "TTreeName\n  AD\nRootFileName\n  ../Analysis/RootA/r0%03d_00*a.root" $run_number >> RunToTreat.txt

	npanalysis -R RunToTreat.txt -C Calibration.txt -D pista.detector -O $run_name
#}

"$@"

