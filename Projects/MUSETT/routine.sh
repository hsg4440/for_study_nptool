#!/bin/bash

DIR="BranchingRatio_test"
SIM_FILE="222Ra_ChangedBR_p3_1e6__energyBeam30"
PHYS_FILE="$SIM_FILE"_physics
MACRO="run2e5"
NSIMU=5
EVENT_SOURCE="Ra222.beam"
CHANGE_RUN_TO_TREAT=1

#npcompilation -a
#./multithread_simu.sh -D "$DIR" -F "$SIM_FILE" -B "$MACRO" -E "$EVENT_SOURCE" -N "$NSIMU"
#./clean.py -SD "$DIR" -SF "$SIM_FILE" -AD "$DIR" -AF "$PHYS_FILE"
./merge_simu.py -D "$DIR" -F "$SIM_FILE" -R "$CHANGE_RUN_TO_TREAT"
#./modify_run_to_treat.py  -D "$DIR/$SIM_FILE" -F "$SIM_FILE"
#./analyse_simu.py -ID "$DIR" -IF "$SIM_FILE" -OF "$PHYS_FILE"
#./sanity_check.py -D "$DIR" -F "$PHYS_FILE"
#./clean_dir.py -SD "$DIR" -SF "$SIM_FILE" -AD "$DIR" -AF "$PHYS_FILE"
#clear
