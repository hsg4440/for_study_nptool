#!/bin/bash

# Launch npsimulation

INPUT_File='./Detector/VariableGeom.detector'
SOURCE_FILE='./sources/electron.source'
G4_MAC_FILE='./RunCondition/run_electron.mac'
OUT_DIR='Electron/FinalElectronGeometries'

PARTICLE='electron'
RANDOM=$$

loop_geo() {
  
  for DIST in 2 3 4 5 6 7 8 9 10 12 15 20 25 30 40 50; do
    
    echo "/random/setSeeds " $RANDOM $RANDOM > $G4_MAC_FILE
    echo "/run/initialize" >> $G4_MAC_FILE
    echo "/run/beamOn 100000" >> $G4_MAC_FILE
    echo "exit" >> $G4_MAC_FILE
    
    sed -i -e '12c \  DSSD1Dist= '$DIST' mm' $INPUT_File
    sed -i -e '13c \  TunnelDist= '$DIST' mm' $INPUT_File
    
    npsimulation -E $SOURCE_FILE -D $INPUT_File -O $OUT_DIR'/Electron_'$DIST'mm.root' -B $G4_MAC_FILE
  done
}

loop_geo

