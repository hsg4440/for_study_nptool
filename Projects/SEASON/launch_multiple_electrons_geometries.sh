#!/bin/bash

# Launch npsimulation

INPUT_File='./Detector/VariableGeom.detector'
SOURCE_FILE='./sources/alpha_electron_cascade_5mm.source'
EMISSION_FILE='./sources/alpha_electron_cascade.txt'
G4_MAC_FILE='./RunCondition/Run.mac'
OUT_DIR='Electron/MultipleEmissionMultipleGeometries'

PARTICLE='electron'
RANDOM=$$

loop_geo() {
  echo "Multiplicity 1" > $EMISSION_FILE
  echo "electron 0.05 0" >> $EMISSION_FILE
  
  for DIST in 2 3 4 5 7 10 15; do
    
    echo "/random/setSeeds " $RANDOM $RANDOM > $G4_MAC_FILE
    echo "/run/initialize" >> $G4_MAC_FILE
    echo "/run/beamOn 100000" >> $G4_MAC_FILE
    echo "exit" >> $G4_MAC_FILE
    
    sed -i -e '12c \  DSSD1Dist= '$DIST' mm' $INPUT_File
    sed -i -e '13c \  TunnelDist= '$DIST' mm' $INPUT_File
    
    npsimulation -E $SOURCE_FILE -D $INPUT_File -O $OUT_DIR"/0Alpha1Electron_"$DIST"mm.root" -B $G4_MAC_FILE
  done
  
  echo "Multiplicity 2" > $EMISSION_FILE
  echo "electron 0.05 0" >> $EMISSION_FILE
  echo "alpha 8 0" >> $EMISSION_FILE
  
  for DIST in 2 3 4 5 7 10 15; do
    
    echo "/random/setSeeds " $RANDOM $RANDOM > $G4_MAC_FILE
    echo "/run/initialize" >> $G4_MAC_FILE
    echo "/run/beamOn 100000" >> $G4_MAC_FILE
    echo "exit" >> $G4_MAC_FILE
    
    sed -i -e '12c \  DSSD1Dist= '$DIST' mm' $INPUT_File
    sed -i -e '13c \  TunnelDist= '$DIST' mm' $INPUT_File
    
    npsimulation -E $SOURCE_FILE -D $INPUT_File -O $OUT_DIR'/1Alpha1Electron_'$DIST'mm.root' -B $G4_MAC_FILE
  done
  
  
  echo "Multiplicity 3" > $EMISSION_FILE
  echo "electron 0.05 0" >> $EMISSION_FILE
  echo "electron 0.15 0" >> $EMISSION_FILE
  echo "alpha 8 0" >> $EMISSION_FILE
  
  for DIST in 2 3 4 5 7 10 15; do
    
    echo "/random/setSeeds " $RANDOM $RANDOM > $G4_MAC_FILE
    echo "/run/initialize" >> $G4_MAC_FILE
    echo "/run/beamOn 100000" >> $G4_MAC_FILE
    echo "exit" >> $G4_MAC_FILE
    
    sed -i -e '12c \  DSSD1Dist= '$DIST' mm' $INPUT_File
    sed -i -e '13c \  TunnelDist= '$DIST' mm' $INPUT_File
    
    npsimulation -E $SOURCE_FILE -D $INPUT_File -O $OUT_DIR"/1Alpha2Electron_"$DIST"mm.root" -B $G4_MAC_FILE
  done
  
  echo "Multiplicity 4" > $EMISSION_FILE
  echo "electron 0.05 0" >> $EMISSION_FILE
  echo "electron 0.15 0" >> $EMISSION_FILE
  echo "electron 0.25 0" >> $EMISSION_FILE
  echo "alpha 8 0" >> $EMISSION_FILE
  
  
  for DIST in 2 3 4 5 7 10 15; do
    
    echo "/random/setSeeds " $RANDOM $RANDOM > $G4_MAC_FILE
    echo "/run/initialize" >> $G4_MAC_FILE
    echo "/run/beamOn 100000" >> $G4_MAC_FILE
    echo "exit" >> $G4_MAC_FILE
    
    sed -i -e '12c \  DSSD1Dist= '$DIST' mm' $INPUT_File
    sed -i -e '13c \  TunnelDist= '$DIST' mm' $INPUT_File
    
    npsimulation -E $SOURCE_FILE -D $INPUT_File -O $OUT_DIR"/1Alpha3Electron_"$DIST"mm.root" -B $G4_MAC_FILE
  done

  echo "Multiplicity 5" > $EMISSION_FILE
  echo "electron 0.05 0" >> $EMISSION_FILE
  echo "electron 0.15 0" >> $EMISSION_FILE
  echo "electron 0.25 0" >> $EMISSION_FILE
  echo "electron 0.5 0" >> $EMISSION_FILE
  echo "alpha 8 0" >> $EMISSION_FILE
  
  for DIST in 2 3 4 5 7 10 15; do
    
    echo "/random/setSeeds " $RANDOM $RANDOM > $G4_MAC_FILE
    echo "/run/initialize" >> $G4_MAC_FILE
    echo "/run/beamOn 100000" >> $G4_MAC_FILE
    echo "exit" >> $G4_MAC_FILE
    
    sed -i -e '12c \  DSSD1Dist= '$DIST' mm' $INPUT_File
    sed -i -e '13c \  TunnelDist= '$DIST' mm' $INPUT_File
    
    npsimulation -E $SOURCE_FILE -D $INPUT_File -O $OUT_DIR"/1Alpha4Electron_"$DIST"mm.root" -B $G4_MAC_FILE
  done

  echo "Multiplicity 6" > $EMISSION_FILE
  echo "electron 0.05 0" >> $EMISSION_FILE
  echo "electron 0.125 0" >> $EMISSION_FILE
  echo "electron 0.15 0" >> $EMISSION_FILE
  echo "electron 0.25 0" >> $EMISSION_FILE
  echo "electron 0.5 0" >> $EMISSION_FILE
  echo "alpha 8 0" >> $EMISSION_FILE
  
  for DIST in 2 3 4 5 7 10 15; do
    
    echo "/random/setSeeds " $RANDOM $RANDOM > $G4_MAC_FILE
    echo "/run/initialize" >> $G4_MAC_FILE
    echo "/run/beamOn 100000" >> $G4_MAC_FILE
    echo "exit" >> $G4_MAC_FILE
    
    sed -i -e '12c \  DSSD1Dist= '$DIST' mm' $INPUT_File
    sed -i -e '13c \  TunnelDist= '$DIST' mm' $INPUT_File
    
    npsimulation -E $SOURCE_FILE -D $INPUT_File -O $OUT_DIR"/1Alpha5Electron_"$DIST"mm.root" -B $G4_MAC_FILE
  done
}

loop_geo

