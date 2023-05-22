#!/bin/bash

for i in {1..10}
do
  npsimulation -D sofia.detector -E sofia.reaction --random-seed $i -O sofia_simu_189Pb_$i -B run.mac &
done

