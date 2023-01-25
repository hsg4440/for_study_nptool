#!/bin/bash

for i in {1..40}
do
  npsimulation -D sofia.detector -E sofia_238U.reaction --random-seed $i -O sofia_simu_238U_$i -B run.mac &
done

