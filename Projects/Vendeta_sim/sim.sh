#!/bin/bash

for i in {1..2}
do
  npsimulation -D Vendeta.detector -E neutron.source --random-seed $i -O vendeta_sim_$i -B run.mac &
done

