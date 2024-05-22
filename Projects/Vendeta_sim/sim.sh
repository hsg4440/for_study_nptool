#!/bin/bash

for i in {1..2}
do
  npsimulation -D Vendeta_inelastic.detector -E 238Uel.reaction --random-seed $i -O vendeta_sim_$i -B run.mac &
done

