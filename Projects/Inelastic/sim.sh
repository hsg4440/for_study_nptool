#!/bin/bash

for i in {1..10}
do
  npsimulation -D Vendeta_inelastic.detector -E 238Uel.reaction --random-seed $i -O vendeta_el_$i -B run.mac &
done

