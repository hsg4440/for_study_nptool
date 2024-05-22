#!/bin/bash

for i in {1..10}
do
 npsimulation -D Vendeta_inelastic.detector -E 238U_el.reaction --random-seed $i -O vendeta_el_$i -B run.mac &
done

for i in {1..10}
do
 npsimulation -D Vendeta_inelastic.detector -E 238U_1st_inel.reaction --random-seed $i -O vendeta_1st_inel_$i -B run_1st.mac &
done

