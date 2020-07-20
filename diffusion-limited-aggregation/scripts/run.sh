#!/bin/bash

# Calculates the mean and standart deviation from a set of Purkinje networks

SEEDS=( 1562002894 1562005513 1562005553 1562006177 1562007596 1562008172 1562008424 1562009134 1562009769 1562002891 )

for SEED in "${SEEDS[@]}"; do
	python mean_std.py /home/berg/Github/Network-Generator/diffusion-limited-aggregation/outputs/01_SRN/03_DLA/seed:${SEED}/segments_length.dat
done

echo

for SEED in "${SEEDS[@]}"; do
    python mean_std.py /home/berg/Github/Network-Generator/diffusion-limited-aggregation/outputs/01_SRN/03_DLA/seed:${SEED}/bifurcations_angle.dat
done

