#!/bin/bash

echo "=============================================="
echo "Squid giant axon"
./bin/CableEquation 1 30 1 500 25
mv output/propagation_velocity.txt output/velocity_squid_axon.txt
echo "=============================================="

echo "=============================================="
echo "Lobster giant axon"
./bin/CableEquation 2 60 1 75
mv output/propagation_velocity.txt output/velocity_lobster_axon.txt
echo "=============================================="

echo "=============================================="
echo "Crab giant axon"
./bin/CableEquation 7 90 1 30
mv output/propagation_velocity.txt output/velocity_crab_axon.txt
echo "=============================================="

echo "=============================================="
echo "Earthworm giant axon"
./bin/CableEquation 12 200 0.3 105
mv output/propagation_velocity.txt output/velocity_earthworm_axon.txt
echo "=============================================="

echo "=============================================="
echo "Marine worm giant axon"
./bin/CableEquation 1.2 57 0.75 560
mv output/propagation_velocity.txt output/velocity_marine_worm_axon.txt
echo "=============================================="

echo "=============================================="
echo "Mammalian cardiac cell"
./bin/CableEquation 7 150 1.2 20
mv output/propagation_velocity.txt output/velocity_mammalian_cardiac.txt
echo "=============================================="

echo "=============================================="
echo "Barnacle muscle fiber"
./bin/CableEquation 0.23 30 20 400
mv output/propagation_velocity.txt output/velocity_barnacle_muscle.txt
echo "=============================================="
