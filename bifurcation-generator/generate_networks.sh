#!/bin/bash

SEGMENT_SIZE=20000

build_network () {
  NUM_BIFF=$1

  echo "[!] Building network !"
  echo "[!] Number of segments leaving bifurcation = $NUM_BIFF"
  ./bin/BifurcationGenerator $SEGMENT_SIZE $NUM_BIFF
  mv output/network.vtk output/network_biff_$NUM_BIFF.vtk
}


for NUM_BIFF in $(seq 1 4); do
  build_network $NUM_BIFF
done
