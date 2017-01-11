#!/bin/bash

module load anaconda

GRID_FILE="/t1/scratch/raphael/CCS1/CCS1-I/CCS_7k_0-360_fred_grd.nc"

./convert_ROMS_grid_to_MOM6.py "${GRID_FILE}"
