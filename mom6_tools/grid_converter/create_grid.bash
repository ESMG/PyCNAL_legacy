#!/bin/bash

module load anaconda

FRETOOLS=/t3/workdir/raphael/fre-nctools/mom/src/tools

ln -s /t1/scratch/raphael/CCS1/CCS1-I/CCS_7k_0-360_fred_grd.nc .

python CreateFMSgridTopo.py

module unload anaconda

module load intel/14.0.0
module load netcdf/4.3.0_intel14.0.0

$FRETOOLS/make_solo_mosaic/make_solo_mosaic --num_tiles 1 --dir ./ --mosaic_name ocean_mosaic --tile_file ocean_hgrid.nc
$FRETOOLS/make_quick_mosaic/make_quick_mosaic --input_mosaic ocean_mosaic.nc --ocean_topog ocean_topog.nc
