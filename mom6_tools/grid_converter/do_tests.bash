#!/bin/bash

function ensure_nc_files_match() {
    # Uses the "nccmp" tool to compare two NetCDF files and ensure the
    # metadata and data values match. (See http://nccmp.sourceforge.net/)

    NCFILE1="$1"
    NCFILE2="$2"

    if [ -z "${NCFILE1}" -o -z "${NCFILE2}" ] ; then
        echo "Usage: $FUNCNAME NCFILE1 NCFILE2"
        echo "Ensure the two NetCDF files match"
        exit 1
    fi

    # Find all the metadata diffrences (except history)
    # "--warn=eos" to ignore the null at the end of empty strings
    # (the original files don't have it, but the new ones do)
    nccmp --metadata --global --globalex code_version --force --warn=eos "${NCFILE1}" "${NCFILE2}"
    status=$?

    # Indicate if the data differs, without printing out every different value
    # NOTE: --Tolerance (captial T) value is in percentage points;
    #       Allowing 0.051% difference allows the xgrid_area (in the
    #       atmos/land/ocean tileXtile files) values from the original
    #       'make_quick_mosaic' tool to match "close enough" to the
    #       ones calculated from the supergrid
    nccmp --data --Tolerance=0.051 --force "${NCFILE1}" "${NCFILE2}"
    status2=$?

    # For nccmp, a return status of 0 means identical; 1 means not
    # identical; >1 means error.  Since we ran it twice, return the
    # largest status value.
    if [ $status2 -gt $status ] ; then
        status=$status2
    fi

    if [ $status -eq 0 ] ; then
	echo "PASSED: ${NCFILE1} vs ${NCFILE2}"
    else
	echo "FAILED ($status): ${NCFILE1} vs ${NCFILE2}"
    fi
    return $status
}

NC_FILES="ocean_hgrid.nc ocean_topog.nc ocean_mosaic.nc land_mask.nc ocean_mask.nc atmos_mosaic_tile1Xland_mosaic_tile1.nc atmos_mosaic_tile1Xocean_mosaic_tile1.nc land_mosaic_tile1Xocean_mosaic_tile1.nc mosaic.nc"

CCS_GRID_FILE="../../../ROMS-Inputs/CCS1/grid/CCS_7k_0-360_fred_grd.nc"

rm -f ${NC_FILES}
./convert_ROMS_grid_to_MOM6.py "${CCS_GRID_FILE}"
for file in ${NC_FILES}; do
    ensure_nc_files_match expected_CCS/${file} ./${file}
done

SUPERCRITICAL_GRID_FILE="expected_Supercritical/grid_Supercritical.nc"

rm -f ${NC_FILES}
./convert_ROMS_grid_to_MOM6.py "${SUPERCRITICAL_GRID_FILE}"
for file in ${NC_FILES}; do
    ensure_nc_files_match expected_Supercritical/${file} ./${file}
done
