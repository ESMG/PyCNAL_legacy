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
    nccmp --metadata --global --force --warn=eos "${NCFILE1}" "${NCFILE2}"
    status=$?

    # Indicate if the data differs, without printing out every different value
    nccmp --data --Tolerance=1e-5 --statistics --force --quiet "${NCFILE1}" "${NCFILE2}"
    status2=$?

    # For nccmp, a return status of 0 means identical; 1 means not
    # identical; >1 means error.  Since we ran it twice, return the
    # largest status value.
    if [ $status2 -gt $status ] ; then
        status=$status2
    fi

    if [ $status -eq 0 ] ; then
	echo "PASSED"
    else
	echo "FAILED ($status): ${NCFILE1} ${NCFILE2}"
    fi
    return $status
}

GRID_FILE="../../../ROMS-Inputs/CCS1/grid/CCS_7k_0-360_fred_grd.nc"

NC_FILES="ocean_hgrid.nc ocean_topog.nc ocean_mosaic.nc land_mask.nc ocean_mask.nc atmos_mosaic_tile1Xland_mosaic_tile1.nc atmos_mosaic_tile1Xocean_mosaic_tile1.nc land_mosaic_tile1Xocean_mosaic_tile1.nc mosaic.nc"

rm -f ${NC_FILES}
./convert_ROMS_grid_to_MOM6.py "${GRID_FILE}"
for file in ${NC_FILES}; do
    ensure_nc_files_match expected_output/${file} ./${file}
done
