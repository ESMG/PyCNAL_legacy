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

    file=`basename "${NCFILE1}"`
    if [ $status -eq 0 ] ; then
	echo "PASSED: ${file}"
    else
	echo "FAILED ($status): ${file}"
    fi
    return $status
}

function run_converter() {
    GRID_FILE="$1"
    OUTPUT_DIR="$2"

    [ -d "${OUTPUT_DIR}" ] && rm -rf "${OUTPUT_DIR}"
    mkdir "${OUTPUT_DIR}"
    pushd "${OUTPUT_DIR}" > /dev/null
    ../convert_ROMS_grid_to_MOM6.py "${GRID_FILE}"
    popd > /dev/null
}

function verify_output() {
    EXPECTED_DIR="$1"
    OUTPUT_DIR="$2"

    echo "===== ${EXPECTED_DIR} vs ${OUTPUT_DIR} ====="
    for file in "${EXPECTED_DIR}"/*; do
	file=`basename "${file}"`
        ensure_nc_files_match "${EXPECTED_DIR}/${file}" "${OUTPUT_DIR}/${file}"
    done
}

# CCS
OUT_DIR="CCS_grid"
run_converter "../../../../ROMS-Inputs/CCS1/grid/CCS_7k_0-360_fred_grd.nc" "${OUT_DIR}"
verify_output "expected/CCS" "${OUT_DIR}"

# Supercritical
OUT_DIR="Supercritical_grid"
run_converter "../expected/grid_Supercritical.nc" "${OUT_DIR}"
verify_output "expected/Supercritical" "${OUT_DIR}"
