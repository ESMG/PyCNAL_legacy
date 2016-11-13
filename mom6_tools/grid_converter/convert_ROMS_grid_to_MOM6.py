#!/bin/env python
# -*- coding: utf-8 -*-
# vim: set fileencoding=utf-8 :

import sys
import numpy
import netCDF4
import os.path
import datetime
import subprocess

import Spherical

# This script converts a ROMS horizontal grid into a MOM6 horizontal
# grid.

# Both ROMS and MOM6 horizontal grids use an Arakawa C-grid, with four
# types of points:
#   rho: the centers of the cells
#   psi: the corners of the cells, located diagonally between the
#        'rho' points
#   u:   the u-velocity points, located between 'rho' points in the
#        east/west direction
#   v:   the v-velocity points, located between 'rho' points in the
#        north/south direction

# The main differences between the two grids are:
#  * the outermost points of the ROMS grid are the 'rho' points, while
#    the outermost points of the MOM6 grid are the 'psi' points (both
#    with interspersed 'u' and 'v' points); and
#  * the MOM6 grid interleaves all four types of points into a single
#    "supergrid", while ROMS stores them as separate grids.

# The ROMS grid looks like this, with an extra layer of 'rho' points
# around the outside:
# (see https://www.myroms.org/wiki/Numerical_Solution_Technique)
#
#       p - p - p - p - p
#    3  | + | + | + | + |     p = rho (center) points
#       p - p - p - p - p     + = psi (corner) points
#    2  | + | + | + | + |     - = u points
#       p - p - p - p - p     | = v points
#    1  | + | + | + | + |
#       p - p - p - p - p
#
#         1   2   3   4

# The MOM6 grid has 'psi' points on the outside, not 'rho':
#
#    3    + | + | + | +       p = rho (center) points
#         - p - p - p -       + = psi (corner) points
#    2    + | + | + | +       - = u points
#         - p - p - p -       | = v points
#    1    + | + | + | +
#
#         1   2   3   4

def read_ROMS_grid(roms_grid_filename):
    """Load the ROMS grid from a NetCDF file."""

    subgrids = ['psi', 'rho', 'u', 'v'] # also available: 'vert'
    fields = ['lon', 'lat', 'x', 'y', 'mask']

    roms_grid = dict()
    with netCDF4.Dataset(roms_grid_filename) as roms_ds:
        # extract fields named based on which grid they are on
        for subgrid in subgrids:
            roms_grid[subgrid] = dict()
            for field in fields:
                var_name = field + '_' + subgrid
                roms_grid[subgrid][field] = roms_ds.variables[var_name][:]
                if (field == 'x') or (field == 'y'):
                    units = roms_ds.variables[var_name].units.lower()
                    assert units.startswith('meter')
                elif (field == 'lat') or (field == 'lon'):
                    units = roms_ds.variables[var_name].units.lower()
                    assert units.startswith('degree')

        # extract fields that don't follow the above naming pattern
        roms_grid['rho']['h'] = roms_ds.variables['h'][:] # on the rho grid, but not called "h_rho"

        roms_grid['metadata'] = dict()

        spherical = roms_ds.variables['spherical'][:]
        if (spherical == 0) or (spherical == 'F') or (spherical == 'f'):
            roms_grid['metadata']['is_spherical'] = False
        elif (spherical == 1) or (spherical == 'T') or (spherical == 't'):
            roms_grid['metadata']['is_spherical'] = True
        else:
            warn('Unrecognized value for spherical in ROMS grid: %s', str(spherical))
        print 'spherical = ', roms_grid['metadata']['is_spherical']

    return roms_grid

def trim_ROMS_grid(old_grid):
    """Remove extraneous points on the outside of the ROMS grid."""

    trim_subgrid = dict()
    # remove the outer:   ( rows,  cols)
    trim_subgrid['psi'] = (False, False) # Cell corners (leave alone)
    trim_subgrid['rho'] = ( True,  True) # Cell centers (remove outer row and column)
    trim_subgrid[ 'u' ] = ( True, False) # U-points (remove outer row)
    trim_subgrid[ 'v' ] = (False,  True) # V-points (remove outer column)

    new_grid = dict()
    for subgrid in old_grid.keys():
        if subgrid == 'metadata':
            new_grid[subgrid] = dict(old_grid[subgrid])
            continue
        new_grid[subgrid] = dict()
        trim_rows,trim_cols = trim_subgrid[subgrid]
        for field in old_grid[subgrid].keys():
            if trim_rows and trim_cols:
                new_grid[subgrid][field] = old_grid[subgrid][field][1:-1,1:-1]
            elif trim_rows:
                new_grid[subgrid][field] = old_grid[subgrid][field][1:-1, :  ]
            elif trim_cols:
                new_grid[subgrid][field] = old_grid[subgrid][field][ :  ,1:-1]
            else:
                new_grid[subgrid][field] = old_grid[subgrid][field][ :  , :  ]

    return new_grid

def get_git_repo_version_info():
    """Describe the current version of this script as known by Git."""
    repo_name = 'ESMG/PyCNAL'
    git_describe =  subprocess.check_output(['git', 'describe', '--all', '--long', '--dirty', '--abbrev=10']).rstrip()
    return repo_name + ': ' + git_describe

def get_history_entry(argv):
    """Construct an entry for the global 'history' attribute of a NetCDF file,
    which is a date and the command used."""
    today = datetime.datetime.strftime(datetime.datetime.now(), '%Y-%m-%d')
    command = ' '.join(argv)
    return today + ': ' + command

def setup_MOM6_grid(argv):
    mom6_grid = dict()

    mom6_grid['filenames'] = dict()
    mom6_grid['filenames']['supergrid']  = 'ocean_hgrid.nc'
    mom6_grid['filenames']['topography'] = 'ocean_topog.nc'
    mom6_grid['filenames']['mosaic']     = 'ocean_mosaic.nc'
    mom6_grid['filenames']['land_mask']  = 'land_mask.nc'
    mom6_grid['filenames']['ocean_mask'] = 'ocean_mask.nc'

    mom6_grid['netcdf_info'] = dict()
    mom6_grid['netcdf_info']['tile_str']      = 'tile1'
    mom6_grid['netcdf_info']['string_length'] = 255
    mom6_grid['netcdf_info']['grid_version']  = '0.2' # taken from make_solo_mosaic
    mom6_grid['netcdf_info']['code_version']  = '$Name: tikal $' ### for testing
    #mom6_grid['netcdf_info']['code_version']  = get_git_repo_version_info() ### for production
    mom6_grid['netcdf_info']['history_entry'] = get_history_entry(argv)

    mom6_grid['supergrid'] = dict()

    return mom6_grid

def convert_ROMS_to_MOM6(mom6_grid, roms_grid):
    """Convert the ROMS grid data into a skeleton MOM6 grid, mainly by
    merging the four sets of point locations from the ROMS grid
    into a single supergrid for MOM6."""

    # Double the size of the *trimmed* ROMS grid to merge the four
    # sets of points.
    num_rows, num_cols = roms_grid['rho']['lon'].shape
    num_rows *= 2
    num_cols *= 2

    # Store the size for later
    mom6_grid['num_rows'] = num_rows
    mom6_grid['num_cols'] = num_cols

    if roms_grid['metadata']['is_spherical']:
        copy_fields = ['lon', 'lat']
    else:
        copy_fields = ['x', 'y']

    # Copy points from ROMS grid
    for field in copy_fields:
        mom6_grid['supergrid'][field] = numpy.zeros((num_rows+1,num_cols+1))
        mom6_grid['supergrid'][field][ ::2, ::2] = roms_grid['psi'][field] # outer
        mom6_grid['supergrid'][field][1::2,1::2] = roms_grid['rho'][field] # inner
        mom6_grid['supergrid'][field][1::2, ::2] = roms_grid[ 'u' ][field] # between e/w
        mom6_grid['supergrid'][field][ ::2,1::2] = roms_grid[ 'v' ][field] # between n/s

    mom6_grid['depth'] = roms_grid['rho']['h'] * roms_grid['rho']['mask']
    mom6_grid['mask']  = roms_grid['rho']['mask']

    return mom6_grid

def _fill_in_MOM6_grid_metrics_spherical(mom6_grid):
    """Fill in missing MOM6 supergrid metrics by computing best guess
    values based on latitude and longitude coordinates."""

    lat = mom6_grid['supergrid']['lat']
    lon = mom6_grid['supergrid']['lon']

    # Approximate edge lengths as great arcs
    R = 6370.e3 # Radius of sphere
    mom6_grid['supergrid']['dx'][:,:] = R * Spherical.angle_through_center( (lat[ :,1:],lon[ :,1:]), (lat[:  ,:-1],lon[:  ,:-1]) )
    mom6_grid['supergrid']['dy'][:,:] = R * Spherical.angle_through_center( (lat[1:, :],lon[1:, :]), (lat[:-1,:  ],lon[:-1,:  ]) )

    # Approximate angles using centered differences in interior, and side differences on left/right edges
    # TODO: Why do something different at the edges when we have extra ROMS points available?
    cos_lat = numpy.cos(numpy.radians(lat))
    mom6_grid['supergrid']['angle'][:,1:-1] = numpy.arctan( (lat[:,2:] - lat[:,:-2]) / ((lon[:,2:] - lon[:,:-2]) * cos_lat[:,1:-1]) )
    mom6_grid['supergrid']['angle'][:, 0  ] = numpy.arctan( (lat[:, 1] - lat[:, 0 ]) / ((lon[:, 1] - lon[:, 0 ]) * cos_lat[:, 0  ]) )
    mom6_grid['supergrid']['angle'][:,-1  ] = numpy.arctan( (lat[:,-1] - lat[:,-2 ]) / ((lon[:,-1] - lon[:,-2 ]) * cos_lat[:,-1  ]) )

    # Approximate cell areas as that of spherical polygon
    mom6_grid['supergrid']['area'][:,:] = R * R * Spherical.quad_area(lat, lon)

    return mom6_grid

def _fill_in_MOM6_grid_metrics_cartesian(mom6_grid):
    """Fill in missing MOM6 supergrid metrics by computing best guess
    values based on x and y coordinates."""

    x = mom6_grid['supergrid']['x']
    y = mom6_grid['supergrid']['y']

    # Compute edge lengths
    mom6_grid['supergrid']['dx'][:,:] = numpy.sqrt( (x[:,1:] - x[:,:-1])**2 + (y[:,1:] - y[:,:-1])**2 )
    mom6_grid['supergrid']['dy'][:,:] = numpy.sqrt( (x[1:,:] - x[:-1,:])**2 + (y[1:,:] - y[:-1,:])**2 )

    # Compute angles using centered differences in interior, and side differences on left/right edges
    # TODO: Why do something different at the edges when we have extra ROMS points available?
    mom6_grid['supergrid']['angle'][:,1:-1] = numpy.arctan2( (y[:,2:] - y[:,:-2]), (x[:,2:] - x[:,:-2]) )
    mom6_grid['supergrid']['angle'][:, 0  ] = numpy.arctan2( (y[:, 1] - y[:, 0 ]), (x[:, 1] - x[:, 0 ]) )
    mom6_grid['supergrid']['angle'][:,-1  ] = numpy.arctan2( (y[:,-1] - y[:,-2 ]), (x[:,-1] - x[:,-2 ]) )

    # Compute cell areas
    mom6_grid['supergrid']['area'][:,:] = mom6_grid['supergrid']['dx'][:-1, :] * mom6_grid['supergrid']['dy'][:, :-1]

    return mom6_grid

def approximate_MOM6_grid_metrics(mom6_grid):
    """Fill in missing MOM6 supergrid metrics by computing best guess
    values."""

    num_rows = mom6_grid['num_rows']
    num_cols = mom6_grid['num_cols']

    # Declare shapes
    mom6_grid['supergrid']['dx']    = numpy.zeros((num_rows+1,num_cols  ))
    mom6_grid['supergrid']['dy']    = numpy.zeros((num_rows,  num_cols+1))
    mom6_grid['supergrid']['angle'] = numpy.zeros((num_rows+1,num_cols+1))
    mom6_grid['supergrid']['area']  = numpy.zeros((num_rows,  num_cols  ))

    if 'lat' in mom6_grid['supergrid']:
        return _fill_in_MOM6_grid_metrics_spherical(mom6_grid)
    else:
        return _fill_in_MOM6_grid_metrics_cartesian(mom6_grid)

def write_MOM6_supergrid_file(mom6_grid):
    """Save the MOM6 supergrid data into its own file."""
    num_rows, num_cols = mom6_grid['supergrid']['area'].shape

    string_len = len(mom6_grid['netcdf_info']['tile_str'])

    with netCDF4.Dataset(mom6_grid['filenames']['supergrid'], 'w', format='NETCDF3_CLASSIC') as hgrid_ds:
        # Dimensions
        hgrid_ds.createDimension('nx',  num_cols)
        hgrid_ds.createDimension('nxp', num_cols+1)
        hgrid_ds.createDimension('ny',  num_rows)
        hgrid_ds.createDimension('nyp', num_rows+1)
        hgrid_ds.createDimension('string', string_len)

        # Variables & Values
        hx = hgrid_ds.createVariable('x', 'f4', ('nyp','nxp',))
        hy = hgrid_ds.createVariable('y', 'f4', ('nyp','nxp',))

        if 'lon' in mom6_grid['supergrid']:
            hx.units = 'degrees' # TODO: degrees_east?
            hx[:] = mom6_grid['supergrid']['lon']

            hy.units = 'degrees' # TODO: degrees_north?
            hy[:] = mom6_grid['supergrid']['lat']
        else:
            hx.units = 'meters'
            hx[:] = mom6_grid['supergrid']['x']

            hy.units = 'meters'
            hy[:] = mom6_grid['supergrid']['y']

        hdx = hgrid_ds.createVariable('dx', 'f4', ('nyp','nx',))
        hdx.units = 'meters'
        hdx[:] = mom6_grid['supergrid']['dx']

        hdy = hgrid_ds.createVariable('dy', 'f4', ('ny','nxp',))
        hdy.units = 'meters'
        hdy[:] = mom6_grid['supergrid']['dy']

        harea = hgrid_ds.createVariable('area', 'f4', ('ny','nx',))
        harea.units = 'meters^2'
        harea[:] = mom6_grid['supergrid']['area']

        hangle = hgrid_ds.createVariable('angle_dx', 'f4', ('nyp','nxp',))
        hangle.units = 'degrees' # TODO: actually computed in radians!
        hangle[:] = mom6_grid['supergrid']['angle']

        htile = hgrid_ds.createVariable('tile', 'c', ('string',))
        htile[:] = mom6_grid['netcdf_info']['tile_str']

def write_MOM6_topography_file(mom6_grid):
    """Save the MOM6 ocean topography field in a separate file."""
    num_rows, num_cols = mom6_grid['depth'].shape

    with netCDF4.Dataset(mom6_grid['filenames']['topography'], 'w', format='NETCDF3_CLASSIC') as topog_ds:
        # Dimensions
        topog_ds.createDimension('nx', num_cols)
        topog_ds.createDimension('ny', num_rows)
        topog_ds.createDimension('ntiles', 1)

        # Variables & Values
        hdepth = topog_ds.createVariable('depth', 'f4', ('ny','nx',))
        hdepth.units = 'm'
        hdepth[:] = mom6_grid['depth']

def _add_global_attributes(mom6_grid, netcdf_dataset):
    netcdf_dataset.grid_version = mom6_grid['netcdf_info']['grid_version']
    netcdf_dataset.code_version = mom6_grid['netcdf_info']['code_version']
    netcdf_dataset.history      = mom6_grid['netcdf_info']['history_entry']

def write_MOM6_solo_mosaic_file(mom6_grid):
    """Write the "solo mosaic" file, which describes to the FMS infrastructure
     where to find the grid file(s).  Based on tools in version 5 of MOM
     (http://www.mom-ocean.org/)."""

    # NOTE: This function is very basic, since we're skipping the
    # finding of "contact regions" between the tiles that the real
    # make_solo_mosaic tool performs.  It's not needed right now,
    # since we only have one (regional) tile, but I think this feature
    # will be needed if we ever use a tripolar grid.

    with netCDF4.Dataset(mom6_grid['filenames']['mosaic'], 'w', format='NETCDF3_CLASSIC') as mosaic_ds:
        # Dimenisons
        mosaic_ds.createDimension('ntiles', 1)
        mosaic_ds.createDimension('string', mom6_grid['netcdf_info']['string_length'])

        # Variables & Values
        hmosaic = mosaic_ds.createVariable('mosaic', 'c', ('string',))
        hmosaic.standard_name = 'grid_mosaic_spec'
        hmosaic.children = 'gridtiles'
        hmosaic.contact_regions = 'contacts'
        hmosaic.grid_descriptor = ''
        dirname,filename = os.path.split(mom6_grid['filenames']['mosaic'])
        filename,ext = os.path.splitext(filename)
        hmosaic[:len(filename)] = filename

        # Always use './' as the directory since FMS always runs from the
        # "INPUT" directory
        hgridlocation = mosaic_ds.createVariable('gridlocation', 'c', ('string',))
        hgridlocation.standard_name = 'grid_file_location'
        this_dir = './'
        hgridlocation[:len(this_dir)] = this_dir

        hgridfiles = mosaic_ds.createVariable('gridfiles', 'c', ('ntiles', 'string',))
        hgridfiles[0, :len(mom6_grid['filenames']['supergrid'])] = mom6_grid['filenames']['supergrid']

        hgridtiles = mosaic_ds.createVariable('gridtiles', 'c', ('ntiles', 'string',))
        hgridtiles[0, :len(mom6_grid['netcdf_info']['tile_str'])] = mom6_grid['netcdf_info']['tile_str']

        # Global attributes
        _add_global_attributes(mom6_grid, mosaic_ds)

def write_MOM6_land_mask_file(mom6_grid):
    num_rows, num_cols = mom6_grid['supergrid']['area'].shape
    with netCDF4.Dataset(mom6_grid['filenames']['land_mask'], 'w', format='NETCDF3_CLASSIC') as land_mask_ds:
        # Dimenisons (of grid cells, not supergrid)
        land_mask_ds.createDimension('nx', num_cols / 2)
        land_mask_ds.createDimension('ny', num_rows / 2)

        # Variables & Values
        hmask = land_mask_ds.createVariable('mask', 'd', ('ny', 'nx'))
        hmask.standard_name = 'land fraction at T-cell centers'
        hmask.units = 'none'
        hmask[:] = numpy.logical_not(mom6_grid['mask'])

        # Global attributes
        _add_global_attributes(mom6_grid, land_mask_ds)

def write_MOM6_ocean_mask_file(mom6_grid):
    num_rows, num_cols = mom6_grid['supergrid']['area'].shape
    with netCDF4.Dataset(mom6_grid['filenames']['ocean_mask'], 'w', format='NETCDF3_CLASSIC') as ocean_mask_ds:
        # Dimenisons (of grid cells, not supergrid)
        ocean_mask_ds.createDimension('nx', num_cols / 2)
        ocean_mask_ds.createDimension('ny', num_rows / 2)

        # Variables & Values
        hmask = ocean_mask_ds.createVariable('mask', 'd', ('ny', 'nx'))
        hmask.standard_name = 'ocean fraction at T-cell centers'
        hmask.units = 'none'
        hmask[:] = mom6_grid['mask']

        # Global attributes
        _add_global_attributes(mom6_grid, ocean_mask_ds)

def main(argv):
    """Take a ROMS grid file and output a set of files to represent the MOM6 grid."""

    mom6_grid = setup_MOM6_grid(argv)

    if len(argv) == 2:
        roms_grid_filename = argv[1]
    else:
        print 'Usage: %s RGRID' % argv[0]
        print ''
        print 'Converts the ROMS horizontal grid stored in the NetCDF file RGRID into'
        print 'a collection of NetCDF files representing the MOM6 horizontal grid:'
        print ' - supergrid file ("%s")' % mom6_grid['filenames']['supergrid']
        print ' - topography file ("%s")' % mom6_grid['filenames']['topography']
        exit(1)

    roms_grid = read_ROMS_grid(roms_grid_filename)
    roms_grid = trim_ROMS_grid(roms_grid)
    mom6_grid = convert_ROMS_to_MOM6(mom6_grid, roms_grid)
    mom6_grid = approximate_MOM6_grid_metrics(mom6_grid)
    write_MOM6_supergrid_file(mom6_grid)
    write_MOM6_topography_file(mom6_grid)
    write_MOM6_solo_mosaic_file(mom6_grid)
    write_MOM6_land_mask_file(mom6_grid)
    write_MOM6_ocean_mask_file(mom6_grid)

if __name__ == "__main__":
    main(sys.argv)
