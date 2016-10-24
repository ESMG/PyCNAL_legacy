# -*- coding: utf-8 -*-
# vim: set fileencoding=utf-8 :

import numpy
import netCDF4
import Spherical

def read_ROMS_grid(roms_grid_filename):
    """Load the ROMS grid from a NetCDF file."""

    subgrids = ['psi', 'rho', 'u', 'v'] # also available: 'vert'
    fields = ['lon', 'lat', 'mask'] # also available: 'x', 'y'

    roms_grid = dict()
    with netCDF4.Dataset(roms_grid_filename) as roms_ds:
        # extract fields named based on which grid they are on
        for subgrid in subgrids:
            roms_grid[subgrid] = dict()
            for field in fields:
                var_name = field + '_' + subgrid
                roms_grid[subgrid][field] = roms_ds.variables[var_name][:]

        # extract fields that don't follow the above naming pattern
        roms_grid['rho']['h'] = roms_ds.variables['h'][:] # on the rho grid, but not called "h_rho"

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

def convert_ROMS_to_MOM6(roms_grid):
    """Convert the ROMS grid data into a skeleton MOM6 grid, mainly by
    merging the four sets of point locations from the ROMS grid
    into a single supergrid for MOM6."""

    # Double the size of the *trimmed* ROMS grid to merge the four
    # sets of points.
    num_rows, num_cols = roms_grid['rho']['lon'].shape
    num_rows *= 2
    num_cols *= 2

    mom6_grid = dict()
    mom6_grid['supergrid'] = dict()

    # Store the size for later
    mom6_grid['num_rows'] = num_rows
    mom6_grid['num_cols'] = num_cols

    # Copy points from ROMS grid
    for field in ['lon', 'lat']:
        mom6_grid['supergrid'][field] = numpy.zeros((num_rows+1,num_cols+1))
        mom6_grid['supergrid'][field][ ::2, ::2] = roms_grid['psi'][field] # outer
        mom6_grid['supergrid'][field][1::2,1::2] = roms_grid['rho'][field] # inner
        mom6_grid['supergrid'][field][1::2, ::2] = roms_grid[ 'u' ][field] # between e/w
        mom6_grid['supergrid'][field][ ::2,1::2] = roms_grid[ 'v' ][field] # between n/s

    mom6_grid['depth'] = roms_grid['rho']['h'] * roms_grid['rho']['mask']

    return mom6_grid

def approximate_MOM6_grid_metrics(mom6_grid):
    """Fill in missing MOM6 supergrid metrics by computing best guess
    values."""

    num_rows = mom6_grid['num_rows']
    num_cols = mom6_grid['num_cols']
    lat = mom6_grid['supergrid']['lat']
    lon = mom6_grid['supergrid']['lon']

    # Declare shapes
    mom6_grid['supergrid']['area']  = numpy.zeros((num_rows,  num_cols  ))
    mom6_grid['supergrid']['dx']    = numpy.zeros((num_rows+1,num_cols  ))
    mom6_grid['supergrid']['dy']    = numpy.zeros((num_rows,  num_cols+1))
    mom6_grid['supergrid']['angle'] = numpy.zeros((num_rows+1,num_cols+1))

    # Approximate edge lengths as great arcs
    R = 6370.e3 # Radius of sphere
    mom6_grid['supergrid']['dx'][:,:] = R * Spherical.angle_through_center( (lat[ :,1:],lon[ :,1:]), (lat[:  ,:-1],lon[:  ,:-1]) )
    mom6_grid['supergrid']['dy'][:,:] = R * Spherical.angle_through_center( (lat[1:, :],lon[1:, :]), (lat[:-1,:  ],lon[:-1,:  ]) )

    # Approximate angles using centered differences in interior, and side differences on left/right edges
    cos_lat = numpy.cos(numpy.radians(lat))
    mom6_grid['supergrid']['angle'][:,1:-1] = numpy.arctan( (lat[:,2:] - lat[:,:-2]) / ((lon[:,2:] - lon[:,:-2]) * cos_lat[:,1:-1]) )
    mom6_grid['supergrid']['angle'][:, 0  ] = numpy.arctan( (lat[:, 1] - lat[:, 0 ]) / ((lon[:, 1] - lon[:, 0 ]) * cos_lat[:, 0  ]) )
    mom6_grid['supergrid']['angle'][:,-1  ] = numpy.arctan( (lat[:,-1] - lat[:,-2 ]) / ((lon[:,-1] - lon[:,-2 ]) * cos_lat[:,-1  ]) )

    # Approximate cell areas as that of spherical polygon
    mom6_grid['supergrid']['area'][:,:] = R * R * Spherical.quad_area(lat, lon)

    return mom6_grid

def write_MOM6_supergrid(mom6_grid):
    """Save the MOM6 supergrid data into its own file."""
    num_rows, num_cols = mom6_grid['supergrid']['area'].shape

    tile_str = 'tile1'
    string_len = len(tile_str)

    with netCDF4.Dataset('ocean_hgrid.nc', 'w', format='NETCDF3_CLASSIC') as hgrid_ds:
        # Dimensions
        hgrid_ds.createDimension('nx',  num_cols)
        hgrid_ds.createDimension('nxp', num_cols+1)
        hgrid_ds.createDimension('ny',  num_rows)
        hgrid_ds.createDimension('nyp', num_rows+1)
        hgrid_ds.createDimension('string', string_len)

        # Variables & Values
        hx = hgrid_ds.createVariable('x', 'float32', ('nyp','nxp',))
        hx.units = 'degrees'
        hx[:] = mom6_grid['supergrid']['lon']

        hy = hgrid_ds.createVariable('y', 'float32', ('nyp','nxp',))
        hy.units = 'degrees'
        hy[:] = mom6_grid['supergrid']['lat']

        hdx = hgrid_ds.createVariable('dx', 'float32', ('nyp','nx',))
        hdx.units = 'meters'
        hdx[:] = mom6_grid['supergrid']['dx']

        hdy = hgrid_ds.createVariable('dy', 'float32', ('ny','nxp',))
        hdy.units = 'meters'
        hdy[:] = mom6_grid['supergrid']['dy']

        harea = hgrid_ds.createVariable('area', 'float32', ('ny','nx',))
        harea.units = 'meters^2'
        harea[:] = mom6_grid['supergrid']['area']

        hangle = hgrid_ds.createVariable('angle_dx', 'float32', ('nyp','nxp',))
        hangle.units = 'degrees'
        hangle[:] = mom6_grid['supergrid']['angle']

        htile = hgrid_ds.createVariable('tile', 'c', ('string',))
        htile[:string_len] = tile_str

def write_MOM6_topography(mom6_grid):
    """Save the MOM6 ocean topography field in a separate file."""
    num_rows, num_cols = mom6_grid['depth'].shape

    with netCDF4.Dataset('ocean_topog.nc', 'w', format='NETCDF3_CLASSIC') as topog_ds:
        # Dimensions
        topog_ds.createDimension('nx', num_cols)
        topog_ds.createDimension('ny', num_rows)
        topog_ds.createDimension('ntiles', 1)

        # Variables & Values
        hdepth = topog_ds.createVariable('depth', 'float32', ('ny','nx',))
        hdepth.units = 'm'
        hdepth[:] = mom6_grid['depth']

# Expecting a ROMS grid file as input
roms_grid_filename = 'CCS_7k_0-360_fred_grd.nc'
roms_grid = read_ROMS_grid(roms_grid_filename)
roms_grid = trim_ROMS_grid(roms_grid)
mom6_grid = convert_ROMS_to_MOM6(roms_grid)
mom6_grid = approximate_MOM6_grid_metrics(mom6_grid)
write_MOM6_supergrid(mom6_grid)
write_MOM6_topography(mom6_grid)
