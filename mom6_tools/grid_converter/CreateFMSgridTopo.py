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

# Open ROMS grid file
roms_grid_filename = 'CCS_7k_0-360_fred_grd.nc'
roms_grid = read_ROMS_grid(roms_grid_filename)

nj,ni = roms_grid['rho']['lon'].shape
nj -= 2
ni -= 2

# Supergrid shape: Smallest useful super-grid has a multiplier of 2
snj = 2 * nj
sni = 2 * ni

# Declare shapes
mom6_grid = dict()
mom6_grid['supergrid'] = dict()
mom6_grid['supergrid']['lon']   = numpy.zeros((snj+1,sni+1))
mom6_grid['supergrid']['lat']   = numpy.zeros((snj+1,sni+1))
mom6_grid['supergrid']['area']  = numpy.zeros((snj,  sni  ))
mom6_grid['supergrid']['dx']    = numpy.zeros((snj+1,sni  ))
mom6_grid['supergrid']['dy']    = numpy.zeros((snj,  sni+1))
mom6_grid['supergrid']['angle'] = numpy.zeros((snj+1,sni+1))

# Copy in data from ROMS file
mom6_grid['supergrid']['lon'][ ::2, ::2] = roms_grid['psi']['lon'][ :  , :  ] # Cell corners
mom6_grid['supergrid']['lon'][1::2,1::2] = roms_grid['rho']['lon'][1:-1,1:-1] # Cell centers (drop outside row and column)
mom6_grid['supergrid']['lon'][1::2, ::2] = roms_grid['u'  ]['lon'][1:-1, :  ] # U-points (drop outside row)
mom6_grid['supergrid']['lon'][ ::2,1::2] = roms_grid['v'  ]['lon'][ :  ,1:-1] # V-points (drop outside column)
mom6_grid['supergrid']['lat'][ ::2, ::2] = roms_grid['psi']['lat'][ :  , :  ] # Cell corners
mom6_grid['supergrid']['lat'][1::2,1::2] = roms_grid['rho']['lat'][1:-1,1:-1] # Cell centers (drop outside row and column)
mom6_grid['supergrid']['lat'][1::2, ::2] = roms_grid['u'  ]['lat'][1:-1, :  ] # U-points (drop outside row)
mom6_grid['supergrid']['lat'][ ::2,1::2] = roms_grid['v'  ]['lat'][ :  ,1:-1] # V-points (drop outside column)

mask = roms_grid['rho']['mask'][1:-1,1:-1]
depth = roms_grid['rho']['h'][1:-1,1:-1] * mask
mom6_grid['depth'] = depth

# Approximate edge lengths as great arcs
R = 6370.e3 # Radius of sphere
lat = mom6_grid['supergrid']['lat']
lon = mom6_grid['supergrid']['lon']
mom6_grid['supergrid']['dx'][:,:] = R * Spherical.angle_through_center( (lat[ :,1:],lon[ :,1:]), (lat[:  ,:-1],lon[:  ,:-1]) )
mom6_grid['supergrid']['dy'][:,:] = R * Spherical.angle_through_center( (lat[1:, :],lon[1:, :]), (lat[:-1,:  ],lon[:-1,:  ]) )

# Approximate angles using centered differences in interior
mom6_grid['supergrid']['angle'][:,1:-1] = numpy.arctan( (lat[:,2:]-lat[:,:-2]) / ((lon[:,2:]-lon[:,:-2]) * numpy.cos(numpy.deg2rad(lat[:,1:-1]))) )
# Approximate angles using side differences on left/right edges
mom6_grid['supergrid']['angle'][:,   0] = numpy.arctan( (lat[:, 1]-lat[:,  0]) / ((lon[:, 1]-lon[:,  0]) * numpy.cos(numpy.deg2rad(lat[:,   0]))) )
mom6_grid['supergrid']['angle'][:,  -1] = numpy.arctan( (lat[:,-1]-lat[:, -2]) / ((lon[:,-1]-lon[:, -2]) * numpy.cos(numpy.deg2rad(lat[:,  -1]))) )

# Approximate cell areas as that of spherical polygon
mom6_grid['supergrid']['area'][:,:] = R * R * Spherical.quad_area(lat, lon)

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

write_MOM6_supergrid(mom6_grid)
write_MOM6_topography(mom6_grid)
