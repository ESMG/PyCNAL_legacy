# -*- coding: utf-8 -*-
# vim: set fileencoding=utf-8 :

import numpy
import netCDF4
import Spherical

# Open ROMS grid file
nc = netCDF4.Dataset('CCS_7k_0-360_fred_grd.nc')

nj,ni = nc.variables['lon_rho'].shape
nj -= 2
ni -= 2

# Supergrid shape: Smallest useful super-grid has a multiplier of 2
snj = 2 * nj
sni = 2 * ni

# Declare shapes
lon   = numpy.zeros((snj+1,sni+1))
lat   = numpy.zeros((snj+1,sni+1))
area  = numpy.zeros((snj,  sni  ))
dx    = numpy.zeros((snj+1,sni  ))
dy    = numpy.zeros((snj,  sni+1))
angle = numpy.zeros((snj+1,sni+1))

# Copy in data from ROMS file
lon[ ::2, ::2] = nc.variables['lon_psi'][ :  , :  ] # Cell corners
lon[1::2,1::2] = nc.variables['lon_rho'][1:-1,1:-1] # Cell centers (drop outside row and column)
lon[1::2, ::2] = nc.variables['lon_u'  ][1:-1, :  ] # U-points (drop outside row)
lon[ ::2,1::2] = nc.variables['lon_v'  ][ :  ,1:-1] # V-points (drop outside column)
lat[ ::2, ::2] = nc.variables['lat_psi'][ :  , :  ] # Cell corners
lat[1::2,1::2] = nc.variables['lat_rho'][1:-1,1:-1] # Cell centers (drop outside row and column)
lat[1::2, ::2] = nc.variables['lat_u'  ][1:-1, :  ] # U-points (drop outside row)
lat[ ::2,1::2] = nc.variables['lat_v'  ][ :  ,1:-1] # V-points (drop outside column)

# Approximate edge lengths as great arcs
R = 6370.e3 # Radius of sphere
dx[:,:] = R * Spherical.angle_through_center( (lat[ :,1:],lon[ :,1:]), (lat[:  ,:-1],lon[:  ,:-1]) )
dy[:,:] = R * Spherical.angle_through_center( (lat[1:, :],lon[1:, :]), (lat[:-1,:  ],lon[:-1,:  ]) )

# Approximate angles using centered differences in interior
angle[:,1:-1] = numpy.arctan( (lat[:,2:]-lat[:,:-2]) / ((lon[:,2:]-lon[:,:-2]) * numpy.cos(numpy.deg2rad(lat[:,1:-1]))) )
# Approximate angles using side differences on left/right edges
angle[:,   0] = numpy.arctan( (lat[:, 1]-lat[:,  0]) / ((lon[:, 1]-lon[:,  0]) * numpy.cos(numpy.deg2rad(lat[:,   0]))) )
angle[:,  -1] = numpy.arctan( (lat[:,-1]-lat[:, -2]) / ((lon[:,-1]-lon[:, -2]) * numpy.cos(numpy.deg2rad(lat[:,  -1]))) )

# Approximate cell areas as that of spherical polygon
area[:,:] = R * R * Spherical.quad_area(lat, lon)

# Create a mosaic file
rg = netCDF4.Dataset('ocean_hgrid.nc', 'w', format='NETCDF3_CLASSIC')

# Dimensions
rg.createDimension('nx' , sni  )
rg.createDimension('nxp', sni+1)
rg.createDimension('ny' , snj  )
rg.createDimension('nyp', snj+1)
rg.createDimension('string', 5)

# Variables
hx     = rg.createVariable('x'       , 'float32', ('nyp','nxp',)) ; hx.units     = 'degrees'
hy     = rg.createVariable('y'       , 'float32', ('nyp','nxp',)) ; hy.units     = 'degrees'
hdx    = rg.createVariable('dx'      , 'float32', ('nyp','nx' ,)) ; hdx.units    = 'meters'
hdy    = rg.createVariable('dy'      , 'float32', ('ny' ,'nxp',)) ; hdy.units    = 'meters'
harea  = rg.createVariable('area'    , 'float32', ('ny' ,'nx' ,)) ; harea.units  = 'meters^2'
hangle = rg.createVariable('angle_dx', 'float32', ('nyp','nxp',)) ; hangle.units = 'degrees'
htile  = rg.createVariable('tile'    ,       'c', ('string',))

# Values
hx[:] = lon
hy[:] = lat
hdx[:] = dx
hdy[:] = dy
harea[:] = area
hangle[:] = angle
htile[:5] = 'tile1'
rg.close()

# Create a topography file
rg = netCDF4.Dataset('ocean_topog.nc', 'w', format='NETCDF3_CLASSIC')

# Dimensions
rg.createDimension('nx', ni)
rg.createDimension('ny', nj)
rg.createDimension('ntiles', 1)

# Variables
hdepth = rg.createVariable('depth', 'float32', ('ny','nx',)) ; hdepth.units = 'm'

# Values
hmask = nc.variables['mask_rho'][1:-1,1:-1]
hdepth[:] = nc.variables['h'][1:-1,1:-1] * hmask
rg.close()
