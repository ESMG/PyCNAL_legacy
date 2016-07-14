import pycnal

# load the grid
srcgrd = pycnal.grid.get_ROMS_grid('PP_COSINE')
dstgrd = pycnal.grid.get_ROMS_grid('NWGOA3')

# make remap grid file for scrip
pycnal.remapping.make_remap_grid_file(srcgrd)
pycnal.remapping.make_remap_grid_file(srcgrd, Cpos='u')
pycnal.remapping.make_remap_grid_file(srcgrd, Cpos='v')
pycnal.remapping.make_remap_grid_file(dstgrd)
pycnal.remapping.make_remap_grid_file(dstgrd, Cpos='u')
pycnal.remapping.make_remap_grid_file(dstgrd, Cpos='v')

# compute remap weights
# input namelist variables for bilinear remapping at rho points
grid1_file = 'remap_grid_PP_COSINE_rho.nc'
grid2_file = 'remap_grid_NWGOA3_rho.nc'
interp_file1 = 'remap_weights_PP_COSINE_to_NWGOA3_bilinear_rho_to_rho.nc'
interp_file2 = 'remap_weights_NWGOA3_to_PP_COSINE_bilinear_rho_to_rho.nc'
map1_name = 'PP_COSINE to NWGOA3 Bilinear Mapping'
map2_name = 'NWGOA3 to PP_COSINE Bilinear Mapping'
num_maps = 1
map_method = 'bilinear'

pycnal.remapping.compute_remap_weights(grid1_file, grid2_file, \
              interp_file1, interp_file2, map1_name, \
              map2_name, num_maps, map_method)

# compute remap weights
# input namelist variables for bilinear remapping at u points
grid1_file = 'remap_grid_PP_COSINE_u.nc'
grid2_file = 'remap_grid_NWGOA3_u.nc'
interp_file1 = 'remap_weights_PP_COSINE_to_NWGOA3_bilinear_u_to_u.nc'
interp_file2 = 'remap_weights_NWGOA3_to_PP_COSINE_bilinear_u_to_u.nc'
map1_name = 'PP_COSINE to NWGOA3 Bilinear Mapping'
map2_name = 'NWGOA3 to PP_COSINE Bilinear Mapping'
num_maps = 1
map_method = 'bilinear'

pycnal.remapping.compute_remap_weights(grid1_file, grid2_file, \
              interp_file1, interp_file2, map1_name, \
              map2_name, num_maps, map_method)

# compute remap weights
# input namelist variables for bilinear remapping at v points
grid1_file = 'remap_grid_PP_COSINE_v.nc'
grid2_file = 'remap_grid_NWGOA3_v.nc'
interp_file1 = 'remap_weights_PP_COSINE_to_NWGOA3_bilinear_v_to_v.nc'
interp_file2 = 'remap_weights_NWGOA3_to_PP_COSINE_bilinear_v_to_v.nc'
map1_name = 'PP_COSINE to NWGOA3 Bilinear Mapping'
map2_name = 'NWGOA3 to PP_COSINE Bilinear Mapping'
num_maps = 1
map_method = 'bilinear'

pycnal.remapping.compute_remap_weights(grid1_file, grid2_file, \
              interp_file1, interp_file2, map1_name, \
              map2_name, num_maps, map_method)


# compute remap weights
# input namelist variables for bilinear remapping at u points
grid1_file = 'remap_grid_PP_COSINE_u.nc'
grid2_file = 'remap_grid_NWGOA3_rho.nc'
interp_file1 = 'remap_weights_PP_COSINE_to_NWGOA3_bilinear_u_to_rho.nc'
interp_file2 = 'remap_weights_NWGOA3_to_PP_COSINE_bilinear_rho_to_u.nc'
map1_name = 'PP_COSINE to NWGOA3 Bilinear Mapping'
map2_name = 'NWGOA3 to PP_COSINE Bilinear Mapping'
num_maps = 1
map_method = 'bilinear'

pycnal.remapping.compute_remap_weights(grid1_file, grid2_file, \
              interp_file1, interp_file2, map1_name, \
              map2_name, num_maps, map_method)


# compute remap weights
# input namelist variables for bilinear remapping at v points
grid1_file = 'remap_grid_PP_COSINE_v.nc'
grid2_file = 'remap_grid_NWGOA3_rho.nc'
interp_file1 = 'remap_weights_PP_COSINE_to_NWGOA3_bilinear_v_to_rho.nc'
interp_file2 = 'remap_weights_NWGOA3_to_PP_COSINE_bilinear_rho_to_v.nc'
map1_name = 'PP_COSINE to NWGOA3 Bilinear Mapping'
map2_name = 'NWGOA3 to PP_COSINE Bilinear Mapping'
num_maps = 1

pycnal.remapping.compute_remap_weights(grid1_file, grid2_file, \
              interp_file1, interp_file2, map1_name, \
              map2_name, num_maps, map_method)
