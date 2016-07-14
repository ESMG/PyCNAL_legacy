import numpy as np
import pycnal

def mld_from_temp(temp,grd):

  z = grd.vgrid.z_r[0]

  sst = temp[-1]

  mld_temp = sst - 0.2 #threshold of 0.2 degree, following de Boyer Montegut et al., 2004

  mld, lon, lat = pycnal.tools.isoslice(z,temp,mld_temp, grd)

  return mld



