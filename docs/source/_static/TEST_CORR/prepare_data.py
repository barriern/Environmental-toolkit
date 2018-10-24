
from netCDF4 import Dataset
from pylab import *
import envtoolkit.ts
from netcdftime import utime
from datetime import datetime

fin = Dataset('data/psl_1d_19571201_20120331_natl.nc')
yyyymmdd = fin.variables['date'][:]
psl = fin.variables['psl'][:]

year = yyyymmdd/10000
month = (yyyymmdd-10000*year)/100
day = (yyyymmdd-10000*year-100*month)

date, pslmon = envtoolkit.ts.make_monthly_means(psl, yyyymmdd)

fin = Dataset('data/regime4_psl_1d_19571201_20100331_natl.nc', 'r')

mocc = fin.variables['month'][:]
yocc = fin.variables['year'][:]
lon = fin.variables['lon'][:]
lat = fin.variables['lat'][:]

mocc, yocc = np.meshgrid(mocc, yocc)
yocc[mocc==12] -=1

docc = yocc*100 + mocc
docc = np.ravel(docc)

dmin = np.max([docc.min(), date.min()])
dmax = np.min([docc.max(), date.max()])

idocc = np.nonzero((docc<=dmax) & (docc>=dmin))[0]
idate = np.nonzero((date<=dmax) & (date>=dmin))[0]

occu = fin.variables['occu_mth_corr'][:]
ncl, ny, nm = occu.shape
occu = np.reshape(occu, (ncl, ny*nm))

occu = occu[:, idocc].T
pslmon = pslmon[idate, :, :]

fout = Dataset('data_for_correlation.nc', 'w')
fout.createDimension('time', None)
fout.createDimension('ncl', 4)
fout.createDimension('lon', len(lon))
fout.createDimension('lat', len(lat))

fout.createVariable('occu', 'f', ('time', 'ncl'))
fout.createVariable('psl', 'f', ('time', 'lat', 'lon'))
fout.createVariable('lon', 'f', ('lon'))
fout.createVariable('lat', 'f', ('lat'))

fout.variables['occu'][:] = occu
fout.variables['lon'][:] = lon
fout.variables['lat'][:] = lat
fout.variables['psl'][:] = pslmon

fout.close()
