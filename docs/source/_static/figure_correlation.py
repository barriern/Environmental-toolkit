from netCDF4 import Dataset
import envtoolkit.ts
import envtoolkit.map
import numpy as np
import pylab as plt

# load the dataset
fin = Dataset("combined_gridded_index_data.nc", "r")
amo = fin.variables['amo'][:]
nao = fin.variables['nao'][:]
sst = fin.variables["msl"][:]
lon = fin.variables["lon"][:]
lat = fin.variables["lat"][:]
fin.close()

# transfors the NAO/AMO index into 
# an array with dims. (nt, 2)
occu = np.vstack([nao, amo]).T


lag, corrpy = envtoolkit.ts.xcorr_ND(occu, sst, maxlag=5, use_covariance=False)


# generates the basemap bounding the lon/lat data domain
bmap = envtoolkit.map.make_bmap(lon, lat)

# meshgrid for contourf with basemap
lon, lat = np.meshgrid(lon, lat)

# correlation contours
levels = np.linspace(-1, 1, 11)

# initialisation of the figure object
plt.figure()
plt.subplots_adjust(wspace=0, hspace=0.06,
                    top=0.95, bottom=0.15, 
                    left=0.1, right=0.9)

# settings for text 
bbox = {'boxstyle':'round, pad=0.3', 'facecolor':'lightgray'}

# counter for subplot
cpt = 1

# loop over three lags we want to plot
for lagtest in [-1, 0, 1]:

    # recovering the index of the lag we want to plot
    ilag = np.nonzero(lag == lagtest)[0][0]

    # plotting the correlation with NAO index
    ax = plt.subplot(3, 2, cpt)
    if cpt==1:
        plt.title('NAO')
    cs1 = bmap.contourf(lon, lat, corrpy[0, :, :, ilag], levels=levels)
    bmap.drawcoastlines()
    cpt += 1
    plt.text(10, 70, "lag %d"%lagtest, bbox=bbox)

    # plotting the correlation with AMO index
    ax = plt.subplot(3, 2, cpt)
    if cpt==2:
        plt.title("AMO")
    cs = bmap.contourf(lon, lat, corrpy[1, :, :, ilag], levels=levels)
    bmap.drawcoastlines()
    cpt +=1 
    plt.text(10, 70, "lag %d"%lagtest, bbox=bbox)

# adding the colormap
cax = plt.axes([0.2, 0.1, 0.6, 0.03])
cb = plt.colorbar(cs, cax, orientation='horizontal')
cb.set_label('Cross-correlation')

plt.savefig("figure_correlation.png", bbox_inches='tight')
