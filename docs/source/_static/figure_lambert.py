
# Example of using masked Lambert Projections

import pylab as plt
from envtoolkit.map import Lambert
import numpy as np

# defining plotting options
plt.rcParams['text.usetex'] = False
plt.rcParams['font.size'] = 9
plt.rcParams['lines.linewidth'] = 0.5

fig = plt.figure()
plt.subplots_adjust(hspace=0.0, wspace=0.15)

# defining nothern projection
lon1 = -80
lon2 = 30
lat1 = 20
lat2 = 80
lamb = Lambert(lon1, lon2, lat1, lat2)

ax = plt.subplot(2, 2, 1)
lamb.bmap.etopo()
lamb.bmap.drawcoastlines(color='k', linewidth=0.5)
lamb.bmap.drawparallels(np.arange(lat1, lat2+10, 10), labels=[1, 0, 0, 0])
lamb.bmap.drawmeridians(np.arange(lon1, lon2+20, 20), labels=[0, 0, 0, 1])

ax = plt.subplot(2, 2, 3)
lamb.bmap.etopo()
lamb.bmap.drawcoastlines(color='k', linewidth=0.5)
lamb.bmap.drawparallels(np.arange(lat1, lat2+10, 10), labels=[0, 0, 0, 0])
lamb.bmap.drawmeridians(np.arange(lon1, lon2+20, 20), labels=[0, 0, 0, 0])
lamb.make_mask()
lamb.add_lc_labels(spacelon=20, spacelat=10, fontsize=6)

# generating southern projection
lat1 = -80
lat2 = -20
lamb = Lambert(lon1, lon2, lat1, lat2)

ax = plt.subplot(2, 2, 2)
lamb.bmap.fillcontinents(color='lightgray')
lamb.bmap.drawcoastlines(color='k', linewidth=0.5)
lamb.bmap.drawparallels(np.arange(lat1, lat2+10, 10), labels=[0, 1, 0, 0])
lamb.bmap.drawmeridians(np.arange(lon1, lon2+20, 20), labels=[0, 0, 1, 0])

ax = plt.subplot(2, 2, 4)
lamb.bmap.fillcontinents(color='lightgray')
lamb.bmap.drawcoastlines(color='k', linewidth=0.5)
lamb.bmap.drawparallels(np.arange(lat1, lat2+10, 10), labels=[0, 0, 0, 0])
lamb.bmap.drawmeridians(np.arange(lon1, lon2+20, 20), labels=[0, 0, 0, 0])
lamb.make_mask()
lamb.add_lc_labels(spacelon=20, spacelat=10, fontsize=6)

plt.savefig('figure_lambert.png', bbox_inches='tight')
