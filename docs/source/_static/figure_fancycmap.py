
# Example of fancy colormaps use 

import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
import numpy as np
import envtoolkit.colors
import matplotlib.cm

# defining datasets to draw
# as in the contourf pyplot example
delta = 0.1
x = np.arange(-3.0, 3.0, delta)
y = np.arange(-2.0, 2.0, delta)
X, Y = np.meshgrid(x, y)
Z1 = mlab.bivariate_normal(X, Y, 1.0, 1.0, 0.0, 0.0)
Z2 = mlab.bivariate_normal(X, Y, 1.5, 0.5, 1, 1)
Z = 10.0 * (Z2 - Z1)
lmin, lmax = -1.5, 2

plt.figure(figsize=(10, 6))

plt.subplot(2, 2, 1)
# creating a Fancy cmap from a svg file
obj = envtoolkit.colors.FancyCmap('njgs/njavgcntyq.svg', reverse=False)
cmap = obj.makecmap()  # generating the cmap
cs = plt.pcolor(x, y, Z, cmap=cmap)
plt.clim(lmin, lmax)
cb = plt.colorbar(cs)
ticks = obj.generate_cbticks(cb)
cb.set_ticks(ticks)

plt.subplot(2, 2, 2)
vect = ['gold', 'lightsalmon', 'Purple',
        'orange', 'firebrick', 'lightskyblue']
# creating a gradient-type cmap from a list of colors
obj = envtoolkit.colors.FancyCmap(vect, smooth=True, reverse=True)
cmap = obj.makecmap()
cs = plt.pcolor(x, y, Z, cmap=cmap)
plt.clim(lmin, lmax)
cb = plt.colorbar(cs)

plt.subplot(2, 2, 3)
vect = np.array([[10, 20, 30],
                 [39, 100, 69],
                 [10, 30, 100],
                 [10, 50, 70],
                 [100, 0, 0]])/100.
# generating a discrete cmap from a numpy array of RGB percentage
obj = envtoolkit.colors.FancyCmap(vect, smooth=False)
cmap = obj.makecmap()
cs = plt.pcolor(x, y, Z, cmap=cmap)
plt.clim(lmin, lmax)
cb = plt.colorbar(cs)
ticks = obj.generate_cbticks(cb)
cb.set_ticks(ticks)

plt.savefig('figure_fancycmap.png', bbox_inches='tight')
