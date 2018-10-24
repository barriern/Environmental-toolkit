
# Example for discrete color mapping

import numpy as np
import pylab as plt
import matplotlib._cm
import envtoolkit.colors 

# function for adding text within the panels
def adding_text():
    """ Function for adding the text """
    for i in range(0, 3):
        for j in range(0, 3):
            plt.text(i, j, str(z[j, i]), bbox=bbox_props)

# Defining some plot parameters for the imshow function
plt.rcParams['image.cmap'] = 'jet'
plt.rcParams['image.interpolation'] = 'none'
plt.rcParams['image.origin'] = 'lower'

# defining the box settings for text display
bbox_props = dict(boxstyle="round,pad=0.3", fc="lightgray", ec="k", lw=1)

# creation of the data array
z = np.arange(0, 9)
z = np.reshape(z, (3,3))

# defining the number of colors and the boundaries
ncolors = 5
boundaries = np.array([0, 1, 3, 5, 7, 8])

# defining of a new jet colormap
cmap_jet_ncolors = envtoolkit.colors.subspan_default_cmap("jet", ncolors)

# defining the normalisation of the imshow function
norm = envtoolkit.colors.make_boundary_norm(ncolors, boundaries)

# Initialisation the figure 
plt.figure()
plt.subplots_adjust(left=0.05, right=0.98, wspace=0.1, hspace=0.2)

# Drawing the data with no change in cmap and norm (UNCORRECT) 
ax = plt.subplot(2, 2, 1)
cs = plt.imshow(z)
plt.title('No settings')
cs.set_clim(0, 9)
plt.colorbar(cs)
adding_text()

# Drawing the data with no change in norm (UNCORRECT)
ax = plt.subplot(2, 2, 2)
cs = plt.imshow(z, cmap=cmap_jet_ncolors)
plt.title('Colormap, no norm')
plt.colorbar(cs)
adding_text()

# Drawing the data with no change in cmap (UNCORRECT)
ax = plt.subplot(2, 2, 3)
cs = plt.imshow(z, norm=norm)
plt.title('Norm, no cmap')
plt.colorbar(cs)
adding_text()

# Drawing the data with changes in both cmap and norm (CORRECT)
ax = plt.subplot(2, 2, 3)
ax = plt.subplot(2, 2, 4)
cs = plt.imshow(z, cmap=cmap_jet_ncolors, norm=norm)
plt.title('Norm, Cmap')
plt.colorbar(cs)
adding_text()

plt.savefig("figure_col_norm.png", bbox_inches='tight')
