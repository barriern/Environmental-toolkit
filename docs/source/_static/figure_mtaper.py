
# Example on the use of multi-taper spectral analysis

import numpy as np
import pylab as plt
import envtoolkit.spectral

# Loading SST data
filein = np.load('sst_taper.npz')
sst = filein['sst']

# defining the delta t and the localisation of
# the error bar
deltat = 30*24*60*60
ferror = 1e-2

fmin_slope = 1e-2
fmax_slope = 1e-1
yinter_slope = 1e-1

fig = plt.figure()
ax = plt.gca()

# plotting the power spectrum (variance type)
spectrum, freq, error = \
        envtoolkit.spectral.plot_spectra(sst, deltat, ferror,
                                      spec_type='variance',
                                      nbandw=3, color='k')

# ploting the slope, computed from freq[0] to 1
# offy = 2, we move the slope by two grid cells on the vertical
slope1 = envtoolkit.spectral.plot_slope(spectrum, freq, fmax=1,
                                     offy=2, color='FireBrick', linestyle='--')

# plotting the slope, computed from 1 to freq[-1]
# offy = 2, we move the slope by two grid cells on the vertical
slope2 = envtoolkit.spectral.plot_slope(spectrum, freq,
                                     fmin=1, offy=2,
                                     color='DarkOrange', linestyle='--')

# plotting reference slope (slope = -2)
envtoolkit.spectral.plot_ref_slope(fmin_slope, fmax_slope,
                                yinter_slope, slope=-2,
                                color='k', linestyle='--')

# plotting reference slope (slope = -3)
envtoolkit.spectral.plot_ref_slope(fmin_slope, fmax_slope,
                                yinter_slope, slope=-3,
                                color='k', linestyle='--')

# Defining the x and y labels
ax.set_ylabel("Power spectrum density ($K^2.cpy^{-1}$)")
ax.set_xlabel("Frequency($cpy$)")

# setting the ylim
plt.ylim(1e-6, 1e2)

plt.savefig('figure_mtaper.png', bbox_inches='tight')

print slope1, slope2

