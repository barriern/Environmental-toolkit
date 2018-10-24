
import numpy as np
import envtoolkit.ts
from pylab import *

x = np.random.random(200) - 0.5
print x.shape[0]

xstan = envtoolkit.ts.standardize(x)

figure()
subplots_adjust(top=0.95, bottom=0.05, left=0.07, right=0.95, hspace=0.07)
ax = subplot(2, 1, 1)
plot(x)
setp(ax.get_xticklabels(), visible=False)
grid(True)

ax = subplot(2, 1, 2)
plot(xstan)
grid(True)

savefig('figure_stan.png')

