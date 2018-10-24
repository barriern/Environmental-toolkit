
from pylab import *

np.random.seed(1)
x = np.random.random(100)
y = np.random.random(100)

x = np.ma.masked_array(x)
y = np.ma.masked_array(y)

figure()
plot(x, linestyle="-", marker=".")
plot(y, linestyle="-", marker=".")
savefig('toto')

print "cov ", np.cov(x, y)[0, 1]
print "cov ma ", np.ma.cov(x, y)[0, 1]

x[10:20] = np.ma.masked
iok = np.nonzero(np.ma.getmaskarray(x)==False)[0]

print "cov iok ", np.cov(x[iok], y[iok])[0, 1]

print "cov ", np.cov(x, y)[0, 1]
print "ma.cov ", np.ma.cov(x, y)[0, 1]


figure()
plot(x, linestyle="-", marker=".")
plot(y, linestyle="-", marker=".")
savefig('toto2')


