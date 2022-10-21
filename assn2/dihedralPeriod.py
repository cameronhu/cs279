from pylab import *
import numpy as np

## TODO:
#  Play around with the weights defined in k
#  and the offsets defined in offset.  See
#  how the form of the potential energy, U, changes.

k = {1:1.0, 2:1.0, 3:1.0}
offset = {1:0.0, 2:pi, 3:0}

x = arange(0,2*pi,0.05)
sinusoidSum = np.array([k[n]*(1-cos(n*x-offset[n])) for n in (1,2,3)])
U = np.sum(sinusoidSum, 0)

plot(degrees(x),U)
show()





