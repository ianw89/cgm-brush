from cgmbrush.cgmbrush import *
import cProfile
import io
import pstats
from memory_profiler import profile
import matplotlib.pyplot as plt

@profile
def setup():
    r = 16
    size = 1000 * r
    a = np.zeros((size,size))
    a[(np.arange(0,size,step=10), np.arange(0,size,step=10))] = 1

    b = np.ones((20*r,20*r))

    return a,b

a,b = setup()

result = my_convolve(a, b)

plt.imshow(result)