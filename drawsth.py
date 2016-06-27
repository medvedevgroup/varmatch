import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

t = np.arange(0., 5., 0.2)
plt.plot(t,t,'r-o')
plt.plot(t,t/2, 'g-^')
plt.savefig('xx.png')
