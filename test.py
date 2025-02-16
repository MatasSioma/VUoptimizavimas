import matplotlib.pyplot as plot
import numpy as np

plot.figure(figsize=(9, 3))
plot.subplot(131)
x = np.arange(0, 2*np.pi, 0.01)
plot.plot(x, np.sin(x),)
plot.subplot(132)
x = np.arange(0, 2*np.pi, 0.01)
plot.plot(x, np.sin(x))
plot.subplot(133)
x = np.arange(0, 2*np.pi, 0.01)
plot.plot(x, np.sin(x))
# plot.plot(x, y)
plot.show()