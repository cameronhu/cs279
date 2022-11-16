import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import norm
import statistics

# Plot between -10 and 10 with .001 steps.
x4 = np.arange(-4, 5, 1)
x5 = np.arange(-5, 6, 1)

# Calculating mean and standard deviation
mean4 = statistics.mean(x4)
sd4 = statistics.stdev(x4)
mean5 = statistics.mean(x5)
sd5 = statistics.mean(x5)

norm4 = norm.pdf(x4, mean4, sd4)
norm5 = norm.pdf(x5, mean5, sd5)
total = 0
for i in norm4:
    total += i

print(total)

plt.plot(x4, norm4)
plt.show()
