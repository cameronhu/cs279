import numpy as np
import matplotlib.pyplot as plt 
import csv
import os

data = []
with open(os.path.dirname(os.path.abspath(__file__))+'/scatter.csv','r') as f:
  points = csv.reader(f)
  for p in points:
    data.append([float(p[0]), float(p[1])])

data = np.array(data).transpose()

plt.axis([1.5,6.5,0,5])
plt.scatter(data[0,:],data[1,:])
plt.show()
