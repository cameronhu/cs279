import sys
import numpy as np
import matplotlib.pyplot as plt
import csv

data = []
with open('scatter.csv','r') as f:
  points = csv.reader(f)
  for p in points:
    data.append([float(p[0]), float(p[1])])
data = np.array(data).transpose()

max_slope = 0
max_variance = 0

for i in np.arange(0, 5, 0.1):
    slope = i
    vec = np.array([0.0,1.0])
    if slope != np.inf:
      vec = np.array([1,slope])
    vec /= np.linalg.norm(vec)

    P = np.outer(vec,vec)
    print("Projection Matrix\n",P)

    slopeAligned = P.dot(data - data.mean(axis=1, keepdims=True)) + data.mean(axis=1, keepdims=True)

    phi = -np.arctan(slope)
    # clockwise rotation by phi
    R = np.array([[np.cos(phi), -np.sin(phi)],[np.sin(phi), np.cos(phi)]])
    xAligned = R.dot(slopeAligned)
    #compute variance here
    variance = np.var(xAligned[0,:])

    print("Empirical Variance is %.4f"%(variance))
    if variance > max_variance:
        max_variance = variance
        max_slope = slope

print(f'Max slope is {max_slope}')
print(f"Max variance is {max_variance}")

plt.figure("Points with PC line shown")
x_mean = np.mean(data[0,:])
y_mean = np.mean(data[1:,])
b = y_mean - x_mean * slope

line = np.array([
  [-20*vec[0],20*vec[0]],
  [-20*vec[1],20*vec[1]]
]) # Just extend the vector length 20 in + and - directions
line += data.mean(axis=1, keepdims=True) # Center the line on the data

plt.axis([1.5, 6.5,0,5])
plt.scatter(data[0,:],data[1,:])
plt.plot(line[0], line[1])        # abline

plt.figure("Projected points")
plt.scatter(slopeAligned[0,:],slopeAligned[1,:])
plt.axis([1.5, 6.5,0,5])
plt.show()
