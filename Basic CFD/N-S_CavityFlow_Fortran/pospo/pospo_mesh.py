import matplotlib.pyplot as plt
import numpy as np
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D

# Path to the mesh file
file_mesh = "../code/fort.1"

# Load fort.1 file

X = np.loadtxt(file_mesh)[:,0] # x points 
Y = np.loadtxt(file_mesh)[:,1]  # y points
Mx_line = np.loadtxt(file_mesh)[0,2] # Mesh points in x direction
My_line = np.loadtxt(file_mesh)[0,3] # Mesh points in y direction

Mx = int(Mx_line)
My = int(My_line)

# Create the plot
fontsize=15

plt.figure(figsize = (8,8))
plt.xlabel("X (m)", fontsize=fontsize)
plt.ylabel("Y (m)" , fontsize=fontsize)

# Plot the points
plt.scatter (X, Y, s=15, color = "darkorange")

# Plot horizontal dashed lines
for i, yi in enumerate(Y):
    linestyle = 'solid' if i == 0 else 'dashed'
    color = 'black' if i == 0 else 'darkgray'
    plt.plot(X, [yi] * len(X), linestyle=linestyle, color=color)

# Plot vertical dashed lines
for i, xi in enumerate(X):
    linestyle = 'solid' if i == 0 or i == Mx-1 else 'dashed'
    color = 'black' if i == 0 or i == Mx-1 else 'darkgray'
    plt.plot([xi] * len(Y), Y, linestyle=linestyle, color=color)

plt.savefig("CavityMeshF.pdf", bbox_inches='tight')

