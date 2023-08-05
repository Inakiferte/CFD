import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import ListedColormap
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm


# Path to output values
file_output = "../code/fort.4"

# Load the data
data = np.loadtxt(file_output)

# Extract x, y, and p arrays from the data
X = data[:, 0]  # Mesh points in x direction
Y = data[:, 1]  # Mesh points in y direction
p = data[:, 2]  # Pressure
u = data[:, 3]  # Velocity in x direction
v = data[:, 4]  # Velocity in y direction

# Reshape X, Y, and p to create a 2D grid
nx = len(np.unique(X))
ny = len(np.unique(Y))
X = X.reshape((ny, nx))
Y = Y.reshape((ny, nx))
p = p.reshape((ny, nx))

# Reshape u and v to create a 2D grid with the same shape as X and Y
u = u.reshape((ny, nx))
v = v.reshape((ny, nx))

fontsize=15
plt.figure(figsize=(8, 8))
# plotting the pressure field as a contour
plt.contourf(X, Y, p, alpha=0.5, cmap=cm.viridis)
plt.colorbar()
# plotting the pressure field outlines
plt.contour(X, Y, p, cmap=cm.viridis)
# plotting velocity field
plt.quiver(X[::2, ::2], Y[::2, ::2], u[::2, ::2], v[::2, ::2],
           color="blue", scale=50)
plt.text(2.5, 1.0, 'Pressure', rotation=-90, ha='center', va='bottom',
         fontsize=fontsize)
plt.text(-0.25, 1.0, 'Y', rotation=90, ha='center', va='bottom',
         fontsize=fontsize)
plt.text(1.0, 1.90, 'Time Step = 700', ha='center', va='bottom',
         fontsize=fontsize, color='green')
plt.xlabel("X", fontsize=fontsize)
# plt.ylabel("Y (m)" , fontsize=fontsize)
plt.show()

