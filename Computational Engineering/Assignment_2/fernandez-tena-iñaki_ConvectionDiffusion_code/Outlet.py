#=======================================
# Code for the outlet values in S-H
# problem.
# Master in Space and Aeronautical
# Engineering.
# Computational Engineering: Assginament
# 2.
# Author: Iñaki Fernandez Tena
# email: inakiphy@gmail.com
# last update: 21/10/2023
#=======================================

#=======================================
# Import modules
#=======================================

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import ListedColormap
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm

#=======================================
# Reference values
#=======================================

x_plot_exp = [0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0]
phi_10     = [1.989, 1.402, 1.146, 0.946, 0.775, 0.621, 0.480, 0.349, 0.227, 0.111, 0.000]
phi_10E3   = [2.0000, 1.9990, 1.9997, 1.9850, 1.8410, 0.9510, 0.1540, 0.0010, 0.0000, 0.0000, 0.0000]
phi_10E6   = [2.000, 2.000, 2.000, 1.999, 1.964, 1.000, 0.036, 0.001, 0.000, 0.000, 0.000]

#=======================================
# Computed values
#=======================================

N   = 200
M   = 100
Pe1 = 10
Pe2 = 10E3
Pe3 = 10E6

# Path to input file
file_input1 = f"Results/S-H/Outlet_{N}x{M}_{Pe1}-EDS.dat"
file_input2 = f"Results/S-H/Outlet_{N}x{M}_{Pe2}-EDS.dat"
file_input3 = f"Results/S-H/Outlet_{N}x{M}_{Pe3}-EDS.dat"


# Load data
x1   = np.loadtxt(file_input1)[:,0]
phi1 = np.loadtxt(file_input1)[:,1]
x2   = np.loadtxt(file_input2)[:,0]
phi2 = np.loadtxt(file_input2)[:,1]
x3   = np.loadtxt(file_input3)[:,0]
phi3 = np.loadtxt(file_input3)[:,1]


#=======================================
# Compute errors
#=======================================
r_ref10   = np.sum(phi_10) / float(len(phi_10))
r_ref10E3 = np.sum(phi_10E3) / float(len(phi_10E3))
r_ref10E6 = np.sum(phi_10E6) / float(len(phi_10E6))
r_phi1    = np.sum(phi1) / float(len(phi1))
r_phi2    = np.sum(phi2) / float(len(phi2))
r_phi3    = np.sum(phi3) / float(len(phi3))

error_10   = (np.abs(r_phi1 - r_ref10)) / r_ref10 * 100
error_10E3 = (np.abs(r_phi2 - r_ref10E3)) / r_ref10E3 * 100
error_10E6 = (np.abs(r_phi3 - r_ref10E6)) / r_ref10E6 * 100

print("=================================")
print("EDS:")
print("rho/Gamma=10 error:", error_10)
print("rho/Gamma=10E3 error:", error_10E3)
print("rho/Gamma=10E6 error:", error_10E6)
print("=================================")
#=======================================
# Generate the plot
#=======================================
output_plot = f"Results/S-H/OutletPlot_{N}x{M}-EDS.pdf"
print(" ")
print("=================")
print("Starting the plot")
print(" ")
print("It might take a few seconds")
print("===========================")

# Figure specifications
fontsize=15
plt.figure(figsize = (8,8))
plt.tick_params(axis='both', which='both',length=3, width=1.0,
labelsize=15, right=True, top=True, direction='in') # For ticks in borders

# Figure labels
plt.xlabel(r"$X(m)$", fontsize=fontsize)
plt.ylabel(r"$\phi$", fontsize=fontsize)


# Plot
plt.plot(x_plot_exp,phi_10,color="blue",linestyle = "dashed", label=r"Ref:$\rho/\Gamma=10$")
plt.plot(x1, phi1, color="blue", linestyle="solid", label=r"$\rho/\Gamma=10$")

plt.plot(x_plot_exp,phi_10E3,color="darkorange",linestyle = "dashed", label=r"Ref:$\rho/\Gamma=10E3$")
plt.plot(x2, phi2, color="darkorange", linestyle="solid", label=r"$\rho/\Gamma=10E3$")

plt.plot(x_plot_exp,phi_10E6,color="green",linestyle = "dashed", label=r"Ref:$\rho/\Gamma=10E6$")
plt.plot(x3, phi3, color="green", linestyle="solid", label=r"$\rho/\Gamma=10E6$")

# Legend specifications
plt.legend(fontsize=fontsize-2, loc='upper right')

# Save figure
plt.savefig(output_plot,bbox_inches='tight')

print("Plot has been done sucesfully!")
