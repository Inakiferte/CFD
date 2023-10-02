#=======================================
# Import modules
#=======================================

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import ListedColormap
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm

#=======================================
# Physical input data
#=======================================

L    = 10.0                                # Channel lenght.
H    = 10.0                                # Channel height.
rho  = 1.225                               # Air density at sea level in kg/m^3. Incompresible = constant density at all mesh points.
T_in = 298.0                               # Normal air temperature in K.
P_in = 1.013E05                            # Air pressure at sea level in N/m^2.
V_in = 0.15                                # Air velocity in m/s.
R    = 8.31                                # Ideal gas constant in J/mol*K

#=======================================
# Numerical input data
#=======================================

N         = 11                             # Control volumes in x direction. Odd number!
M         = 11                             # Control volumes in y direction. Odd number!
delta_X   = L / N                          # Control volume lenght in x direction.
delta_Y   = H / M                          # Control volume height in y direction.
eps       = 1.0E-06                        # Gauss-Seidel method convergence parameter.
psi_B     = 0.0                            # Stream function at the bottom of the channel.
psi_T     = V_in * H                       # Stream function at the top of the channel.
psi_s     = (psi_B + psi_T) / 2.0          # Stream function start value.
t_max     = 100000                         # G-S loop steps.

# Mesh generation

x_cv = np.linspace(0, L, N+1)              # Generate x points with the same spacing
y_cv = np.linspace(0, H, M+1)              # Generate y points with the same spacing
x_p  = np.zeros(N+2)                       # Vectors for the centered control volumes with N+2 elements
y_p  = np.zeros(M+2)                       # Vectors for the centered control volumes with M+2 elements

# Fill x_p and y_p

for i in range(1, N+1):
    x_p[i] = (x_cv[i] + x_cv[i-1]) / 2.0

for j in range(1, M+1):
    y_p[j] = (y_cv[j] + y_cv[j-1]) / 2.0

# Set boundary points at the ends of the domain
x_p[0]  = x_cv[0]                         # Left boundary
x_p[-1] = x_cv[-1]                        # Right boundary
y_p[0]  = y_cv[0]                         # Bottom boundary
y_p[-1] = y_cv[-1]                        # Top boundary

#=======================================
# Define matrixes
#=======================================

psi       = np.zeros((N+2,M+2))           # Stream function matrix
psi_start = np.zeros((N+2,M+2))           # Start stream function
psi_ax    = np.zeros((N,M))               # Auxiliary matrix for G-S solver
rho       = np.zeros((N+2,M+2))           # Density matrix
a_P       = np.zeros((N+2,M+2))           # Auxiliar matrix at point p
a_E       = np.zeros((N+2,M+2))           # Auxiliar matrix at point east
a_S       = np.zeros((N+2,M+2))           # Auxiliar matrix at point south
a_W       = np.zeros((N+2,M+2))           # Auxiliar matrix at point west
a_N       = np.zeros((N+2,M+2))           # Auxiliar matrix at point north
b_P       = np.zeros((N+2,M+2))           # Generation term
conver    = np.zeros((N+2,M+2))           # Matrix that will be used in G-S algorithm to observe the convergence

#=======================================
# Initialize stream function and
# density matrixes
#=======================================
rho_in = P_in / R * T_in
for i in range(N+2):
    for j in range(M+2):
        psi[i,j]       = psi_s
        psi_start[i,j] = psi_s
        rho[i,j]       = rho_in

#=======================================
# Evaluate internal nodes
#=======================================

rho_ref = rho_in                         # Set the reference value as the input (rho_ref/rho=1)
for i in range(1,N+1):
    for j in range(1,M+1):
        a_E[i,j] = (rho_ref / rho[i + 1,j]) * (delta_Y / np.abs(x_p[i] - x_p[i + 1]))
        a_W[i,j] = (rho_ref / rho[i - 1,j]) * (delta_Y / np.abs(x_p[i] - x_p[i - 1]))
        a_N[i,j] = (rho_ref / rho[i,j + 1]) * (delta_X / np.abs(y_p[j] - y_p[j + 1]))
        a_S[i,j] = (rho_ref / rho[i,j - 1]) * (delta_X / np.abs(y_p[j] - y_p[j - 1]))
        a_P[i,j] = a_E[i,j] + a_W[i,j] + a_N[i,j] + a_S[i,j]

for i in range(N+2):
    for j in range(M+2):
        psi[i,0]  = 0.0                  # Bottom nodes
        psi[i,-1] = V_in * H             # Top nodes
for i in range(N+2):
    for j in range(M+2):
        a_E[i,0]  = 0.0                  # Bottom nodes
        a_W[i,0]  = 0.0                  # Bottom nodes
        a_N[i,0]  = 0.0                  # Bottom nodes
        a_S[i,0]  = 0.0                  # Bottom nodes
        a_P[i,0]  = 1.0                  # Bottom nodes
        b_P[i,0]  = psi_B                # Bottom nodes
        a_E[i,-1] = 0.0                  # Top nodes
        a_W[i,-1] = 0.0                  # Top nodes
        a_N[i,-1] = 0.0                  # Top nodes
        a_S[i,-1] = 0.0                  # Top nodes
        a_P[i,-1] = 1.0                  # Top nodes
        b_P[i,-1] = psi_T                # Top nodes

#=======================================
# Inlet Flow
#=======================================

for j in range(1, M+1):
    psi[0,j] = V_in * y_p[j]
    b_P[0,j] = V_in * y_p[j]
    a_P[0,j] = 1.0
#=======================================
# Here we start the G-S algorithm
# loop since the outlet of psi must
# be also solved interatively, so that
# the stream function converge.
#=======================================


psi_ax = psi                                 # Set the initial value of the axuliary stream function
for t in range(t_max):
    #=======================================
    # Outlet Flow
    #=======================================
    for j in range(1, M+1):
        psi[-1,j] = psi[-2,j]
        a_W[-1,j] = 1.0
        a_P[-1,j] = 1.0

    #=======================================
    # Compute G-S for non boundary points
    #=======================================

    r = np.sum(psi_ax)                       # Sum the values of the psi auxiliary matrix
    for i in range(1, N+1):
        for j in range(1, M+1):
            psi[i, j] = (a_E[i, j] * psi_ax[i + 1, j] + a_W[i, j] * psi_ax[i - 1, j] + a_N[i, j] * psi_ax[i, j + 1] + a_S[i,j] * psi_ax[i, j - 1] + b_P[i, j]) / a_P[i, j]

    sum = np.sum(psi)                        # Sum the values of psi matrix
    if np.abs(sum-r) <= eps:                 # Watch if the |psi - psi_aux| < precision
        print("=======================")
        print("G-S algorithm converged")
        print("The requiered steps has been:" +  " " + str(t))
        print("=======================")
        break                                # If corveges just break the main G-S loop
    else:
        psi_ax = psi                         # If it does not converge, reset the valu of psi_aux and start again the loop


#=======================================
# Generate the plot
#=======================================
output_plot = "Vertical_Paralel_Flow.pdf"
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
plt.xlabel(r"$X$", fontsize=fontsize)
plt.ylabel(r"$Y$", fontsize=fontsize)

# Plot
plt.contour(x_p, y_p, psi, levels=60, colors='blue', linewidths=2.0)

# Save figure
plt.savefig(output_plot,bbox_inches='tight')

print("Plot has been done sucesfully!")
print(" ")
print("You should find it on:" + " " + str(output_plot))
print("==============================")

