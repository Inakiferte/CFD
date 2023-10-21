#=======================================
# Code for the Diagonal Flow EDS
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
# Physical input data
#=======================================

L        = 1.0                           # Channel lenght.
H        = L                             # Channel height.
Z        = 1.0                            # Channel depth.
V_in     = 5.15                           # Air velocity in m/s.
T_in     = 298.0                          # Normal air temperature in K.
P_in     = 1.013E05                       # Air pressure at sea level in N/m^2.
R        = 287                            # Ideal gas constant in J/mol*K
rho_in   = P_in / (R * T_in)              # Air density at sea level in kg/m^3. Incompresible = constant density at all mesh points.
delta_t  = 0.1
Pe       = 10E7
phi_high = 1.0
phi_low  = 0.0
phi_in   = 10.0
rhoGamma = 10
Gamma    = (rho_in * V_in * L) / Pe


#=======================================
# Numerical input data
#=======================================

N         = 140                             # Control volumes in x direction.
M         = N                              # Control volumes in y direction.
delta_X   = L / N                          # Control volume lenght in x direction.
delta_Y   = H / M                          # Control volume height in y direction.
eps       = 1.0E-08                        # Gauss-Seidel method convergence parameter.
t_max     = 1000000                        # G-S loop steps.
time_max  = 1000000                        # Time loop steps.
delta_Z   = 1.0                            # Just one step into z direction


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

phi       = np.zeros((M+2,N+2))           # Phi function
phi_zero  = np.zeros((M+2,N+2))           # Phi at t=0 function
phi_ax    = np.zeros((M+2,N+2))           # Phi auxiliar for G-S
rho       = np.zeros((M+2,N+2))           # Density matrix
a_P       = np.zeros((M+2,N+2))           # Auxiliar matrix at point p
a_E       = np.zeros((M+2,N+2))           # Auxiliar matrix at point east
a_S       = np.zeros((M+2,N+2))           # Auxiliar matrix at point south
a_W       = np.zeros((M+2,N+2))           # Auxiliar matrix at point west
a_N       = np.zeros((M+2,N+2))           # Auxiliar matrix at point north
b_P       = np.zeros((M+2,N+2))           # Generation term
v_xP      = np.zeros((M+2,N+2))           # Velocity in x direction
v_yP      = np.zeros((M+2,N+2))           # Velocity in y direction
m_e       = np.zeros((M+2,N+2))           # east mass flow rate
m_s       = np.zeros((M+2,N+2))           # south mass flow rate
m_w       = np.zeros((M+2,N+2))           # west mass flow rate
m_n       = np.zeros((M+2,N+2))           # north mass flow rate
D_e       = np.zeros((M+2,N+2))
D_w       = np.zeros((M+2,N+2))
D_n       = np.zeros((M+2,N+2))
D_s       = np.zeros((M+2,N+2))

#=======================================
# Initialize velocities
#=======================================

for i in range(N+2):
    for j in range(M+2):
        v_xP[j,i] = V_in * np.cos(np.pi / 4.0)
        v_yP[j,i] = V_in * np.sin(np.pi / 4.0)


#=======================================
# Initialize the map of \phi_{0}
#=======================================

# Internal nodes
for i in range(1, N+1):
    for j in range(1, M+1):
        phi_zero[j,i] = phi_in
        phi[j,i] = phi_in
# BC
for i in range(N+2):
    for j in range(M+2):
        phi_zero[j,0]  = phi_high  # Left
        phi_zero[j,-1] = phi_low   # Right
        phi_zero[0,i]  = phi_low   # Bottom
        phi_zero[-1,i] = phi_high  # Top
        phi[j,0]  = phi_high  # Left
        phi[j,-1] = phi_low   # Right
        phi[0,i]  = phi_low   # Bottom
        phi[-1,i] = phi_high  # Top


#=======================================
# Compute mass flow rates
#=======================================
for i in range(1,N+1):
    for j in range(1,M+1):
        m_e[j,i] = rho_in * (v_xP[j,i + 1] + v_xP[j,i]) * delta_Y * delta_Z / 2.0
        m_w[j,i] = rho_in * (v_xP[j,i - 1] + v_xP[j,i]) * delta_Y * delta_Z / 2.0
        m_n[j,i] = rho_in * (v_yP[j + 1,i] + v_yP[j,i]) * delta_X * delta_Z / 2.0
        m_s[j,i] = rho_in * (v_yP[j - 1,i] + v_yP[j,i]) * delta_X * delta_Z / 2.0

#=======================================
# Compute Di
#=======================================

# Internal nodes
for i in range(1, N+1):
    for j in range(1, M+1):
        D_e[j,i] = Gamma * delta_Y * delta_Z / np.abs(x_p[i] - x_p[i + 1])
        D_w[j,i] = Gamma * delta_Y * delta_Z / np.abs(x_p[i] - x_p[i - 1])
        D_n[j,i] = Gamma * delta_X * delta_Z / np.abs(y_p[j] - y_p[j + 1])
        D_s[j,i] = Gamma * delta_X * delta_Z / np.abs(y_p[j] - y_p[j - 1])


#=======================================
# Compute auxiliar matrixes
#=======================================

# Internal nodes
for i in range(1, N+1):
    for j in range(1,M+1):
        Pe_e = m_e[j,i] / D_e[j,i]
        Pe_w = m_w[j,i] / D_w[j,i]
        Pe_n = m_n[j,i] / D_n[j,i]
        Pe_s = m_s[j,i] / D_s[j,i]
        a_E[j,i] = D_e[j,i] * (np.abs(Pe_e) / (np.exp(np.abs(Pe_e)) - 1.0)) - ((m_e[j,i] - np.abs(m_e[j,i])) / 2.0)
        a_W[j,i] = D_w[j,i] * np.abs(Pe_w) / (np.exp(np.abs(Pe_w)) - 1.0) + ((m_w[j,i] + np.abs(m_w[j,i])) / 2.0)
        a_N[j,i] = D_n[j,i] * np.abs(Pe_n) / (np.exp(np.abs(Pe_n)) - 1.0) - ((m_n[j,i] - np.abs(m_n[j,i])) / 2.0)
        a_S[j,i] = D_s[j,i] * np.abs(Pe_s) / (np.exp(np.abs(Pe_s)) - 1.0) + ((m_s[j,i] + np.abs(m_s[j,i])) / 2.0)
        a_P[j,i] = a_E[j,i] + a_W[j,i] + a_N[j,i] + a_S[j,i] + ((rho_in * delta_X * delta_Y * delta_Z) / delta_t)
# BC
for i in range(N+2):
    for j in range(M+2):
        a_E[j,0]  = 0.0            # Left
        a_W[j,0]  = 0.0            # Left
        a_N[j,0]  = 0.0            # Left
        a_S[j,0]  = 0.0            # Left
        b_P[j,0]  = phi_high       # Left
        a_P[j,0]  = 1.0            # Left
        a_E[j,-1] = 0              # Right
        a_W[j,-1] = 0              # Right
        a_N[j,-1] = 0              # Right
        a_S[j,-1] = 0              # Right
        a_P[j,-1] = 1              # Right
        b_P[j,-1] = phi_low        # Right
        a_E[0,i]  = 0.0            # Bottom nodes
        a_W[0,i]  = 0.0            # Bottom nodes
        a_N[0,i]  = 0.0            # Bottom nodes
        a_S[0,i]  = 0.0            # Bottom nodes
        a_P[0,i]  = 1.0            # Bottom nodes
        b_P[0,i]  = phi_low        # Bottom nodes
        a_E[-1,i] = 0.0            # Top nodes
        a_W[-1,i] = 0.0            # Top nodes
        a_N[-1,i] = 0.0            # Top nodes
        a_S[-1,i] = 0.0            # Top nodes
        a_P[-1,i] = 1.0            # Top nodes
        b_P[-1,i] = phi_high       # Top nodes

print("=========================")
print("Starting loops")
print("=========================")
#=======================================
# Begin time step (t+delta_t)
#=======================================
for time in range(time_max):
    r2 = np.sum(phi_zero)                       # Sum the values for t = 0
    #=======================================
    # Begin Gauss-Seidel
    #=======================================
    phi_ax = phi_zero                             # Set the initial value of the axuliary stream function
    for t in range(t_max):
        # b_P internal nodes
        for i in range(1, N+1):
            for j in range(1, M+1):
                 b_P[j,i] = ((rho_in * delta_X * delta_Y * delta_Z) / delta_t) * phi_ax[j,i]

        r = np.sum(phi_ax)                        # Sum the values of the psi auxiliary matrix
        for i in range(1, N+1):
            for j in range(1, M+1):
                phi[j, i] = (a_E[j, i] * phi_ax[j, i + 1] + a_W[j, i] * phi_ax[j, i - 1] + a_N[j, i] * phi_ax[j + 1, i] + a_S[j,i] * phi_ax[j - 1, i] + b_P[j, i]) / a_P[j, i]

        sum = np.sum(phi)                        # Sum the values of psi matrix
        if np.abs(sum-r) <= eps:                 # Watch if the |psi - psi_aux| < precision
            print("=======================")
            print("G-S algorithm converged")
            print("The requiered steps has been:" +  " " + str(t))
            print("=======================")
            break                                # If corveges just break the main G-S loop
        else:
            phi_ax = phi
    sum2 = np.sum(phi)                           # Sum the values for i+1
    if np.abs(sum2 - r2)<=eps:
        print("=======================")
        print("Time loop converged")
        print("The requiered steps has been:" +  " " + str(time))
        print("=======================")
        break                                # If corveges just break the main G-S loop
    else:
        phi_zero = phi


#=======================================
# Generate the plot
#=======================================
output_plot1 = f"Results/Diagonal/Phi_{N}x{M}_{Pe}-EDS.pdf"
output_plot2 = f"Results/Diagonal/Vel_{N}x{M}_{Pe}-EDS.pdf"
output_plot3 = f"Results/Diagonal/Phi3D_{N}x{M}_{Pe}-EDS.pdf"
print(" ")
print("=================")
print("Starting the plot")
print(" ")
print("It might take a few seconds")
print("===========================")

# Figure specifications
fontsize=15

# Start first plot
plt.figure(1)
plt.figure(figsize = (8,8))
plt.tick_params(axis='both', which='both',length=3, width=1.0,
labelsize=15, right=True, top=True, direction='in') # For ticks in borders

# Figure labels
plt.xlabel(r"$X(m)$", fontsize=fontsize)
plt.ylabel(r"$Y(m)$", fontsize=fontsize)

# Plot
plt.contour(x_p, y_p, phi, levels=20, colors='r', linewidths=0.5)
plt.imshow(phi, cmap= 'cool', extent=(x_p.min(), x_p.max(), y_p.min(), y_p.max()), origin='lower')
plt.colorbar()

# Get the current axis
ax = plt.gca()

# Add text next to the color bar
text_x = 1.18  # Adjust the x-coordinate as needed
text_y = 0.5   # Adjust the y-coordinate as needed
text = r"$\phi$"
ax.text(text_x, text_y, text, transform=ax.transAxes, rotation=270, va='center', fontsize=fontsize)

# Save figure
plt.savefig(output_plot1,bbox_inches='tight')

# Start second plot
plt.figure(2)
plt.figure(figsize = (8,8))
plt.tick_params(axis='both', which='both',length=3, width=1.0,
labelsize=15, right=True, top=True, direction='in') # For ticks in borders

# Figure labels
plt.xlabel(r"$X(m)$", fontsize=fontsize)
plt.ylabel(r"$Y(m)$", fontsize=fontsize)

# Plot
P,Z = np.meshgrid(x_p, y_p)
plt.quiver(P[::2, ::2], Z[::2, ::2], v_xP[::2, ::2], v_yP[::2, ::2], color="green")


# Save figure
plt.savefig(output_plot2,bbox_inches='tight')

# Start third plot
plt.figure(3)
plt.figure(figsize = (8,8))
plt.tick_params(axis='both', which='both',length=3, width=1.0,
labelsize=15, right=True, top=True, direction='in') # For ticks in borders

# Figure labels
plt.xlabel(r"$X(m)$", fontsize=fontsize)
plt.ylabel(r"$Y(m)$", fontsize=fontsize)

# Plot
# Create a 3D figure
fig = plt.figure(figsize=(12, 12))
ax = fig.add_subplot(111, projection='3d')

# Create the surface plot
surf = ax.plot_surface(P, Z, phi, cmap='cool')

# Add a colorbar
cbar = fig.colorbar(surf)

# Save figure
plt.savefig(output_plot3,bbox_inches='tight')

plt.close(1)
plt.close(2)
plt.close(3)
