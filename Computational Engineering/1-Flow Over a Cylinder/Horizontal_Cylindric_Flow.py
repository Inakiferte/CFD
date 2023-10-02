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
V_in = 1.15                                # Air velocity in m/s.
R    = 8.31                                # Ideal gas constant in J/mol*K


#=======================================
# Numerical input data
#=======================================

N         = 41                             # Control volumes in x direction. Odd number!
M         = 41                             # Control volumes in y direction. Odd number!
delta_X   = L / N                          # Control volume lenght in x direction.
delta_Y   = H / M                          # Control volume height in y direction.
eps       = 1.0E-08                        # Gauss-Seidel method convergence parameter.
psi_B     = 0.0                            # Stream function at the bottom of the channel.
psi_T     = V_in * H                       # Stream function at the top of the channel.
psi_s     = (psi_B + psi_T) / 2.0          # Stream function start value.
t_max     = 1000000                        # G-S loop steps.
psi_body  = psi_T / 2.0                    # Stream function inside the body. Add +- 1.5 for clockwise and anti-clockwise rotation.

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
# Define the cylinder
#=======================================

# Cylinder definition for the plot

D = 2.0                                  # Diameter
x0 = L / 2.0                             # X position of the center of the cylinder
y0 = H / 2.0                             # Y position of the center of the cylinder

# Cylinder definition

points = 100
X = np.linspace (-D / 2.0, D / 2.0, points)
Y = np.zeros(points)
for i in range(points):
    Y[ i ] = np.sqrt((D / 2.0)**2 - X[ i ] ** 2)


#=======================================
# Define matrixes
#=======================================

psi       = np.zeros((N+2,M+2))           # Stream function matrix
psi_start = np.zeros((N+2,M+2))           # Start stream function
psi_ax    = np.zeros((N+2,M+2))               # Auxiliary matrix for G-S solver
rho       = np.zeros((N+2,M+2))           # Density matrix
a_P       = np.zeros((N+2,M+2))           # Auxiliar matrix at point p
a_E       = np.zeros((N+2,M+2))           # Auxiliar matrix at point east
a_S       = np.zeros((N+2,M+2))           # Auxiliar matrix at point south
a_W       = np.zeros((N+2,M+2))           # Auxiliar matrix at point west
a_N       = np.zeros((N+2,M+2))           # Auxiliar matrix at point north
b_P       = np.zeros((N+2,M+2))           # Generation term
conver    = np.zeros((N+2,M+2))           # Matrix that will be used in G-S algorithm to observe the convergence
I_body    = np.zeros((N+2,M+2))           # Solid or fluid indentificator matrix

#=======================================
# Compute I_body matrix
#=======================================

for i in range(N+2):
    for j in range(M+2):
        x_aux = x_p[i]                                         # Compute the p cell x position
        y_aux = y_p[j]                                         # Compute the p cell y position
        distance = np.sqrt((x_aux - x0)**2 + (y_aux - y0)**2)  # Compute the p cell r position
        if distance <= (D / 2.0):  # Is the r vector inside the cylinder?
            I_body[i,j] = 1        # yes ==> I_body = 1
        else:
            I_body[i,j] = 0        # no ==> I_body = 0


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
        psi[0,j]  = 0.0                  # Bottom nodes
        psi[-1,j] = V_in * H             # Top nodes
for i in range(N+2):
    for j in range(M+2):
        a_E[0,j]  = 0.0                  # Bottom nodes
        a_W[0,j]  = 0.0                  # Bottom nodes
        a_N[0,j]  = 0.0                  # Bottom nodes
        a_S[0,j]  = 0.0                  # Bottom nodes
        a_P[0,j]  = 1.0                  # Bottom nodes
        b_P[0,j]  = psi_B                # Bottom nodes
        a_E[-1,j] = 0.0                  # Top nodes
        a_W[-1,j] = 0.0                  # Top nodes
        a_N[-1,j] = 0.0                  # Top nodes
        a_S[-1,j] = 0.0                  # Top nodes
        a_P[-1,j] = 1.0                  # Top nodes
        b_P[-1,j] = psi_T                # Top nodes

#=======================================
# Inlet Flow
#=======================================

for j in range(1, M+1):
    psi[j,0] = V_in * y_p[j]
    b_P[j,0] = V_in * y_p[j]
    a_P[j,0] = 1.0
#=======================================
# Here we start the G-S algorithm
# loop since the outlet of psi must
# be also solved interatively, so that
# the stream function converge.
#=======================================

print("========================")
print("Starting G-S algorithm")
print("======================")
psi_ax = psi                             # Set the initial value of the axuliary stream function
for t in range(t_max):
    #=======================================
    # Outlet Flow
    #=======================================
    for j in range(1, M+1):
        psi[j,-1] = psi[j,-2]
        a_W[j,-1] = 1.0
        a_P[j,-1] = 1.0
      
    #=======================================
    # Compute G-S for non boundary points
    #=======================================
    
    r = np.sum(psi_ax)                       # Sum the values of the psi auxiliary matrix
    for i in range(1, N+1):
        for j in range(1, M+1):
            if I_body[i,j] == 1:             # Is the p cell inside the body?
                psi[i,j] = psi_body          # yes, psi=psi_body
            else:                            # no, psi= G-S alsogirthm
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
# Compute the Velocity
#=======================================

# Define the requiered velocity matrix
v_xn = np.zeros((N+2,M+2))
v_xs = np.zeros((N+2,M+2))
v_ye = np.zeros((N+2,M+2))
v_yw = np.zeros((N+2,M+2))
v_xP = np.zeros((N+2,M+2))
v_yP = np.zeros((N+2,M+2))
v_P  = np.zeros((N+2,M+2))

for i in range(N+2):
    for j in range(M+2):
        v_xP[0,j]  = V_in             # Bottom nodes
        v_yP[0,j]  = 0.0              # Bottom nodes
        v_xP[-1,j] = V_in             # Top nodes
        v_yP[-1,j] = 0.0              # Top nodes

# Inlet flow
for j in range(1,M+1):
    v_xP[j,0] = V_in
    v_yP[j,0] = 0.0

# Fill the velocities
for i in range(1,N+1):
    for j in range(1,M+1):
        if I_body[i,j] == 1:
                v_xP[i,j] = 0.0
                v_yP[i,j] = 0.0
        else:
            v_xn[i,j] = (rho_ref / rho[i,j + 1]) * ((psi[i + 1,j] - psi[i,j]) / np.abs(y_p[j] - y_p[j + 1]))
            v_xs[i,j] = -(rho_ref / rho[i,j - 1]) * ((psi[i - 1,j] - psi[i,j]) / np.abs(y_p[j] - y_p[j - 1]))
            v_ye[i,j] = -(rho_ref / rho[i + 1,j]) * ((psi[i,j +1] - psi[i,j]) / np.abs(x_p[i] - x_p[i + 1]))
            v_yw[i,j] = (rho_ref / rho[i - 1,j]) * ((psi[i,j -1] - psi[i,j]) / np.abs(x_p[i] - x_p[i - 1]))
            v_xP[i,j] = (v_xn[i,j] + v_xs[i,j]) / 2.0
            v_yP[i,j] = (v_ye[i,j] + v_yw[i,j]) / 2.0
            v_P[i,j]  = np.sqrt(v_xP[i,j]**2 + v_yP[i,j]**2)

# Outlet flow
for j in range(1, M+1):
        v_xP[j,-1] = v_xP[j,-2]
        v_yP[j,-1] = v_yP[j,-2] 

#=======================================
# Compute the Temperature
#=======================================

# Input data
T_ref = T_in
V_ref = V_in
cp    = 1005 # J / Kg K

# Define temperature vector
T_P   = np.zeros((N+2, M+2))

# Compute temperature values
for i in range(N+2):
    for j in range(M+2):
        T_P[i,j] = T_ref + (V_ref**2 - v_P[i,j]**2) / (2.0 * cp)


#=======================================
# Compute the Pressure
#=======================================

# Physical values
P_ref = P_in
gamma = 1.4
P_P   = np.zeros((N+2, M+2))

for i in range(N+2):
    for j in range(M+2):
        P_P[i,j] = P_ref * (T_P[i,j] / T_ref) ** (gamma / (gamma - 1))


#=======================================
# Compute Lift and Dragg
#=======================================

L_plus = 0.0
L_minus = 0.0
D_plus = 0.0
D_minus=0.0

# Compute Lift
for i in range(1, N+1):
    for j in range(1, M+1):
        if I_body[i,j] == 1 and I_body[i + 1,j] == 0:
            L_plus = L_plus - P_P[i + 1,j] * delta_X 
            
        elif I_body[i,j] ==1 and I_body[i - 1, j] == 0:
            L_minus = L_minus + P_P[i - 1,j] * delta_X


# Compute Lift
for i in range(1, N+1):
    for j in range(1, M+1):
        if I_body[i,j] == 1 and I_body[i,j + 1] == 0:
            D_plus = D_plus - P_P[i,j + 1] * delta_Y 
        elif I_body[i,j] ==1 and I_body[i, j - 1] == 0:
            D_minus = D_minus + P_P[i,j - 1] * delta_Y

L = np.abs(L_minus + L_plus)
D = np.abs(D_minus + D_plus)

# Write the results forces
forces_file = f"Results/Static/ForcesS_{N}x{M}.txt"
with open(forces_file, 'w') as file:
    file.write("=================\n")
    file.write(f' The total lift over the cylinder is: L = {L} N\n')
    file.write(" \n")
    file.write("=================\n")
    file.write(f' The total dragg over the cylinder is: D = {D} N\n')

#=======================================
# Compute The Circulation
#=======================================

circ = 0.0

# Compute Circulation
for i in range(1, N+1):
    for j in range(1, M+1):
        if I_body[i,j] == 1 and I_body[i + 1,j] == 0:
            circ = circ - v_xs[i + 1,j] * delta_X 
            
        elif I_body[i,j] ==1 and I_body[i - 1, j] == 0:
            circ = circ + v_xn[i - 1,j] * delta_X


# Compute Circulation
for i in range(1, N+1):
    for j in range(1, M+1):
        if I_body[i,j] == 1 and I_body[i,j + 1] == 0:
            circ = circ - v_ye[i,j + 1] * delta_Y 
        elif I_body[i,j] ==1 and I_body[i, j - 1] == 0:
            circ = circ + v_yw[i,j - 1] * delta_Y

circ_abs = np.abs(circ)

# Write the results forces
circulation_file = f"Results/Static/CirculationS_{N}x{M}.txt"
with open(circulation_file, 'w') as file:
    file.write("=================\n")
    file.write(f' The circulation over the cylinder is {circ_abs} m^2 / s')

#=======================================
# Generate the plot
#=======================================
output_plot1 = f"Results/Static/Static_SL_{N}x{M}.pdf"
output_plot2 = f"Results/Static/Static_Vel_{N}x{M}.pdf"
output_plot3 = f"Results/Static/Static_Temp_{N}x{M}.pdf"
output_plot4 = f"Results/Static/Static_Pres_{N}x{M}.pdf"
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
plt.contour(x_p, y_p, psi, levels=60, colors='blue', linewidths=1.0)

# Cylinder plot
plt.plot(X + x0, Y + y0, color = "black")
plt.plot(X + x0,-Y + y0, color = "black")
plt.fill_between(X + x0, Y + y0, -Y + y0, color="black")

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
P,Z = np.meshgrid(x_p, y_p)      # We generate the mesh points (for plotting, not computing)
plt.quiver(P[::2, ::2], Z[::2, ::2], v_xP[::2, ::2], v_yP[::2, ::2], color="green") 

# Cylinder plot
plt.plot(X + x0, Y + y0, color = "black")
plt.plot(X + x0,-Y + y0, color = "black")
plt.fill_between(X + x0, Y + y0, -Y + y0, color="black")

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
plt.imshow(T_P, cmap= 'inferno', extent=(x_p.min(), x_p.max(), y_p.min(), y_p.max()), origin='lower')
plt.colorbar()

# Get the current axis
ax = plt.gca()

# Add text next to the color bar
text_x = 1.25  # Adjust the x-coordinate as needed
text_y = 0.5   # Adjust the y-coordinate as needed
text = r"$T~(K)$"
ax.text(text_x, text_y, text, transform=ax.transAxes, rotation=270, va='center', fontsize=fontsize)

# Cylinder plot
plt.plot(X + x0, Y + y0, color = "black")
plt.plot(X + x0,-Y + y0, color = "black")
plt.fill_between(X + x0, Y + y0, -Y + y0, color="black")

# Save figure
plt.savefig(output_plot3,bbox_inches='tight')

# Start fourth plot
plt.figure(4)
plt.figure(figsize = (8,8))
plt.tick_params(axis='both', which='both',length=3, width=1.0,
labelsize=15, right=True, top=True, direction='in') # For ticks in borders

# Figure labels
plt.xlabel(r"$X(m)$", fontsize=fontsize)
plt.ylabel(r"$Y(m)$", fontsize=fontsize)

# Plot
plt.imshow(P_P, cmap= 'viridis', extent=(x_p.min(), x_p.max(), y_p.min(), y_p.max()), origin='lower')
plt.colorbar()

# Get the current axis
ax = plt.gca()

# Add text next to the color bar
text_x = 1.25  # Adjust the x-coordinate as needed
text_y = 0.5   # Adjust the y-coordinate as needed
text = r"$P~(Pa)$"
ax.text(text_x, text_y, text, transform=ax.transAxes, rotation=270, va='center', fontsize=fontsize)

# Cylinder plot
plt.plot(X + x0, Y + y0, color = "black")
plt.plot(X + x0,-Y + y0, color = "black")
plt.fill_between(X + x0, Y + y0, -Y + y0, color="black")

# Save figure
plt.savefig(output_plot4,bbox_inches='tight')


print("Plots has been done sucesfully!")
print(" ")
print("You should find them on:" + " " + str(output_plot1) +  " " + "and" +  " " + str(output_plot2) + " " + "and" + " " + str(output_plot3) + " " + "and" + " " + str(output_plot4))
print("==============================")


plt.close(1)
plt.close(2)
plt.close(3)
plt.close(4)
