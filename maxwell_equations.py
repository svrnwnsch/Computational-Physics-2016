# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm

# Calculate numeric solution acc. to finite difference method, using backward derivative approx.

# Variables
#A = 1 # Amplitude
#omega = 20
#k     = 10
c = 0.7#omega / k
t_max = 1 #2 * np.pi / omega # Show one period
x_max = 1 #2 * np.pi / k     # Show one wavelength
sampling_rate_t = 4
sampling_rate_x = 4

# Generate mesh
x = np.linspace(0, x_max, sampling_rate_x)
t = np.linspace(0, t_max, sampling_rate_t)

# Set initial- and boundary-conditions
E = np.zeros([sampling_rate_t, sampling_rate_x])
B = np.zeros([sampling_rate_t, sampling_rate_x])
E[1, 0] = 1

for i in range(1, sampling_rate_x):
    for j in range(1, sampling_rate_t):
        E[i, j] = 1/(1 - 1/c**2) * (E[i-1, j] - E[i, j-1]/c**2 + B[i, j-1] - B[i-1, j])
print(E)
# Calculate exact solution
#def exact_solution(x, t, omega, k):
#    return np.real(A * np.exp(1j*(k * x)) * np.exp(1j*( - omega * t)))

# Plotting Routine
T, X = np.meshgrid(t, x)
plt.figure(1)
plt.pcolormesh(T, X, E, cmap = cm.seismic)
plt.xlabel('t')
plt.ylabel('x')
plt.colorbar()
plt.show()

for i in range(1, sampling_rate_x):
    for j in range(1, sampling_rate_t):
        E[i, j] = 1/(1 - 1/c**2) * (E[i-1, j] - E[i, j-1]/c**2 + B[i, j-1] - B[i-1, j])
print(E)        
plt.figure(2)
plt.pcolormesh(T, X, E, cmap = cm.seismic)
plt.xlabel('t')
plt.ylabel('x')
plt.colorbar()
plt.show()
# Show projections on t- and x-axis as well
#y_t = exact_solution(0, t, omega, k)
#fig2 = plt.plot(t, y_t)
#y_x = exact_solution(x, 0, omega, k)
#fig2 = plt.plot(x, y_x)

