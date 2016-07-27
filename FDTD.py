# -*- coding: utf-8 -*-
"""
Created on Thu Jul 21 13:29:58 2016

@author: domi
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
#import scipy.constants as c

# Source
def f(x):
    return np.exp(-10*(x-4)*(x-4))

# Physical constants

c = 0.1 # Speed of wave

# Space domain: 0..z_max
z_max = 8

# Time domain: 0..t_max
t_max = 100

# Cell size and time stepping
dt = 1E-3
dz = 1E-2

# Spacial index
n = int(z_max / dz)

# Time index
m = int(t_max / dt)

z = np.zeros(n + 1)
for i in range(0, n + 1):
    z[i] = i * dz

t = np.zeros(m + 1)
for j in range(0, m + 1):
    t[j] = j * dt


# Constants
cb = dt / dz
ce = cb * c**2 

# Initialize vectors
E = np.zeros([n+1, m+1])
B = np.zeros([n+1, m+1])

# Boundary condition for all time-nodes j in 0..m
for j in range(0, m + 1):
        E[0, j] = 0
        E[n, j] = 0
        
# Initial Condition
for i in range(1, n):
    E[i, 0] = f(z[i])

eps = 0   
# Generate E- and B-fields with FDTD-Method
for j in range(0, m):
    if(j%100==0):
        print(float(j)/m)
    # E field loop    
    for i in range(1, n):
       
        E[i, j+1] = E[i, j] + ce * (B[i-1, j] - B[i, j])  
        if abs(E[i, j+1]) < eps:
            E[i, j + 1] = 0
    # B field loop
    for i in range(0, n):
        B[i, j+1] = B[i, j] + cb * (E[i, j] - E[i+1, j])    
'''
# Generate E- and B-fields with FDTD-Method (simplified)
for j in range(0, m):
    # E field loop
    for i in range(0, n-1):
        E[i+1, j+1] = E[i, j] + ce * B[i-1, j] - (1/cb - ce) * B[i, j]  - B[i, j+1] / cb
'''
#print(E)
#print(B)
# Plotting routine    
# ----------------

# Animation routine

def getE(time):
    line.set_ydata(E[:, time])
    time_text.set_text('t = {:5.2f}'.format(time * dt))
    return line

fig = plt.figure('Animation for {} s'.format(t_max))

ax = fig.add_subplot(111)
ax.set_ylim(-1, 1)
ax.set_xlabel('z')
ax.set_ylabel('E_x')

line, = ax.plot(z, E[:, 0])
time_text = ax.text(0.02, 0.95, '', transform=ax.transAxes)

ani = animation.FuncAnimation(fig, getE, range(0, E.shape[1], int(E.shape[1]/1000)),
                              interval=10, blit=False)

## Method to save animation in movie file needs ffmpg installed

#ani.save('Gaus_dz{}dt{}z{}t{}.mp4'.format(dt, dz, z_max, t_max), fps=30, extra_args=['-vcodec', 'libx264'])

# Produce snapshot at time t_hat

for t_hat in range(0, 5001, 1000):
    plt.figure('Snapshot at time t = {} s'.format(t_hat*dt))
    plt.plot(z, E[:, t_hat])
    plt.xlabel('z / m')
    plt.ylabel('E_x[t_hat = {}]'.format(t_hat*dt))
    plt.title('Snapshot at time t_hat = {} s'.format(t_hat*dt))
    plt.show()


'''
# Produce time evaluation at location z_hat
z_hat = 5
i_hat = int( z_hat / dz)  
print(i_hat)
plt.figure('Time evolution at location z_hat = {}'.format(z_hat)) 
plt.plot(t, E[i_hat, :])
plt.xlabel('t / s')
plt.ylabel('E_x[z_hat = {}]'.format(z_hat))
plt.title('Time evolution at location z_hat = {}'.format(z_hat))
plt.show()
'''

# Source    
    #e[z_source] = np.exp(-(t - 8)**2/4**2) #np.sin(2 * np.pi * 278E6 * t  * dt)      
    