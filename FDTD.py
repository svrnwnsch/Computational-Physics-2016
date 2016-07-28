# -*- coding: utf-8 -*-
"""
Created on Thu Jul 21 13:29:58 2016

@authors: domi, severin
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

# Source
def f(x):
    return np.exp(-0.5*(x-5)**2/(1)**2)

# Alternative Source (has to be included between E- and B-field for loops)    
#e[z_source] = np.exp(-(t - 8)**2/4**2) #np.sin(2 * np.pi * 278E6 * t  * dt) 

# Speed of wave
c = .1

# Space domain: 0..z_max
z_max = 10

# Time domain: 0..t_max
t_max = 60

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
for j in range(0, m):
    E[0, j] = 0
    E[n, j] = 0
    #B[0, j] = 0
    #B[n, j] = 0

# Initial Condition
for i in range(1, n):
    E[i, 0] = f(z[i])

# Generate E- and B-fields with FDTD-Method
for j in range(0, m):
    if(j%100==0):
        print(float(j)/m)
    # E field loop    
    for i in range(1, n):
        E[i, j+1] = E[i, j] - ce * (B[i, j] - B[i-1, j])  
    # B field loop
    for i in range(0, n):
        B[i, j+1] = B[i, j] - cb * (E[i+1, j+1] - E[i, j+1])    
    

# Animation routine
# Show time evolution of E- and B-field in separate plots

# Create figure with 2 subplots
fig, (ax1, ax2) = plt.subplots(2, 1)

# Intialize two line objects (one for each axes)
line = [ax1.plot(z, E[:, 0], color='r'), ax2.plot(z, B[:, 0], color='b')] 
    
# Axes initialization
ax1.set_ylim(-1, 1)
ax1.set_xlim(0, z_max)
ax1.set_ylabel('E_x')
ax1.grid()
ax2.set_ylim(-10, 10)
ax2.set_xlim(0, z_max)
ax2.set_ylabel('B_y')
ax2.grid()

def getFields(time):
    line[0][0].set_ydata(E[:, time])    
    line[1][0].set_ydata(B[:, time])
    time_text.set_text('t = {:5.2f}'.format(time * dt))
    return line

#fig = plt.figure('Animation for {} s'.format(t_max))

#ax = fig.add_subplot(111)
#ax.set_ylim(-5, 5)
#ax.set_xlabel('z')
#ax.set_ylabel('B_y')

#line, = ax.plot(z, B[:, 0])
time_text = ax1.text(0.02, 0.95, '', transform=ax1.transAxes)

ani = animation.FuncAnimation(fig, getFields, range(0, B.shape[1], int(B.shape[1]/1000)),
                              interval=10, blit=False)


## Method to save animation in movie file needs ffmpg installed

#ani.save('Gaus_dz{}dt{}z{}t{}.mp4'.format(dt, dz, z_max, t_max), fps=30, extra_args=['-vcodec', 'libx264'])

'''
# Plotting routine
# Produce snapshot at time t_hat
for t_hat in range(0, int(m-1), int((m-1)/5)):
    plt.figure('Snapshot at time t = {} s'.format(t_hat*dt))
    plt.plot(z, E[:, t_hat])
    plt.plot(z, B[:, t_hat])    
    plt.xlabel('z / m')
    plt.ylabel('B_x[t_hat = {}]'.format(t_hat*dt))
    plt.title('Snapshot at time t_hat = {} s'.format(t_hat*dt))
    #plt.savefig('foo.pdf')
    plt.show()

# Produce time evolution at location z_hat
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

     
    