# 1D Wave eqn solver for zero boundary condition
# Space domain: [xmin, xmax]
# Time domain: [0, tmax]
# Space mesh size: dx
# Time grid size: dt
# Speed of wave: c
# Initial wave: f(x)
# Initial speed: g(x)
# Special notice:
#   1. To compile, simply cd to the dir of this file and type "python2.7 wave1d.py" in command line
#   2. Here the exact solution assumes initial speed g(x) = 0
#   3. Please set dx and dt so that the number of space and time grids are integer 
#   4. For 1d Maxwell, don't use SI units (c = 1e8).
#       Please rescale everything so that c ~ O(1)
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm

def f(x):
    return np.exp(-10*(x-4)*(x-4))

def g(x):
    return 0

def uexact(x,t,c):
    return 0.5*(f(x-c*t)+f(x+c*t))

# Set up parameters for wave eqn:
xmin = 0
xmax = 8
tmax = 1
dx = 1e-2
dt = 1e-2
c = 0.5
# Boundary condition:
# bc = 0: zero boundary condition
bc = 0 

# Check stability condition
a = c*dt/dx;
if a > 1:
    print('Unstable Solver! c*dt/dx > 1. Try finer time grid (smaller dt) or rougher space grid (larger dx). ')
    quit()


n = int((xmax - xmin)/dx) # Spacial index: 0, 2, .., n. Interial: 1, 2, .., n-1
m = int(tmax/dt)          # Tempal index: 0, 2, .., m

# Initial Condition.
# u(x,0) = f(x); u'(x,0) = g(x)
# Now, f(x) is Gaussian function; g = 0.
u = np.empty((n+1,m+1))
x = np.empty(n+1)


for i in range(0, n+1):
    x[i] = i * dx

if bc == 0:
    for j in range(0, m+1):
        u[0,j] = 0
        u[n,j] = 0

for i in range(1, n):
    u[i,0] = f(x[i])
    u[i,1] = dt*g(x[i]) + (1-a*a)*f(x[i]) + 0.5*a*a*(f(x[i-1]) + f(x[i+1]))

# Compute solution using explicit central difference (for both time and space)
for j in range(2, m+1):
    for i in range(1, n):
        u[i,j] = -u[i,j-2]+ 2*(1-a*a)*u[i,j-1] + a*a*(u[i+1,j-1]+u[i-1,j-1])
                



# Compute error (g = 0)
err = 0
for i in range(1, n):
    tmperr = (u[i,m] - uexact(x[i], tmax, c))
    err += tmperr*tmperr
err = np.sqrt(err)

print(err)


# Plotting Routine
t = np.linspace(0,tmax,m+1)
T, X = np.meshgrid(t, x)
plt.figure(1)
plt.pcolormesh(T, X, u, cmap = cm.seismic)
plt.xlabel('t')
plt.ylabel('x')
plt.colorbar()
plt.show()

for j in range(2, m+1):
    for i in range(1, n):
        u[i,j] = -u[i,j-2]+ 2*(1-a*a)*u[i,j-1] + a*a*(u[i+1,j-1]+u[i-1,j-1])
print(u)        
plt.figure(2)
plt.pcolormesh(T, X, u, cmap = cm.seismic)
plt.xlabel('t')
plt.ylabel('x')
plt.colorbar()
plt.show()