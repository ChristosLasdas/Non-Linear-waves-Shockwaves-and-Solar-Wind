import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import scienceplots

plt.style.use(['science', 'notebook'])
plt.rcParams.update({"text.usetex" : True})
n = 2000 #number of iterations

nx = 1001 #x resolution
    
    
x_min = 0.
x_max = 10.
    
v = np.zeros(nx)
p = np.ones(nx)
rho = np.ones(nx)
mach = np.zeros(nx)

# initial conditions
gamma = 5/3
u0 = 1.5
v[0] = u0
rho[0] = 1.
p[0] = 1.
eta = 0.01

e = 0.5 * rho * v**2 + ( 1 / ( gamma - 1 ) ) * p    

drho_dt = np.zeros(nx)
dv_dt = np.zeros(nx)
dp_dt = np.zeros(nx)
de_dt = np.zeros(nx)


x = np.zeros(nx)

t = 0.

dx = ( x_max - x_min ) / ( nx - 1 )

dt = dx * 0.1


for i in range(0, nx):
    x[i] = i * dx



for k in range(0, n):
    
    for i in range (1, nx - 1):
        
        drho_dt[i] = - v[i] * ( ( rho[i] - rho[i - 1] ) / dx ) - rho[i] * ( ( v[i + 1] - v[i - 1] ) / ( 2 * dx ) )
        
        dv_dt[i] = - v[i] * ( ( v[i] - v[i - 1] ) / dx ) - ( 1 / rho[i] ) * ( ( p[i + 1] - p[i - 1] ) / ( 2 * dx ) )\
                   + eta * ( ( v[i + 1] - 2 * v[i] + v[i - 1] ) / dx**2 )      

        de_dt[i] =  - v[i] * ( ( e[i] - e[i - 1] ) / dx ) - ( e[i] + p[i] ) * ( ( v[i + 1] - v[i - 1] ) / ( 2 * dx ) )\
                    - v[i] * ( ( p[i + 1] - p[i - 1] ) / ( 2 * dx ) )


    for i in range (1, nx - 1):
        
        rho[i] = rho[i] + drho_dt[i] * dt
        v[i] = v[i] + dv_dt[i] * dt
        mach[i] = v[i] / np.sqrt( gamma * ( p[i] / rho[i] ) )
        e[i] = e[i] + de_dt[i] * dt

    # Implement steady boundary conditions at x = 0
    v[0] = u0
    rho[0] = 1.
    p[0] = 1.
    e[0] = 0.5 * rho[0] * v[0]**2 + (1 / (gamma - 1)) * p[0]

    # Implement continuous boundary conditions at x = 10
    v[-1] = v[-2]
    rho[-1] = rho[-2]
    p[-1] = p[-2]
    e[-1] = e[-2]

    p = (e - 0.5*rho*v**2)*(gamma-1)
        
    t = t + dt

print(f"t={t}")           

#Here we plot the results
        
fig = plt.figure(dpi = 500)
ax1 = fig.add_subplot(1,1,1, xlim = [x_min, x_max])
plt.plot(x, v, '.')
plt.xlabel('x')
plt.ylabel('Velocity')
plt.title(f"t={round(t,1)}, v = {u0}")

plt.show()

fig = plt.figure(dpi = 500)
ax1 = fig.add_subplot(1,1,1, xlim = [x_min, x_max])
plt.plot(x, rho)
plt.xlabel('x')
plt.ylabel('Density')
plt.title(f"t={round(t,1)}, v = {u0}")

plt.show()

fig = plt.figure(dpi = 500)
ax1 = fig.add_subplot(1,1,1, xlim = [x_min, x_max])
plt.plot(x, p)
plt.xlabel('x')
plt.ylabel('Pressure')
plt.title(f"t={round(t,1)}, v = {u0}")

plt.show()

fig = plt.figure(dpi = 500)
ax1 = fig.add_subplot(1,1,1, xlim = [x_min, x_max])
plt.plot(x, mach)
plt.xlabel('x')
plt.ylabel('Mach Number')
plt.title(f"t={round(t,1)}, v = {u0}")

plt.show()
    
