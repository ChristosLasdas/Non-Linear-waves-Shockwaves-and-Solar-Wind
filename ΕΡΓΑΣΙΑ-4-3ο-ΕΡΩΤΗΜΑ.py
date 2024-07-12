import numpy as np
  
import matplotlib.pyplot as plt
    
#for n in [50, 100, 1000, 2000, 3000, 3500]: #number of iterations, for 10000 iterations the run takes about 1 minute.
n = 1000

nr=101 #x resolution

r_min=1
r_max=11

rho = np.zeros((nr))
vel=np.zeros((nr))

rho[0] = 1.
for i in range (1,nr):
    rho[i] = 0.01
    
η = 0.1
K = 1.
G_M = 5.

eq_rho = np.zeros((nr))
eq_vel = np.zeros((nr))

r=np.ones(nr)

t = 0.

dr=(r_max-r_min)/(nr-1)

dt=dr*0.01

for i in range(1,nr):
    r[i]=i*dr + 1

#υπολογισμος vel[0]
vel[0] = 1.564

for k in range(0, n):
    for i in range (1, nr-1):
        # 2η παραγωγος ταχύτητας
        #drrvel=(vel[i+1]-2*vel[i]+vel[i-1])/dr**2
        
        #1η παράγωγο πυκνότητας με μέθοδο upwind
        #dr_rho = (rho[i] - rho[i-1])/dr
        #drvel= (vel[i] - vel[i-1])/dr
    
        #eq_rho[i] = - (1/r[i]**2)*(2*r[i]*rho[i]*vel[i] + r[i]**2*vel[i]*dr_rho + r[i]**2*rho[i]*(vel[i+1]-vel[i-1])/(2*dr))
        
        eq_rho[i] = - ( 1 / r[i]**2 ) * ( 2 * r[i] * rho[i] * vel[i]+ r[i]**2 * vel[i]* ( ( rho[i] - rho[i - 1] ) / dr ) \
                             + r[i]**2 * rho[i] * ( vel[i + 1] - vel[i - 1] ) / ( 2 * dr )  )
        
        #eq_vel[i] =  - (vel[i]*drvel + K*(1/rho[i])*(rho[i+1]-rho[i-1])/(2*dr) + G_M/r[i]**2 - η*((2/r[i])*(vel[i+1]-vel[i-1])/(2*dr)) + drrvel)
    
        eq_vel[i] = - vel[i]* ( ( vel[i]- vel[i - 1] ) / dr ) - K / rho[i] * ( ( rho[i + 1] - rho[i - 1] ) / ( 2 * dr ) ) \
                           - G_M / r[i]**2 + η * (2 / r[i] * (vel[i]- vel[i - 1]) / (dr) + ( vel[i + 1] - 2 * vel[i]+ vel[i - 1] ) / ( dr**2 ))
    
    for i in range (2, nr-1):    
        rho[i] = rho[i] + dt*eq_rho[i]
        vel[i] = vel[i] + dt*eq_vel[i]
     
    vel[nr-1] = vel[nr-2]
    rho[nr-1] = rho[nr-2]
    t=t+dt
print(t)           


#Here we plot the results
        
fig = plt.figure(dpi=1000)
ax1 = fig.add_subplot(1,1,1, aspect=1, xlim=[r_min, r_max])

plt.plot(r, vel)
plt.xlabel('r')
plt.ylabel('Velocity')

plt.show()


