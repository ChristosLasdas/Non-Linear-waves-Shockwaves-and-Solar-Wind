
import numpy as np

import matplotlib.pyplot as plt



for n in [0, 10, 20, 30, 60, 90, 120, 150, 180, 210]: 

    nx=100 #x resolution
    
    
    x_min=0.
    x_max=10
    
    vel=np.zeros((nx))
    
    Eq=np.zeros((nx))
    
    #Result=np.zeros((2,nx))
    
    σ = 2.
    xm = 5.
    η = 0.01
    
    x=np.zeros(nx)
    
    t = 0.
    
    dx=(x_max-x_min)/(nx)
    
    
    dt=dx*0.1
    
    
    for i in range(0,nx):
        x[i]=i*dx
    
    
    for i in range(0,nx-2):
        vel[i]= np.exp(-((x[i]-xm)/σ)**2)
    #print (vel)
    
    vel[99] = vel[1]
    vel[98] = vel[0]
    
    for k in range(0, n):
        
        for i in range (1, nx-1):
            if vel[i] > 0:
    
                #1η παράγωγος ταχύτητας
                dxvel= (vel[i] - vel[i-1])/(dx)
    
                # 2η παραγωγος ταχύτητας
                dxxvel=(vel[i+1]-2*vel[i]+vel[i-1])/dx**2
                #dxxvel=(vel[i] - 2*vel[i-1] + vel[i-2])/dx**2
            elif vel[i] < 0:
            
                dxvel = (vel[i+1] - vel[i])/dx
                dxxvel=(vel[i+1]-2*vel[i]+vel[i-1])/dx**2
                #dxxvel = (vel[i+2] - 2*vel[i+1] + vel[i])/dx**2
               # print ('αρνητικο')
                #print (i)
                #print (k)
                
            Eq[i]= - vel[i]*dxvel + η*dxxvel
        
           
        for i in range (1, nx-1):
            vel[i]=vel[i]+Eq[i]*dt
        
        vel[nx-1] = vel[1]
        vel[0] = vel[nx-2]   
        t=t+dt
    #In this loop we save the results in a single array
    print(t)           
    
    #Here we plot the results
            
    fig = plt.figure(dpi=1000)
    ax1 = fig.add_subplot(1,1,1, aspect=1, xlim=[x_min, x_max], ylim=[0 , 2])
    plt.plot(x, vel, '.')
    plt.xlabel('x')
    plt.ylabel('Velocity')
    
    plt.show()
    
