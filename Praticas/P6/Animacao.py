# -*- coding: utf-8 -*-
"""
Created on Thu Mar 23 16:31:04 2023

@author: pedro
"""

#%% Importar libraries e diretorio para figuras

import matplotlib.pyplot as plt
import numpy as np
import os
import matplotlib.animation as animation
import time

dir_figs = './figures/'
if not os.path.exists(dir_figs): 
    os.mkdir(dir_figs)


#%% Método de Euler + Força de Magnus (Complexo com RA) - 3 dimensão


## ANIMATION

# The maximum z-range of ball's trajectory to plot.
ZMAX = 5
# The coefficient of restitution for bounces (-v_up/v_down).
cor = 0.65





################


PI = np.pi
p_ar = 1.225  # kg/m3

t0 = 0
tf = 1.5
dt = 0.001

g = 9.8

x0 = 0
y0 = 0.11
z0 = 45


v0x = 25
v0y = 7
v0z = -100


vt = 100 * 1000/3600  # 100 km/h --> m/s

m = 0.45  # kg
r = 0.11  # m   # Raio
A = PI * r**2


D = g/(vt**2)

Nt = int(np.ceil((tf - t0) / dt) + 1)

t = np.linspace(t0, tf, Nt)

y = np.zeros((Nt,)) # np.zeros((Nt,2))  faria com 2 colunas
y[0] = y0

x = np.zeros((Nt,)) # np.zeros((Nt,2))  faria com 2 colunas
x[0] = x0  

z = np.zeros((Nt,)) # np.zeros((Nt,2))  faria com 2 colunas
z[0] = z0  

vx = np.zeros((Nt,))
vx[0] = v0x
vy = np.zeros((Nt,))
vy[0] = v0y
vz = np.zeros((Nt,))
vz[0] = v0z

ax = np.zeros((Nt,))
ay = np.zeros((Nt,))
az = np.zeros((Nt,))

c1 = 0
c2 = 0
ang_w = (0,200,0)
# Método de Euler
for i in range(Nt - 1):
    # 𝐹_𝑀𝑎𝑔𝑛𝑢s = 1/2 𝐴 𝜌_𝑎𝑟 𝑟 𝜔⃗ × 𝑣 
    # ang_w é o ang da rotacao
    
    
    # F_magnus = (1/2) * A * p_ar * r * (np.cross(ang_w,(vx[i],vy[i],vz[i])))
    
    # Neste caso !!!
    F_magnus = (1/2) * A * p_ar * r * (np.cross(ang_w,(vx[i],vy[i],vz[i])))
    
    ax_magnus = F_magnus[0]/m
    ay_magnus = F_magnus[1]/m  # 0
    az_magnus = F_magnus[2]/m
    
    vv=np.sqrt(vx[i]**2 + vy[i]**2 + vz[i]**2) # Intensidade
    
    ax[i]=-D*vv*vx[i] + ax_magnus
    ay[i]=-g-D*vv*vy[i] + ay_magnus
    az[i]=-D*vv*vz[i] + az_magnus

    
    vx[i+1]=vx[i]+ax[i]*dt
    vy[i+1]=vy[i]+ay[i]*dt
    vz[i+1]=vz[i]+az[i]*dt
    
    x[i+1]=x[i]+vx[i]*dt
    y[i+1]=y[i]+vy[i]*dt
    z[i+1]=z[i]+vz[i]*dt
    
    
    
    if x[i] < 0 and c1 <= 10 :
        c1 += 1
        print("a")
        # bounce!
        x[i+1] = 0
        vx[i+1] = -vx[i] * cor
        ang_w = (-20,-20,20)
    
    if y[i] < 0 and c2 <= 10:
        c2 += 1
        # bounce!
        y[i+1] = 0
        vy[i+1] = -vy[i] * cor
        ang_w = (-50,-1,-50)
    
    
        
    




# Plotting

# fig, axs = plt.subplots(3, 1,figsize=(5,7)) # sharex=True faria com que o x fosse partilhado


# axs[0].plot(t, x, '-r', linewidth=1) # Reta
# axs[0].set_xlabel('t (s)')
# axs[0].set_ylabel('x (m)')
# axs[0].grid()

# axs[1].plot(t, y, '-r', linewidth=1) # Reta
# axs[1].set_xlabel('t (s)')
# axs[1].set_ylabel('y (m)')
# axs[1].grid()

# axs[2].plot(t, z, '-r', linewidth=1) # Reta
# axs[2].set_xlabel('t (s)')
# axs[2].set_ylabel('z (m)')
# axs[2].grid()



plt.figure(figsize=(8,8))
ax = plt.axes(projection='3d')
goalx = [0,0,0,0]
goaly = [0,2.4,2.4,0]
goalz = [-7.3,-7.3,7.3,7.3]
ax.plot3D(goalx,goalz,goaly, 'k')
ax.set_xlim3d(0, 15)
ax.set_ylim3d(-45, 45)
ax.set_zlim3d(0, 5)
ax.set_box_aspect((2,4,1))
ax.set_xlabel('x')
ax.set_ylabel('z')
ax.set_zlabel('y')

a = 0
while True:
    
    
    if(a == 0):
        i = 0
        point1 = ax.scatter(x0,-z0, y0, c="black",s=1)
        for i in range(Nt - 1):
            
            i += 1
            if(i == 1):
                plt.pause(0.5)
            if(x[i]>=0 and i%20 == 0):
                
                point1.remove()
                point1 = ax.scatter(x[i],-z[i], y[i], c='black',s=50)
                plt.pause(0.05)
        a = 1
    if(a == 1):
        i = 0
        point2 = ax.scatter(x[-1],-z[-1], y[-1], c="red",s=1)
        for i in range(Nt - 1):
            
            i += 1
            if(i == 1):
                plt.pause(0.5)
            if(x[i]>=0 and i%20 == 0):
                
                point2.remove()
                point2 = ax.scatter(x[-1-i],-z[-1-i], y[-1-i], c='red',s=50)
                plt.pause(0.0005)        
        a = 0
    
    point1.remove()
    point2.remove()
    





#%% Gravar em PNG 


plt.savefig(dir_figs+'2.png') # Save the graph in png

plt.show() 

# Abrir um separador com o grafico
# Tools > Preferences > IPython console > Graphics > Graphics Backend > Backend: Automatic

