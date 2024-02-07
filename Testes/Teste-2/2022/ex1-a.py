# -*- coding: utf-8 -*-
"""
Created on Tue May  9 22:15:21 2023

@author: pedro
"""
#%% Importar libraries e diretorio para figuras

import matplotlib.pyplot as plt
import numpy as np
import sympy as sy
import os

dir_figs = './figures/'
if not os.path.exists(dir_figs): 
    os.mkdir(dir_figs)
    
#%% Método de Euler + Força de Magnus + Resistencia do Ar (Complexo com RA) - 3 dimensão


PI = np.pi
p_ar = 1.225  # kg/m3

t0 = 0
tf = 1.5
dt = 0.00001

g = 9.8

x0 = 0
y0 = 0

ang = 16

v0 = 100 * 1000 / 3600 # v0 = 100km/h --> m/s

v0x = v0 * np.cos(np.deg2rad(ang))
v0y = v0 * np.sin(np.deg2rad(ang))


vt = 100 * 1000/3600  # 100 km/h --> m/s


D = g/(vt**2)

Nt = int(np.ceil((tf - t0) / dt) + 1)

t = np.linspace(t0, tf, Nt)

y = np.zeros((Nt,)) # np.zeros((Nt,2))  faria com 2 colunas
y[0] = y0

x = np.zeros((Nt,)) # np.zeros((Nt,2))  faria com 2 colunas
x[0] = x0  

z = 0


vx = np.zeros((Nt,))
vx[0] = v0x
vy = np.zeros((Nt,))
vy[0] = v0y


ax = np.zeros((Nt,))
ay = np.zeros((Nt,))
#az = np.zeros((Nt,))

entrou = False;
# Método de Euler
for i in range(Nt - 1):
    
    
    vv=np.sqrt(vx[i]**2 + vy[i]**2)# Intensidade
    
    ax[i]=-D*vv*vx[i] #+ ax_magnus
    ay[i]=-g-D*vv*vy[i] #+ ay_magnus
    

    
    vx[i+1]=vx[i]+ax[i]*dt
    vy[i+1]=vy[i]+ay[i]*dt
    
    
    x[i+1]=x[i]+vx[i]*dt
    y[i+1]=y[i]+vy[i]*dt
   
    if(not entrou) and x[i] > 20 :
        print("x = ",x[i])
        print("y = ",y[i])
        print("z = 0")
        if (-3.75 < z < 3.75 and 0 < y[i] < 2.4): 
            print("A bola entrou!")
            entrou = True;
        else:
            print("A bola nao entrou!")
        break;
        
        

# Plotting

fig, axs = plt.subplots(2, 1,figsize=(5,7)) # sharex=True faria com que o x fosse partilhado


axs[0].plot(t, x, '-r', linewidth=1) # Reta
axs[0].set_xlabel('t (s)')
axs[0].set_ylabel('x (m)')
axs[0].grid()

axs[1].plot(t, y, '-r', linewidth=1) # Reta
axs[1].set_xlabel('t (s)')
axs[1].set_ylabel('y (m)')
axs[1].grid()



# plt.figure(figsize=(8,8))
# ax = plt.axes(projection='3d')
# ax.plot3D(x[x>=0],-z[x>=0],y[x>=0], 'r')
# goalx = [0,0,0,0]
# goaly = [0,2.4,2.4,0]
# goalz = [-3.66,-3.66,3.66,3.66]
# ax.plot3D(goalx,goalz,goaly, 'k')
# ax.set_xlim3d(0, 15)
# ax.set_ylim3d(-25, 25)
# ax.set_zlim3d(0, 5)
# ax.set_box_aspect((2,6,2))
# ax.set_xlabel('x')
# ax.set_ylabel('z')
# ax.set_zlabel('y')

#%% Gravar em PNG 


plt.savefig(dir_figs+'1-a.png') # Save the graph in png

plt.show() 

# Abrir um separador com o grafico
# Tools > Preferences > IPython console > Graphics > Graphics Backend > Backend: Automatic


