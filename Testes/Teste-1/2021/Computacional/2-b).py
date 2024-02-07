# -*- coding: utf-8 -*-
"""
Created on Tue Mar 21 12:34:35 2023

@author: pedro
"""


#%% Importar libraries e diretorio para figuras

import matplotlib.pyplot as plt
import numpy as np
import os

dir_figs = './figures/'
if not os.path.exists(dir_figs): 
    os.mkdir(dir_figs)


#%% Método de Euler (Complexo com RA) - 2 dimensão

PI = np.pi
p_ar = 1.225  # kg/m3

t0 = 0
tf = 2
dt = 0.001

g = 9.8
x0 = -10
y0 = 1
v0 = 130 * 1000/3600  # 130 km/h --> m/s
ang = 10  # 10 graus
v0x = v0 * np.cos(np.radians(ang))
v0y = v0 * np.sin(np.radians(ang))


vt = 100 * 1000/3600  # 100 km/h --> m/s

m = 0.057  # kg
d = 0.067  # m   # Diametro
r = d / 2  # m   # Raio
A = PI * r**2


D = g/(vt**2)

Nt = int(np.ceil((tf - t0) / dt) + 1)

t = np.linspace(t0, tf, Nt)

y = np.zeros((Nt,)) # np.zeros((Nt,2))  faria com 2 colunas
y[0] = y0

x = np.zeros((Nt,)) # np.zeros((Nt,2))  faria com 2 colunas
x[0] = x0  

vx = np.zeros((Nt,))
vx[0] = v0x
vy = np.zeros((Nt,))
vy[0] = v0y

ax = np.zeros((Nt,))
ay = np.zeros((Nt,))


# Método de Euler
for i in range(Nt - 1):
    # 𝐹_𝑀𝑎𝑔𝑛𝑢s = 1/2 𝐴 𝜌_𝑎𝑟 𝑟 𝜔⃗ × 𝑣 
    # ang_w é o ang da rotacao
    
    ang_w = (0,0,100)
    F_magnus = (1/2) * A * p_ar * r * (np.cross(ang_w,(vx[i],vy[i],0)))
    ax_magnus = F_magnus[0]/m
    ay_magnus = F_magnus[1]/m
    
    vv=np.sqrt(vx[i]**2 +vy[i]**2) # Intensidade
    
    ax[i]=-D*vv*vx[i] + ax_magnus
    ay[i]=-g-D*vv*vy[i] + ay_magnus
    
    vx[i+1]=vx[i]+ax[i]*dt
    vy[i+1]=vy[i]+ay[i]*dt
    
    x[i+1]=x[i]+vx[i]*dt
    y[i+1]=y[i]+vy[i]*dt


# Max position

if(len(y == np.max(y)) >= 0):
    y_max = np.max(y)
    inters_max = np.where(y == np.max(y))[0][0]
    t_inters_max = t[inters_max]
    print("Altura máxima da bola: y = ",y_max,"em t =", )

if(len(np.where(y<=0)[0]) > 1):
    inters_ground = np.where(y<=0)[0][0]
    alcance = x[inters_ground]
    print("Alcance da bola: ",alcance)



# Plotting

plt.plot(x, y, '-r', linewidth=1) # Reta


plt.xlabel('x (m)')
plt.ylabel('y (m)')
plt.grid()

plt.legend(['posicao sem RA','posicao com RA'])


#%% Gravar em PNG 


plt.savefig(dir_figs+'2-b).png') # Save the graph in png

plt.show() 

# Abrir um separador com o grafico
# Tools > Preferences > IPython console > Graphics > Graphics Backend > Backend: Automatic
