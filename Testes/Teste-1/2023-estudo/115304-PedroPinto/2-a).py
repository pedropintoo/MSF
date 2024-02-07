# -*- coding: utf-8 -*-
"""
Created on Wed Mar 22 17:08:30 2023

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


t0 = 0
tf = 2
dt = 0.0001

g = 9.8
x0 = 0
y0 = 2
v0 = 15  # 130 km/h --> m/s
ang = 30  # 10 graus
v0x = v0 * np.cos(np.radians(ang))
v0y = v0 * np.sin(np.radians(ang))


vt = 20 # m/s



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
   
    
    vv=np.sqrt(vx[i]**2 +vy[i]**2) # Intensidade
    
    ax[i]=-D*vv*vx[i]
    ay[i]=-g-D*vv*vy[i]
    
    vx[i+1]=vx[i]+ax[i]*dt
    vy[i+1]=vy[i]+ay[i]*dt
    
    x[i+1]=x[i]+vx[i]*dt
    y[i+1]=y[i]+vy[i]*dt



inters_ground = np.where(y>=3)[0][-1]
alcance = x[inters_ground]
print("Alcance da bola: ",alcance)

tempoal = t[inters_ground]
print("No momento: ",tempoal)



# Plotting

plt.plot(x, y, '-r', linewidth=1) # Reta


plt.xlabel('x (m)')
plt.ylabel('y (m)')
plt.grid()



#%% Gravar em PNG 


plt.savefig(dir_figs+'2-a).png') # Save the graph in png

plt.show() 

# Abrir um separador com o grafico
# Tools > Preferences > IPython console > Graphics > Graphics Backend > Backend: Automatic

