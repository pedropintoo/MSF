
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


#%% Método de Euler (Complexo sem RA) - 2 dimensão


t0 = 0
tf = 2
dt = 0.1

g = 9.8
x0 = 0
y0 = 0
v0 = 100 * 1000 / 3600  # 130 km/h --> m/s
ang = 10  # 10 graus
v0x = v0 * np.cos(np.radians(ang))
v0y = v0 * np.sin(np.radians(ang))



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
   
    ax[i]=0
    ay[i]=-g
    
    vx[i+1]=vx[i]+ax[i]*dt
    vy[i+1]=vy[i]+ay[i]*dt
    
    x[i+1]=x[i]+vx[i]*dt
    y[i+1]=y[i]+vy[i]*dt



# Plotting

fig, axs = plt.subplots(1, 2,figsize=(5,7)) # sharex=True faria com que o x fosse partilhado


axs[0].plot(x, t, '-r', linewidth=1) # Reta
axs[0].set_xlabel('t (s)')
axs[0].set_ylabel('x (m)')
axs[0].grid()

axs[1].plot(y, t, '-r', linewidth=1) # Reta
axs[1].set_xlabel('t (s)')
axs[1].set_ylabel('y (m)')
axs[1].grid()





#%% Gravar em PNG 


plt.savefig(dir_figs+'1-a).png') # Save the graph in png

plt.show() 

# Abrir um separador com o grafico
# Tools > Preferences > IPython console > Graphics > Graphics Backend > Backend: Automatic

