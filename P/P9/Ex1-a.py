# -*- coding: utf-8 -*-
"""
Created on Thu May  4 16:14:19 2023

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
tf = 5
dt = 0.01

g = 9.8
x0 = 0
y0 = 0
v0 = 100 * 1000/3600  # 130 km/h --> m/s
ang = 10  # 10 graus
v0x = v0 * np.cos(np.radians(ang))
v0y = v0 * np.sin(np.radians(ang))


vt = 100 * 1000/3600  # 100 km/h --> m/s

m = 0.057  # kg


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


Ec = np.zeros((Nt,)) # Energia cinética
Ec[0] = (1/2) * m * np.sqrt(v0x**2 + v0y**2)**2

Ep = np.zeros((Nt,))# Energia potencial
Ep[0] = 0 # partiu do solo

Em = np.zeros((Nt,)) # Energia mecanica
Em[0] = Ec[0] + Ep[0]


# Método de Euler Cromer
for i in range(Nt - 1):
        
    vv = np.sqrt(vx[i]**2 + vy[i]**2) #intensidade
    
    
    ax[i]=0
    ay[i]=-g
    
    vx[i+1]=vx[i]+ax[i]*dt
    vy[i+1]=vy[i]+ay[i]*dt
    
    x[i+1]=x[i]+vx[i+1]*dt
    y[i+1]=y[i]+vy[i+1]*dt
    
    
    Ec[i+1] = (1/2) * m * vv**2
    
    Ep[i+1] = m * g * y[i]
    
    # Deveria ser constante pois só tem forcas conservativas
    Em[i+1] = Ec[i+1] + Ep[i+1]



# Plotting
plt.plot(t, Em, '-r', linewidth=1) # Reta


plt.xlabel('t (s)')
plt.ylabel('Em (J)')
plt.grid()



#%% Gravar em PNG 


plt.savefig(dir_figs+'1-a).png') # Save the graph in png

plt.show() 

# Abrir um separador com o grafico
# Tools > Preferences > IPython console > Graphics > Graphics Backend > Backend: Automatic


