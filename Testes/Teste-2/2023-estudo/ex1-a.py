# -*- coding: utf-8 -*-
"""
Created on Wed May 10 17:17:46 2023

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
tf = 1.2
dt = 0.0001

g = 9.8

x0 = 0
y0 = 3
z0 = 0


v0x = 30 # 30m/s
v0y = 0
v0z = 0


vt = 20 # 20 m/s


D = g/(vt**2)

Nt = int(np.ceil((tf - t0) / dt) + 1)

t = np.linspace(t0, tf, Nt)

z = np.zeros((Nt,)) # np.zeros((Nt,2))  faria com 2 colunas
z[0] = z0

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

entrou = False;
# Método de Euler
for i in range(Nt - 1):
    
    
    vv=np.sqrt(vx[i]**2 + vy[i]**2)# Intensidade
    
    ax[i]=-D*vv*vx[i]
    ay[i]=-g-D*vv*vy[i]
    az[i]= 0
    

    
    vx[i+1]=vx[i]+ax[i]*dt
    vy[i+1]=vy[i]+ay[i]*dt
    vz[i+1]=0
    
    
    x[i+1]=x[i]+vx[i]*dt
    y[i+1]=y[i]+vy[i]*dt
    z[i+1]=0
    
    
        
# Como z = 0 em todo o tempo, nao consideramos
    
inters_ground = np.where(x >= 12)[0][0]      
if(y[inters_ground] < 1):
    entrou = False # Se bater na rede



inters_ground = np.where(y <= 0)[0][0]
alcance = x[inters_ground]

if (not entrou) and (12 < alcance < 18.4):
    entrou = True # se quando bater no chao estiver no limite é valido

print("Alcance da bola: ",alcance,"m")

if(entrou):
    print("O serviço é valido!")
else:
    print("O serviço não é valido!")


plt.plot(x, y, '-r', linewidth=1) # Reta
plt.xlabel('x (m)')
plt.ylabel('y (m)')
plt.grid()




#%% Gravar em PNG 


plt.savefig(dir_figs+'1-a.png') # Save the graph in png

plt.show() 

# Abrir um separador com o grafico
# Tools > Preferences > IPython console > Graphics > Graphics Backend > Backend: Automatic

