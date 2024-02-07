# -*- coding: utf-8 -*-
"""
Created on Tue May  9 23:07:23 2023

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
    
#%% Método de Euler Cromer (Complexo com RA) - Movimento Harmónico Simples


t0 = 0
tf = 60
dt = 0.001

g = 9.8
x0 = 0 
v0 = 0  


m = 1.5  # kg
k = 1.2   # Constante elástica
c = -0.01 # N/m^2

Nt = int(np.ceil((tf - t0) / dt) + 1)

t = np.linspace(t0, tf, Nt)


x = np.zeros((Nt,)) # np.zeros((Nt,2))  faria com 2 colunas
x[0] = x0  

v = np.zeros((Nt,))
v[0] = v0

a = np.zeros((Nt,))

Ep = np.zeros((Nt,))


# Método de Euler Cromer
for i in range(Nt - 1):


    a[i]=-k* x[i] / m  - 3 * c * (x[i]**2) / m  # Fx = -k * x(t) --> a = -k/m * x(t)

    
    v[i+1]=v[i]+a[i]*dt
    
    
    x[i+1]=x[i]+v[i+1]*dt
    
    
    Ep[i] = (1/2) * k * (x[i]**2) + c * (x[i]**3)
    


# Plotting

plt.plot(t, x, '-r', linewidth=1) # Reta



plt.xlabel('t (s)')
plt.ylabel('v (m/s^2)')
plt.grid()


#%% Gravar em PNG 


plt.savefig(dir_figs+'2.png') # Save the graph in png

plt.show() 

# Abrir um separador com o grafico
# Tools > Preferences > IPython console > Graphics > Graphics Backend > Backend: Automatic

