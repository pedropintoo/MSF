# -*- coding: utf-8 -*-
"""
Created on Thu Apr 20 17:23:58 2023

@author: pedro
"""

#%% Importar libraries e diretorio para figuras

import matplotlib.pyplot as plt
import numpy as np
import os

dir_figs = './figures/'
if not os.path.exists(dir_figs): 
    os.mkdir(dir_figs)


#%% Método de Euler Cromer (Complexo com RA) - Evolução temporal de ciclista com empurao

par = 1.225  # densidade do ar

t0 = 0
tf = 200
dt = 0.001

g = 9.8
x0 = 0
v0 = 1  # empurao de 1m/s
a0 = 0  

pot = 298 # potencia (W)

m = 75  # kg
Catr = 0.004    # coeficiente de resistência 𝜇 de um piso liso de alcatrão
Cres = 0.9  # coeficiente de resistência do ar
A = 0.60 # area frontal

D = (Cres/2) * A * par

Nt = int(np.ceil((tf - t0) / dt) + 1)

t = np.linspace(t0, tf, Nt)


x = np.zeros((Nt,)) # np.zeros((Nt,2))  faria com 2 colunas
x[0] = x0  

v = np.zeros((Nt,))
v[0] = v0

a = np.zeros((Nt,))
a[0] = a0


# Método de Euler Cromer
for i in range(Nt - 1):
    
    v_norma = np.sqrt(v[i]**2 - 0) # Norma da velocidade

    # D = D = (Cres/2) * A * par
    # Px = m * g * sin(ang)  # Peso segundo Ox (com inclinação de ang)
    # N = m * g * cos(ang)   # Forca normal (com inclinação de ang)
    
    # a = 1/m * (Potencia/velocidade - N * Catr - Px - D * velocidade * v_norma)

    a[i] = 1 / m * (pot/v[i] - m * g * np.cos( np.deg2rad(0)) * Catr - D * v[i] * v_norma )  # Px = 0 neste caso e N = m*g porque ang = 0
 
    
    v[i+1]=v[i]+a[i]*dt
    
    
    x[i+1]=x[i]+v[i+1]*dt
    
# Plotting

plt.plot(t, v, '-r', linewidth=1) # Reta



plt.xlabel('t (s)')
plt.ylabel('v (m/s)')
plt.grid()


#%% Resposta das perguntas 


# P = F * V

v = 30 * 1000 / 3600
P = 298
F = P / v
print("a) ",F)

v = 40 * 1000 / 3600
P = 298
F = P / v
print("b) ",F)




#%% Gravar em PNG 


plt.savefig(dir_figs+'3.png') # Save the graph in png

plt.show() 

# Abrir um separador com o grafico
# Tools > Preferences > IPython console > Graphics > Graphics Backend > Backend: Automatic