# -*- coding: utf-8 -*-
"""
Created on Tue May  9 22:35:49 2023

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
    
#%% Método de Euler Cromer (Complexo com RA) - (considerando inclinação) Evolução temporal de ciclista com empurao

par = 1.225  # densidade do ar

t0 = 0
tf = 400
dt = 0.001

g = 9.8
x0 = 0
v0 = 0.5  # empurao de 1m/s
a0 = 0  

ang = 4 # graus

pot = 0.48 *  735.4975 # potencia (W)

m = 72  # kg
Catr = 0.01    # coeficiente de resistência 𝜇 de um piso liso de alcatrão
Cres = 0.9  # coeficiente de resistência do ar
A = 0.50 # area frontal

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
    
    # Fx = Femp - Nx
    
    # a = 1/m * (Potencia/velocidade - N * Catr - Px - D * velocidade * v_norma)
    
    Px = m * g * np.sin(np.deg2rad(ang)) # Peso (com inclinação de 5 graus)
    N = m * g * np.cos( np.deg2rad(ang))
    
    a[i] = 1 / m * (pot/v[i] - N * Catr - Px - D * v[i] * v_norma ) 
 
    
    v[i+1]=v[i]+a[i]*dt
    
    
    x[i+1]=x[i]+v[i+1]*dt
    
    if(x[i] > 1500):
        ang = -1    # Desce com inclinacao de 1 grau
        # a = 1/m * (Potencia/velocidade - N * Catr + Px - D * velocidade * v_norma)
    
    
    
# Plotting

plt.plot(t, v, '-r', linewidth=1) # Reta



plt.xlabel('t (s)')
plt.ylabel('v (m/s)')
plt.grid()

tempo_v = np.where(x>=2000)[0][0]
t_tempo_v = t[tempo_v]
print("b) Percorre 2000 metros depois de",t_tempo_v," segundos")



#%% Gravar em PNG 


plt.savefig(dir_figs+'2-c.png') # Save the graph in png

plt.show() 

# Abrir um separador com o grafico
# Tools > Preferences > IPython console > Graphics > Graphics Backend > Backend: Automatic

