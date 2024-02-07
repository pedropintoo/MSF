# -*- coding: utf-8 -*-
"""
Created on Wed May 10 17:43:03 2023

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
#%% Aproximação retangular

 # Aproximação retangular:
 # I = dx * np.sum(f[0:N]) = dt * np.sum(f[0:N] * vx[0:N])

def integral_retangular(f, v, iInicial, iFinal, dt):
    N = int((iFinal - iInicial) / dt)
    
    if N < 1:
        raise ValueError('Número de subintervalos inválido')
        
    integral = dt * sum(f[0:N]*v[0:N])
    return integral

#%% Método de Euler Cromer (Complexo com RA) - (considerando inclinação) Evolução temporal 

par = 1.225  # densidade do ar

t0 = 0
tf = 150
dt = 0.001

g = 9.8
x0 = 0
v0 = 1  # empurao de 1m/s
a0 = 0  

ang = 5 # graus

pot = 40000 # 40k potencia (W)

m = 2000  # kg
Catr = 0.04    # coeficiente de resistência 𝜇 de um piso liso de alcatrão
Cres = 0.25  # coeficiente de resistência do ar
A = 2 # area frontal

D = (Cres/2) * A * par



Nt = int(np.ceil((tf - t0) / dt) + 1)

t = np.linspace(t0, tf, Nt)


x = np.zeros((Nt,)) # np.zeros((Nt,2))  faria com 2 colunas
x[0] = x0  

v = np.zeros((Nt,))
v[0] = v0

a = np.zeros((Nt,))
a[0] = a0


Fmot = np.zeros((Nt,)) # Forca de resistencia do ar segundo X



# Método de Euler Cromer
for i in range(Nt - 1):
    
    v_norma = np.sqrt(v[i]**2 - 0) # Norma da velocidade

    # D = D = (Cres/2) * A * par
    # Px = m * g * sin(ang)  # Peso segundo Ox (com inclinação de ang)
    # N = m * g * cos(ang)   # Forca normal (com inclinação de ang)
        
    # a = 1/m * (Potencia/velocidade - N * Catr - Px - D * velocidade * v_norma)
    
    Px = m * g * np.sin(np.deg2rad(ang)) # Peso (com inclinação de 5 graus)
    N = m * g * np.cos( np.deg2rad(ang))
    Fmot[i] = pot/v[i]
    
    a[i] = 1 / m * (pot/v[i] - N * Catr - Px - D * v[i] * v_norma ) 
 
    
    v[i+1]=v[i]+a[i]*dt
    
    
    x[i+1]=x[i]+v[i+1]*dt
    

    
    
    
    
tempo_v = np.where(x>=2000)[0][0]
t_tempo_v = t[tempo_v]

print("b) Percorre 2000 metros depois de",t_tempo_v," segundos")

# Calcular o trabalho da Forca do motor até aos 2km

iInicial = 0
iFinal = t[tempo_v]

# Aprox. retangular
Wmotor = integral_retangular(Fmot, v, iInicial, iFinal, dt)


print("c) Trabalho realizado é: ",Wmotor)


    
    
    
# Plotting

fig, axs = plt.subplots(2, 1,figsize=(5,7)) # sharex=True faria com que o x fosse partilhado


axs[0].plot(t, x, '-r', linewidth=1) # Reta
axs[0].set_xlabel('t (s)')
axs[0].set_ylabel('x (m)')
axs[0].grid()

axs[1].plot(t, v, '-r', linewidth=1) # Reta
axs[1].set_xlabel('t (s)')
axs[1].set_ylabel('v (m/s)')
axs[1].grid()



#%% Gravar em PNG 


plt.savefig(dir_figs+'2.png') # Save the graph in png

plt.show() 

# Abrir um separador com o grafico
# Tools > Preferences > IPython console > Graphics > Graphics Backend > Backend: Automatic

