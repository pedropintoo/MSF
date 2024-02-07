# -*- coding: utf-8 -*-
"""
Created on Thu Mar  2 16:07:54 2023

@author: pedro
"""
#%% Importar libraries e diretorio para figuras

import matplotlib.pyplot as plt
import numpy as np
import os

dir_figs = './figures/'
if not os.path.exists(dir_figs): 
    os.mkdir(dir_figs)


#%% EX1: Carro vs Patrulha

t0 = 0
tf = 30
dt = 0.01

v0_a = 70 * 1000/3600 # Velocidade constante de 70km/h para m/s
a0_p = 2              # m/s^2

Nt = int(np.ceil((tf - t0) / dt) + 1)

t = np.linspace(t0, tf, Nt)
# t = np.arange(t0,tf+dt,dt)

# Carro A
x_a = np.zeros((Nt,)) # np.zeros((Nt,2))  faria com 2 colunas
v_a = np.zeros((Nt,))
a_a = np.zeros((Nt,))

# Carro patrulha
x_p = np.zeros((Nt,))
v_p = np.zeros((Nt,))
a_p = np.zeros((Nt,))


# Método de Euler
for i in range(Nt - 1):
    # Carro A
    v_a[i + 1] = v0_a # Porque é constante (acelaracao = 0)
    x_a[i + 1] = x_a[i] + v_a[i]*dt
    
    # Carro P
    a_p[i + 1] = a0_p
    v_p[i + 1] = v_p[i] + a_p[i]*dt
    x_p[i + 1] = x_p[i] + v_p[i]*dt
    


# Carro A
plt.plot(t, x_a, '-g', linewidth=1) # Reta

# Carro patrulha
plt.plot(t, x_p, linewidth=1) # Reta

# Interseção

inters = np.where(x_p > x_a)[0][0]
t_inters = t[inters]

print("Os gráficos intersetam-se em t = ",t_inters)


plt.xlabel('t (s)')
plt.ylabel('x (m)')
plt.grid()

plt.legend(['Carro A','Carro patrulha'])

plt.title("Ex1 - Carro vs Patrulha")




#------------  Gravar em PNG  ------------#


plt.savefig(dir_figs+'Ex1_carro_vs_patrulha.png') # Save the graph in png

plt.show() 

# Abrir um separador com o grafico
# Tools > Preferences > IPython console > Graphics > Graphics Backend > Backend: Automatic




