# -*- coding: utf-8 -*-
"""
Created on Thu Mar 16 16:34:55 2023

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


#%% EX3: Vetor posicao e velocidade


# (x0,y0) = ponto inicial do vetor
# v1 = (x,y) = comprimentos x e y do vetor

x0,y0 = 0, 0

time = np.arange(1, 5, 1)
r = np.array([2*time, time**2])
v = np.array([2*np.ones(len(time)), 2*time])


fig, axs = plt.subplots(1, 2)
axs[0].plot(x0,y0, 'o', markersize = 5)

for i in range(len(time)):
    
    print(r[0,i])
    axs[0].arrow(x0, y0, r[0,i], r[1,i], color='k', width=0.1)
    axs[0].plot(r[0,i],r[1,i], 'ro', markersize = 5)
    axs[1].plot(r[0,i],r[1,i], 'ro', markersize = 5)
    axs[1].arrow(r[0,i], r[1,i], v[0,i], v[1,i], color='k', width=0.1)


axs[0].set_xlim(0,15)
axs[0].set_ylim(0,25)
axs[1].set_xlim(0,15)
axs[1].set_ylim(0,25)
    
axs[0].set_aspect('equal')



plt.show()




#%% Gravar em PNG 


plt.savefig(dir_figs+'Ex6_vetor_posicao_velocidade.png') # Save the graph in png

plt.show() 

# Abrir um separador com o grafico
# Tools > Preferences > IPython console > Graphics > Graphics Backend > Backend: Automatic
