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


#%% EX3: Vetor perpendicular


# (x0,y0) = ponto inicial do vetor
# v1 = (x,y) = comprimentos x e y do vetor

x0,y0 = 0, 0
v1 = [3, 4] 
v2 = [1, -3/4]

fig, ax = plt.subplots()
ax.plot(x0,y0, 'o', markersize = 5)
ax.arrow(x0,y0,v1[0],v1[1],color='r',width=0.05)
ax.arrow(x0,y0,v2[0],v2[1],color='r',width=0.05)
ax.set_aspect('equal')
plt.show()



#%% Gravar em PNG 


plt.savefig(dir_figs+'Ex1_bola_verticalmente.png') # Save the graph in png

plt.show() 

# Abrir um separador com o grafico
# Tools > Preferences > IPython console > Graphics > Graphics Backend > Backend: Automatic
