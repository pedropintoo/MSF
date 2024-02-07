# -*- coding: utf-8 -*-
"""
Created on Thu May 18 17:07:29 2023

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
    
#%% Método de Euler Cromer (Complexo com RA) - Oscilador Harmónico Duplo


x = np.linspace(-6,6,100)
k = 1
xeq = 1.5
Ep = 1/2 * k * (x**2 - xeq**2)**2

plt.plot(x,Ep)
plt.ylim(0,10)
plt.xlabel("x [m]")
plt.ylabel("Energia potencial [J]")


#%% Gravar em PNG 


plt.savefig(dir_figs+'ex2-a.png') # Save the graph in png

plt.show() 

# Abrir um separador com o grafico
# Tools > Preferences > IPython console > Graphics > Graphics Backend > Backend: Automatic

