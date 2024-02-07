# -*- coding: utf-8 -*-
"""
Created on Wed Jun 21 10:11:43 2023

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
    
#%% Oscilador quartico



x = np.linspace(-3,3,100)
k = 2
xeq = 0
alpha = -0.1
b = 0.02
Ep = 1/2 * k * x**2 + alpha * x**3 - b * x**4

plt.plot(x,Ep)
plt.ylim(0,10)
plt.xlabel("x [m]")
plt.ylabel("Energia potencial [J]")

  
    
    
#%% Gravar em PNG 


plt.savefig(dir_figs+'ex1-a).png') # Save the graph in png

plt.show() 

# Abrir um separador com o grafico
# Tools > Preferences > IPython console > Graphics > Graphics Backend > Backend: Automatic