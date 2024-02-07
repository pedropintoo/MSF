# -*- coding: utf-8 -*-
"""
Created on Wed Jun 21 14:48:40 2023

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
k = 1.2
alpha = -0.01
Ep = 1/2 * k * x**2 + alpha * x**3

print("xi =",x[np.where(Ep <= 2)[0][0]])
print("xf =",x[np.where(Ep <= 2)[0][-1]])


plt.plot(x,Ep)
plt.plot(x[np.where(Ep <= 2)[0][0]], 2, "-or")
plt.plot(x[np.where(Ep <= 2)[0][-1]], 2, "-or")
plt.ylim(0,5)
plt.xlabel("x [m]")
plt.ylabel("Energia potencial [J]")
plt.grid()

  
    
    
#%% Gravar em PNG 


plt.savefig(dir_figs+'ex4-a).png') # Save the graph in png

plt.show() 

# Abrir um separador com o grafico
# Tools > Preferences > IPython console > Graphics > Graphics Backend > Backend: Automatic