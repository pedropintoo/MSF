# -*- coding: utf-8 -*-
"""
Created on Fri Jun 23 09:38:49 2023

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
    
#%% Oscilador, analise Energia Potencial



x = np.linspace(-8,4,100) # Alterar consoante o problema e o Xeq
k = 1
alpha = 0.05
Ep = 1/2 * k * x**2 + alpha * x**3

EP_MENOR_1 = 7

print("xf (7) =",x[np.where(Ep <= EP_MENOR_1)[0][-1]])

EP_MENOR_2 = 8

print("xf (8) =",x[np.where(Ep <= EP_MENOR_2)[0][-1]])

plt.plot(x,Ep)
plt.plot(x[np.where(Ep <= EP_MENOR_1)[0][-1]], EP_MENOR_1, "-or")
plt.plot(x[np.where(Ep <= EP_MENOR_2)[0][-1]], EP_MENOR_2, "-og")

plt.ylim(0,10)   # Alterar consoante o problema
plt.xlabel("x [m]")
plt.ylabel("Energia potencial [J]")
plt.grid()

# Resposta:
    
# 𝐸𝑝 < 7 J:
# O corpo oscila periodicamente entre as posições em que a 𝐸𝑝 = 7 J. Como a energia potencial não é simétrica à
# volta da posição de equilíbrio, o movimento oscilatório tem uma posição média (por período) < 0.   
    
# 𝐸𝑝 > 8 J:
# O corpo não oscila, nem tem um movimento periodico, pois para x<0 nunca é atingida a
# Ep > 8 J. 




#%% Gravar em PNG 


plt.savefig(dir_figs+'ex1-a).png') # Save the graph in png

plt.show() 

# Abrir um separador com o grafico
# Tools > Preferences > IPython console > Graphics > Graphics Backend > Backend: Automatic



