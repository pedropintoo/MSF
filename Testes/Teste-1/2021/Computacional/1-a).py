# -*- coding: utf-8 -*-
"""
Created on Tue Mar 21 11:50:07 2023

@author: pedro
"""

#%% Importar libraries e diretorio para figuras

import matplotlib.pyplot as plt
import numpy as np
import os

dir_figs = './figures/'
if not os.path.exists(dir_figs): 
    os.mkdir(dir_figs)
    

#%% Mostrar gráficamente apenas

massa = np.array([ 0.15, 0.20, 0.16, 0.11, 0.25, 0.32, 0.40, 0.45, 0.50, 0.55 ])  # X
periodo = np.array([ 1.21, 1.40, 1.26, 1.05, 1.60, 1.78, 2.00, 2.11, 2.22, 2.33]) # Y


plt.plot(massa, periodo, marker='s', linestyle='none', color='purple')
plt.xlabel('M (kg)')
plt.ylabel('T (s)')
plt.grid()





#%% Gravar em PNG 


plt.savefig(dir_figs+'1-a).png') # Save the graph in png

plt.show() 

# Abrir um separador com o grafico
# Tools > Preferences > IPython console > Graphics > Graphics Backend > Backend: Automatic