# -*- coding: utf-8 -*-
"""
Created on Tue Mar 21 12:01:23 2023

@author: pedro
"""

#%% Importar libraries e diretorio para figuras

import matplotlib.pyplot as plt
import numpy as np
import os

dir_figs = './figures/'
if not os.path.exists(dir_figs): 
    os.mkdir(dir_figs)

#%% Função para regressão linear

def lin_reg(x, y, print_res=0):
    
    if len(x) != len(y):
        raise ValueError('ERROR: x and y must be the same length in lin_reg(x, y)!')
    N = len(x)
    if N < 3:
        raise ValueError('ERROR: N must be higher than 2')
    
    # summations
    s_x = np.sum(x)
    s_y = np.sum(y)
    s_xy = np.sum(x * y)
    s_x2 = np.sum(x ** 2)
    s_y2 = np.sum(y ** 2)
    
    # linear regression parameters
    m = (N * s_xy - s_x * s_y) / (N * s_x2 - s_x ** 2)
    b = (s_x2 * s_y - s_x * s_xy) / (N * s_x2 - s_x ** 2)
    
    # squared correlation coefficient
    r2 = ((N * s_xy - s_x * s_y)**2) / ((N * s_x2 - s_x ** 2)*(N * s_y2 - s_y ** 2))

    # errors
    delta_m = np.abs(m) * np.sqrt((1 / r2 - 1) / (N - 2))
    delta_b = delta_m * np.sqrt(s_x2 / N)
    
    if print_res == 1:    # For printing
        print(f"m = {m}")
        print(f"b = {b}")
        print(f"r² = {r2}")
        print(f"delta_m = {delta_m}")
        print(f"delta_b = {delta_b}")
        
    
    return m, b, r2, delta_m, delta_b

#%% Regressão linear (normal) Método dos mínimos quadrados

massa = np.array([ 0.15, 0.20, 0.16, 0.11, 0.25, 0.32, 0.40, 0.45, 0.50, 0.55 ])  # X
periodo = np.array([ 1.21, 1.40, 1.26, 1.05, 1.60, 1.78, 2.00, 2.11, 2.22, 2.33]) # Y

log_massa = np.log(massa)
log_periodo = np.log(periodo)


plt.plot(log_massa, log_periodo, marker='s', linestyle='none', color='blue')
plt.xlabel('log(M)')
plt.ylabel('log(T)')
plt.grid()


m, b, r2, delta_m, delta_b = lin_reg(log_massa, log_periodo, 1)
log_fit = m*log_massa + b

plt.plot(log_massa, log_fit,'g', linewidth=1)


#%% Gravar em PNG 


plt.savefig(dir_figs+'1-b).png') # Save the graph in png

plt.show() 

# Abrir um separador com o grafico
# Tools > Preferences > IPython console > Graphics > Graphics Backend > Backend: Automatic
