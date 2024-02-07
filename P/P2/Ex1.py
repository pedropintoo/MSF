# -*- coding: utf-8 -*-
"""

MSF - Aula prática 1
Pedro Pinto

"""

#%% Importar libraries e diretorio para figuras

import matplotlib.pyplot as plt
import numpy as np
import os

dir_figs = './figures/'
if not os.path.exists(dir_figs): 
    os.mkdir(dir_figs)

#%% Funcoes

# linear regression function
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


#%% EX1: EX1: Difração fenda única



L = np.array([222.0, 207.5, 194.0, 171.5, 153.0, 133.0, 113.0, 92.0]) # for X
X = np.array([2.3  , 2.2  , 2.0  , 1.8  , 1.6  , 1.4  , 1.2  , 1.0 ]) # for Y


#print(np.corrcoef(X,L)[0][1]**2) # Com o numpy o [0][1] == [0,1]


plt.plot(L, X, marker='s', linestyle='none', color='purple') # Pontos
plt.xlabel('L (cm)')
plt.ylabel('X (cm)')
plt.grid()

plt.title("Ex1 - Difração fenda única")

m, b, r2, delta_m, delta_b = lin_reg(L,X)
y = m*L + b

plt.plot(L, y, '-g', linewidth=1) # Reta


#------------  Colocar caixas de texto e gravar em PNG  ------------#

textstr = f'$X = {m:.3f} \cdot L + ({b:.3f}) $ \n'+f'$r^2={r2:.3f}$'
props = dict(boxstyle='round', facecolor='lightgreen', alpha=.5)
plt.text(min(L), max(X), textstr,
        verticalalignment='top', bbox=props)


plt.savefig(dir_figs+'Ex1_fenda_unica_fit.png') # Save the graph in png

plt.show()

# Abrir um separador com o grafico
# Tools > Preferences > IPython console > Graphics > Graphics Backend > Backend: Automatic


