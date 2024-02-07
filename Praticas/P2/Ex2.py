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


#%% EX2: Ciclista


distancia = np.array([ 0.00, 0.735, 1.363, 1.739, 2.805, 3.814, 4.458, 4.955, 5.666, 6.329 ])
tempo = np.array(     [    0,     1,     2,     3,     4,     5,     6,     7,     8,     9 ])

plt.plot(tempo, distancia, marker='s', linestyle='none', color='purple')
plt.xlabel('tempo (min)')
plt.ylabel('distancia (km)')
plt.grid()


# Limite de Plot in X and Y
# plt.xlim(ax, bx) 
# plt.ylim(cy, dy)

plt.title("Ex2 - Ciclista")


m, b, r2, delta_m, delta_b = lin_reg(tempo,distancia)
normal_fit = m*tempo + b

plt.plot(tempo, normal_fit, '-g', linewidth=1)


#------------  With polyfit function  ------------#

#a, b = np.polyfit(tempo,distancia,1)
#y = a*tempo + b
#plt.plot(tempo,y,'y')



#------------  Colocar caixas de texto e gravar em PNG  ------------#

textstr = f'$distancia = {m:.3f} \cdot tempo + ({b:.3f}) $ \n'+f'$r^2={r2:.3f}$'
props = dict(boxstyle='round', facecolor='lightgreen', alpha=.5)
plt.text(min(tempo), max(distancia), textstr,
        verticalalignment='top', bbox=props)


plt.savefig(dir_figs+'Ex2_ciclista.png') # Save the graph in png
plt.show()

# Abrir um separador com o grafico
# Tools > Preferences > IPython console > Graphics > Graphics Backend > Backend: Automatic


