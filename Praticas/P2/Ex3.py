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


#%% EX3: Potencia de Corpo Negro



temperatura = np.array([200,   300,   400,   500,   600,   700,   800,   900,  1000,  1100 ])   # X
energia = np.array([ 0.6950, 4.363, 15.53, 38.74, 75.08, 125.2, 257.9, 344.1, 557.4, 690.7])    # Y

log_temperatura = np.log(temperatura)
log_energia = np.log(energia)


fig, axs = plt.subplots(2, 1,figsize=(5,7)) # sharex=True faria com que o x fosse partilhado


m, b, r2, delta_m, delta_b = lin_reg(log_temperatura, log_energia, 1)
log_fit = m*log_temperatura + b

axs[1].plot(log_temperatura, log_fit,'g', linewidth=1)

# Converting log_fit --> normal_fit

temperatura = temperatura**m

m, b, r2, delta_m, delta_b = lin_reg(temperatura, energia, 1)
normal_fit = m*temperatura + b

axs[0].plot(temperatura, normal_fit,'g', linewidth=1)


axs[0].plot(temperatura, energia, marker='s', linestyle='none', color='purple')
axs[0].set_xlabel('temperatura^m (J^m)')
axs[0].set_ylabel('energia (s)')
axs[0].grid()


axs[1].plot(log_temperatura, log_energia, marker='s', linestyle='none', color='blue')
axs[1].set_xlabel('log(temperatura)')
axs[1].set_ylabel('log(energia)')
axs[1].grid()


axs[0].set_title("Ex3 - Potencia de Corpo Negro")

#------------  Colocar caixas de texto e gravar em PNG  ------------#


# place a text box in upper right in axes coords
textstr = f'$energia = {np.exp(b):.3f} \cdot exp_t({m:.3f})$ \n'
props = dict(boxstyle='round', facecolor='wheat', alpha=.5)
axs[0].text(min(temperatura), max(energia)+100, textstr,
        verticalalignment='top', bbox=props)

textstr = f'$ln(energia) = ln({np.exp(b):.3f}) +{m:.3f} \cdot ln(t) $ \n'+f'$r^2={r2:.3f}$'
props = dict(boxstyle='round', facecolor='wheat', alpha=.5)
axs[1].text(min(log_temperatura), max(log_energia), textstr,
        verticalalignment='top', bbox=props)


plt.savefig(dir_figs+'Ex3_potencia_corpo_negro.png') # Save the graph in png
plt.show()

# Abrir um separador com o grafico
# Tools > Preferences > IPython console > Graphics > Graphics Backend > Backend: Automatic
