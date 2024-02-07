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



#%% EX4: Amostra do Isotopo Radioativo


tempo =     np.array([0,      5,     10,    15,    20,    25,    30,     35,     40,     45])     # X
atividade = np.array([9.676 , 6.355, 4.261, 2.729, 1.862, 1.184, 0.7680, 0.4883, 0.3461, 0.2119]) # Y

log_atividade = np.log(atividade)


fig, axs = plt.subplots(2, 1,figsize=(5,7), sharex=True) # sharex=True faria com que o x fosse partilhado

axs[0].plot(tempo, atividade, marker='s', linestyle='none', color='purple')
axs[0].set_ylabel('activity (mCi)')
axs[0].set_xlabel('time (s)')
axs[0].grid()


axs[1].plot(tempo, log_atividade, marker='s', linestyle='none', color='blue')
axs[1].set_ylabel('log(activity) (mCi)')
axs[1].set_xlabel('time (s)')
axs[1].grid()


m, b, r2, delta_m, delta_b = lin_reg(tempo, log_atividade, 1)
log_fit = m*tempo + b          # Fit for log function
axs[1].plot(tempo, log_fit,'g', linewidth=1)


normal_fit = np.exp(log_fit)   # Convert to normal function
axs[0].plot(tempo, normal_fit,'g', linewidth=1)


axs[0].set_title("Ex4 - Amostra do Isotopo Radioativo")


#------------  Colocar caixas de texto e gravar em PNG  ------------#

# place a text box in upper right in axes coords
textstr = f'$a = {np.exp(b):.3f} \cdot exp({m:.3f} \cdot t)$ \n'
props = dict(boxstyle='round', facecolor='wheat', alpha=.5)
axs[0].text(max(tempo)-22, max(atividade), textstr,
        verticalalignment='top', bbox=props)

textstr = f'$ln(a) = ln({np.exp(b):.3f}) + ({m:.3f}) \cdot t $ \n'+f'$r^2={r2:.3f}$'
props = dict(boxstyle='round', facecolor='wheat', alpha=.5)
axs[1].text(max(tempo)-26, max(log_atividade), textstr,
        verticalalignment='top', bbox=props)



plt.savefig(dir_figs+'Ex4_amostra_isotipo_reativo.png') # Save the graph in png
plt.show()

# Abrir um separador com o grafico
# Tools > Preferences > IPython console > Graphics > Graphics Backend > Backend: Automatic



