# -*- coding: utf-8 -*-
"""
Created on Tue Mar 21 22:42:56 2023

@author: pedro
"""

#%% Importar libraries e diretorio para figuras

import matplotlib.pyplot as plt
import numpy as np
import os

dir_figs = './figures/'
if not os.path.exists(dir_figs): 
    os.mkdir(dir_figs)


#%% Queda livre com RA (método de Euler)


t0 = 0
tf = 1
dt = 0.001

g = 9.8
y0 = 0      # Starting in y0 = 500 m
v0 = -200 * 1000 / 3600            # -200 km/h --> m/s
vt = 6.80    # m/s

D = g/(vt**2)

Nt = int(np.ceil((tf - t0) / dt) + 1)

t = np.linspace(t0, tf, Nt)
# t = np.arange(t0,tf+dt,dt)


y_RA = np.zeros((Nt,)) # np.zeros((Nt,2))  faria com 2 colunas
y_RA[0] = y0
v_RA = np.zeros((Nt,))
v_RA[0] = v0
a_RA = np.zeros((Nt,))


# Método de Euler
for i in range(Nt - 1):

    # Com RA
    vv = np.sqrt(v_RA[i]**2 + 0) # Se quisesse fazer um vx e vy teria de trocar o 0 --> vx_RA
    a_RA[i] = -g-D*vv*v_RA[i]  # coloca-se o '-g' porque é em funcao de y
    v_RA[i + 1] = v_RA[i] + a_RA[i]*dt 
    y_RA[i + 1] = y_RA[i] + v_RA[i]*dt
    

# Max position

# if(len(y == np.max(y)) >= 0):
#     y_max = np.max(y)
#     inters_max = np.where(y == np.max(y))[0][0]
#     t_inters_max = t[inters_max]
#     print("Máximo do objeto:",y_max,"em t =",t_inters_max)

# if(len(np.where(y<=0)[0]) > 1):
#     inters_ground = np.where(y<=0)[0][1]
#     t_inters_ground = t[inters_ground]
#     print("Posição inicial em t =",t_inters_ground)

inters_ground = np.where(y_RA<=-4)[0][0]
t_inters_ground = t[inters_ground]
print("t =",t_inters_ground)



# Plotting

plt.plot(t, y_RA, '-k', linewidth=1) # Reta






#%% Gravar em PNG 


plt.savefig(dir_figs+'2-c).png') # Save the graph in png

plt.show() 

# Abrir um separador com o grafico
# Tools > Preferences > IPython console > Graphics > Graphics Backend > Backend: Automatic
