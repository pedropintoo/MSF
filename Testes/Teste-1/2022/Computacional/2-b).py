# -*- coding: utf-8 -*-
"""
Created on Sat Mar 18 10:55:41 2023

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
tf = 100
dt = 0.0001

g = 9.8
y0 = 800
v0 = 0            # m/s
vt = 60     # m/s
vt_aberto = 5

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
    if(t[i] == 10):
        D = g/(vt_aberto**2)  #vt = 5.0 m/s quando abre o paraquedas

# Max position

if(len(np.where(y_RA <= 0)[0]) >= 1):
    
    inters_solo = np.where(y_RA <= 0)[0][0]
    t_inters_solo = t[inters_solo]
    print("Objeto chegou ao solo em t =",t_inters_solo,"s. Com velociade =",-1*v_RA[inters_solo],"m/s")
    plt.plot(t_inters_solo, 0, 'ko', linewidth=1) # Reta



# Plotting

plt.plot(t, y_RA, '-k', linewidth=1) # Reta


plt.xlabel('t (s)')
plt.ylabel('y (m)')
plt.grid()

plt.legend(['interseção com o solo','posicao com RA'])


#%% Gravar em PNG 


plt.savefig(dir_figs+'2-b).png') # Save the graph in png

plt.show() 

# Abrir um separador com o grafico
# Tools > Preferences > IPython console > Graphics > Graphics Backend > Backend: Automatic

