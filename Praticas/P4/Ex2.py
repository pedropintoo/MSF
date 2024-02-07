# -*- coding: utf-8 -*-
"""
Created on Thu Mar  9 17:34:20 2023

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


#%% EX2: Bola verticalmente


t0 = 0
tf = 2
dt = 0.01

g = 9.8
y0 = 0
v0 = 10            # m/s
a0 = -g
vt = 100 * 1000/3600

D = g/(vt**2)

Nt = int(np.ceil((tf - t0) / dt) + 1)

t = np.linspace(t0, tf, Nt)
# t = np.arange(t0,tf+dt,dt)


y = np.zeros((Nt,)) # np.zeros((Nt,2))  faria com 2 colunas
y[0] = y0
v = np.zeros((Nt,))
v[0] = v0
a = np.zeros((Nt,))

y_RA = np.zeros((Nt,)) # np.zeros((Nt,2))  faria com 2 colunas
y_RA[0] = y0
v_RA = np.zeros((Nt,))
v_RA[0] = v0
a_RA = np.zeros((Nt,))


# Método de Euler
for i in range(Nt - 1):
    # Sem RA
    a[i + 1] = a0 # Porque é constante (acelaracao = g)
    v[i + 1] = v[i] + a[i]*dt 
    y[i + 1] = y[i] + v[i]*dt
    
    # Com RA
    a_RA[i] = -g-D*v_RA[i]**2  # coloca-se o '-g' porque é em funcao de y
    v_RA[i + 1] = v_RA[i] + a_RA[i]*dt 
    y_RA[i + 1] = y_RA[i] + v_RA[i]*dt
    

# Max position

if(len(y == np.max(y)) >= 0):
    y_max = np.max(y)
    inters_max = np.where(y == np.max(y))[0][0]
    t_inters_max = t[inters_max]
    print("Máximo do objeto:",y_max,"em t =",t_inters_max)

if(len(np.where(y<=0)[0]) > 1):
    inters_ground = np.where(y<=0)[0][1]
    t_inters_ground = t[inters_ground]
    print("Posição inicial em t =",t_inters_ground)



# Plotting

plt.plot(t, y, '-r', linewidth=1) # Reta
plt.plot(t, y_RA, '-k', linewidth=1) # Reta




plt.xlabel('t (s)')
plt.ylabel('y (m)')
plt.grid()

plt.legend(['posicao sem RA','posicao com RA'])

plt.title("Ex1 - Queda livre")



#%% Gravar em PNG 


plt.savefig(dir_figs+'Ex1_bola_verticalmente.png') # Save the graph in png

plt.show() 

# Abrir um separador com o grafico
# Tools > Preferences > IPython console > Graphics > Graphics Backend > Backend: Automatic

