# -*- coding: utf-8 -*-
"""
Created on Thu Mar 16 16:34:55 2023

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


#%% EX5: Representar vetores


# (x0,y0) = ponto inicial do vetor
# v1 = (x,y) = comprimentos x e y do vetor

PI = np.pi
F = 5
x0 = 0
y0 = 0


conv = PI/180 # converter graus em radianos G*conv = R
theta = np.array([PI/2, 60*conv, -7*PI/6, 310*conv])

vetores = np.array([F*np.cos(theta), F*np.sin(theta)])


fig, ax = plt.subplots()
ax.plot(x0,y0, 'o', markersize = 5)

for i in range(len(theta)):
    ax.arrow(x0, y0, vetores[0,i], vetores[1,i], color='r', width=0.05)
ax.set_aspect('equal')


# Confirmar com uma circunferencia de raio 5

# Creating equally spaced 100 data in range 0 to 2*pi
theta = np.linspace(0, 2 * np.pi, 100)

# Setting radius
radius = 5

# Generating x and y data
x = radius * np.cos(theta)
y = radius * np.sin(theta)

# Plotting
ax.plot(x, y)

plt.show()

#%% Gravar em PNG 


plt.savefig(dir_figs+'Ex5_vetor_representacao.png') # Save the graph in png

plt.show() 

# Abrir um separador com o grafico
# Tools > Preferences > IPython console > Graphics > Graphics Backend > Backend: Automatic
