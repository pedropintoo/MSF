# -*- coding: utf-8 -*-
"""
Created on Thu Mar  9 16:14:06 2023

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


#%% EX1: Queda livre


t0 = 0
tf = 4
dt = 0.01

g = 9.8
a0 = g              # m/s^2

Nt = int(np.ceil((tf - t0) / dt) + 1)

t = np.linspace(t0, tf, Nt)
# t = np.arange(t0,tf+dt,dt)


y = np.zeros((Nt,)) # np.zeros((Nt,2))  faria com 2 colunas
v = np.zeros((Nt,))
a = np.zeros((Nt,))


# Método de Euler
for i in range(Nt - 1):
    a[i + 1] = a0 # Porque é constante (acelaracao = g)
    v[i + 1] = v[i] + a[i]*dt 
    y[i + 1] = y[i] + v[i]*dt


# Plotting
fig, axs = plt.subplots(2, 1,figsize=(5,7)) # sharex=True faria com que o x fosse partilhado


axs[0].plot(t, y, '-r', linewidth=1) # Reta
axs[0].plot(t, v, '-b', linewidth=1) # Reta



axs[0].set_xlabel('t (s)')
axs[0].set_ylabel('v (m/s)')
axs[0].grid()

axs[0].legend(['posicao','velocidade'])

axs[0].set_title("Ex1 - Queda livre")

#%% Erro Aproximado consoante os passos relativamente á posição

y_exato = (1/2)*g*t**2 # y = y0 + v0*t + (1/2)*g*t**2

axs[1].plot( y_exato - y, 'k-', linewidth=1)

axs[1].set_xlabel('Passos')
axs[1].set_ylabel('Erro')
axs[1].grid()

axs[1].legend(['Erro da Posição'])



#%% Pontos de interseção (t == 3 e t == 2)

# t == 3

v_sec3 = g * 3
print("Quando t = 3 a velocidade é igual a",v_sec3,"m/s")
print("Passo = 0.01 --> v =",v[t==3][0],"m/s")
print("Erro:",(v_sec3-v[t==3][0]))
print("Passo = 0.10 --> v =",28.42,"m/s") # 28.42 é o valor que dá com o passo = 0.10
print("Erro:",(v_sec3-28.42))
print("Ao multiplicar por 10 o delta t, o erro multipla por 10!")
axs[0].plot(3,v[t==3][0],'bo')

print("\n")
# t == 2

y_sec2 = (1/2) * g * 2**2
print("Quando t = 2 a posicao é igual a",y_sec2,"m")
print("Passo = 0.01 --> x =",y[t==2][0],"m")
print("Erro:",(y_sec2-y[t==2][0]))
print("Passo = 0.10 --> x =",16.758,"m") # 28.42 é o valor que dá com o passo = 0.10
print("Erro:",(y_sec2-16.758))
print("Ao multiplicar por 10 o delta t, o erro multipla por 10!")
axs[0].plot(2,y[t==2][0],'ro')





#------------  Gravar em PNG  ------------#


plt.savefig(dir_figs+'Ex1_queda_livre.png') # Save the graph in png

plt.show() 

# Abrir um separador com o grafico
# Tools > Preferences > IPython console > Graphics > Graphics Backend > Backend: Automatic


