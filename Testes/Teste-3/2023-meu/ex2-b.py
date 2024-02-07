# -*- coding: utf-8 -*-
"""
Created on Fri Jun 23 10:06:45 2023

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
    


#%% Função Runge-Kutta 4ª ordem (x,v)

def acelera(t,x,v):                     # Altera consoante o problema!
    F_x = - 4 * alpha * x**3           
    F_res = -b * v
    F_ext = F0 * np.cos(wf*t)
    return  (F_x + F_res + F_ext) / m

def rk4(t,x,vx,acelera,dt):
    """
    Integração numérica de equação diferencial de 2ª ordem:
    d2x/dt2 = ax(t,x,vx)    com dx/dt= vx    de valor inicial
    Erro global:  proporcional a dt**4
    acelera=dvx/dt=Força(t,x,vx)/massa      : acelera é uma FUNÇÃO
    input:  t = instante de tempo
            x(t) = posição
            vx(t) = velocidade
            dt = passo temporal
    output: xp = x(t+dt)
    vxp = vx(t+dt)
    """
    ax1=acelera(t,x,vx)
    c1v=ax1*dt
    c1x=vx*dt
    ax2=acelera(t+dt/2.,x+c1x/2.,vx+c1v/2.)
    c2v=ax2*dt
    c2x=(vx+c1v/2.)*dt # predicto:  vx(t+dt) * dt
    ax3=acelera(t+dt/2.,x+c2x/2.,vx+c2v/2.)
    c3v=ax3*dt
    c3x=(vx+c2v/2.)*dt
    ax4=acelera(t+dt,x+c3x,vx+c3v)
    c4v=ax4*dt
    c4x=(vx+c3v)*dt
     
    xp=x+(c1x+2.*c2x+2.*c3x+c4x)/6.
    vxp=vx+(c1v+2.*c2v+2.*c3v+c4v)/6.
    return xp,vxp


#%% Oscilador quartico (Complexo) - Função Runge-Kutta 4ª ordem (x,v)



t0 = 0
tf = 80
dt = 0.001
 
k = 1
alpha = 0.15
b = 0.02
F0 = 7.5
wf = 1.0

m = 1

g = 9.8

x0_1 = 1.999
v0_1 = 0

x0_2 = 2.000
v0_2 = 0

x0_3 = 2.001
v0_3 = 0

Nt = int(np.ceil((tf - t0) / dt) + 1)

t = np.linspace(t0, tf, Nt)


x_1 = np.zeros((Nt,)) # np.zeros((Nt,2))  faria com 2 colunas
x_1[0] = x0_1
v_1 = np.zeros((Nt,))
v_1[0] = v0_1
a_1 = np.zeros((Nt,))

x_2 = np.zeros((Nt,)) # np.zeros((Nt,2))  faria com 2 colunas
x_2[0] = x0_2
v_2 = np.zeros((Nt,))
v_2[0] = v0_2
a_2 = np.zeros((Nt,))

x_3 = np.zeros((Nt,)) # np.zeros((Nt,2))  faria com 2 colunas
x_3[0] = x0_3
v_3 = np.zeros((Nt,))
v_3[0] = v0_3
a_3 = np.zeros((Nt,))




# Método de Runge-Kutta de 4ª ordem
for i in range(Nt-1):
    #x0_1 = 1.999
    x_1[i+1],v_1[i+1] = rk4(t[i], x_1[i],v_1[i],acelera,dt)
    
    #x0_2 = 2.000
    x_2[i+1],v_2[i+1] = rk4(t[i], x_2[i],v_2[i],acelera,dt)
    
    #x0_3 = 2.001
    x_3[i+1],v_3[i+1] = rk4(t[i], x_3[i],v_3[i],acelera,dt)
    


# Plotting

# Posicao

plt.plot(t, x_1, '-b', linewidth=2)
plt.plot(t, x_2, '-g', linewidth=2)
plt.plot(t, x_3, '-r', linewidth=2)


plt.xlabel('t (s)')
plt.ylabel('x (m)')
plt.title("Lei do Movimento - Runge-Kutta de 4ª ordem")

plt.grid()

plt.legend([f'x0 = {x0_1} e v0= {v0_1} ',
            f'x0 = {x0_2} e v0= {v0_2} ',
            f'x0 = {x0_3} e v0= {v0_3} '])

plt.savefig(dir_figs+'ex2-b)-1.png') # Save the graph in png

plt.show() 


# Analisando mais de Perto
t_estudo = np.where(t >= 50)[0][0]

plt.plot(t[t_estudo:], x_1[t_estudo:], '-b', linewidth=2)
plt.plot(t[t_estudo:], x_2[t_estudo:], '-g', linewidth=2)
plt.plot(t[t_estudo:], x_3[t_estudo:], '-r', linewidth=2)


plt.xlabel('t (s)')
plt.ylabel('x (m)')
plt.title("Lei do Movimento - Runge-Kutta de 4ª ordem")

plt.grid()

plt.legend([f'x0 = {x0_1} e v0= {v0_1} ',
            f'x0 = {x0_2} e v0= {v0_2} ',
            f'x0 = {x0_3} e v0= {v0_3} '])

#Resposta:
    
# Caraterísticas:
# - Soluções com condições iniciais extramente próximas DIVERGEM!
# - Apesar de a solução ser ÚNICA (para uma condição inicial), 
# não se consegue calcular (prever) a evolução a médio prazo (CAOS)

print("Resposta: até ~59s")
# até ~59s


#%% Gravar em PNG 


plt.savefig(dir_figs+'ex2-b)-2.png') # Save the graph in png

plt.show() 

# Abrir um separador com o grafico
# Tools > Preferences > IPython console > Graphics > Graphics Backend > Backend: Automatic

