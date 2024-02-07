# -*- coding: utf-8 -*-
"""
Created on Thu May 25 16:36:23 2023

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

def acel(t,g,vt,v):
    return g - (g / vt **2) * np.abs(v) * v


def rk4(t,x,v,g,vt):
    N = len(t)
    dt = t[1]-t[0]
    
    for i in range(N-1):
        a1=acel(t[i],g,vt,v[i])
        c1v=a1*dt
        c1x=v[i]*dt
        a2=acel(t[i]+dt/2,g,vt,v[i]+c1v/2.)
        c2v=a2*dt
        c2x=(v[i]+c1v/2.)*dt # predicto: v(t+dt) * dt
        a3=acel(t[i]+dt/2,g,vt,v[i]+c2v/2.)
        c3v=a3*dt
        c3x=(v[i]+c2v/2.)*dt
        a4=acel(t[i]+dt,g,vt,v[i]+c3v)
        c4v=a4*dt
        c4x=(v[i]+c3v)*dt
        x[i+1]=x[i]+(c1x+2.*c2x+2.*c3x+c4x)/6.
        v[i+1]=v[i]+(c1v+2.*c2v+2.*c3v+c4v)/6.
    return x,v



#%% Método de Runge-Kutta de 4ª ordem - Oscilador não Harmónico forçado



t0 = 0
tf = 3
dt = 0.01
 
vt = 6.8

g = 9.8


Nt = int(np.ceil((tf - t0) / dt) + 1)

t = np.linspace(t0, tf, Nt)


y_rungeKutta = np.zeros((Nt,)) # np.zeros((Nt,2))  faria com 2 colunas


v_rungeKutta = np.zeros((Nt,))


a_rungeKutta = np.zeros((Nt,))


# Defini a formula da acelaracao em acel(...)
# Método de Runge-Kutta de 4ª ordem
y_rungeKutta,v_rungeKutta = rk4(t,y_rungeKutta,v_rungeKutta,g,vt)




y_eulerC = np.zeros((Nt,)) # np.zeros((Nt,2))  faria com 2 colunas


v_eulerC = np.zeros((Nt,))


a_eulerC = np.zeros((Nt,))

# Método de Euler Cromer
for i in range(Nt - 1):

    v_norma = np.sqrt(v_eulerC[i]**2)
    
    a_eulerC[i]= g - (g * v_norma * v_eulerC[i] ) / (vt**2)

    
    v_eulerC[i+1]=v_eulerC[i]+a_eulerC[i]*dt
    
    
    y_eulerC[i+1]=y_eulerC[i]+v_eulerC[i+1]*dt
    

v_exata = vt * np.tanh(g * t / vt)

# Resposta

print("Runge-Kutta - v(2s) = ",v_rungeKutta[t==2])
print("Euler Cromer - v(2s) = ",v_eulerC[t==2])
print("Exata - v(2s) = ",v_exata[t==2])


#%% Diferença/erro entre Euler Cromer e Runge-Kutta na velocidade

erroV_rungeKutta = np.zeros((Nt,))
erroV_eulerC = np.zeros((Nt,))

for i in range(Nt - 1):
    
    erroV_rungeKutta[i] = np.abs(v_rungeKutta[i] -  v_exata[i])
    
    erroV_eulerC[i] = np.abs(v_eulerC[i] -  v_exata[i])
    

# Plotting

plt.plot(t, v_rungeKutta, '-b', linewidth=6)
plt.plot(t, v_eulerC, '-r', linewidth=2)

plt.xlabel('t (s)')
plt.ylabel('y (m)')

plt.legend(["Runge-Kutta de 4ª ordem","Euler Cromer"])

plt.show()

fig, axs = plt.subplots(2, 1,figsize=(5,7), sharex=True) # sharex=True faria com que o x fosse partilhado

axs[0].plot(t, erroV_rungeKutta, '-b', linewidth=2)
axs[1].plot(t, erroV_eulerC, '-r', linewidth=2)

axs[0].legend(["Erro - Runge-Kutta de 4ª ordem"])
axs[1].legend(["Erro - Euler Cromer"])



#%% Gravar em PNG 


plt.savefig(dir_figs+'ex2.png') # Save the graph in png

plt.show() 

# Abrir um separador com o grafico
# Tools > Preferences > IPython console > Graphics > Graphics Backend > Backend: Automatic


