# -*- coding: utf-8 -*-
"""
Created on Thu May 18 17:07:29 2023

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
    
#%% Maximo e minimo (3 pontos)
  
def getMaxMin(x0,x1,x2,y0,y1,y2):
    # Máximo ou mínimo usando o polinómio de Lagrange
    # Dados (input): (x0,y0), (x1,y1) e (x2,y2)
    # Resultados (output): xm, yMaxMin
    
    xab=x0-x1
    xac=x0-x2
    xbc=x1-x2
    
    a=y0/(xab*xac)
    b=-y1/(xab*xbc)
    c=y2/(xac*xbc)
    
    xmla=(b+c)*x0+(a+c)*x1+(a+b)*x2
    xm=0.5*xmla/(a+b+c)
    
    xta=xm-x0
    xtb=xm-x1
    xtc=xm-x2
    
    yMaxMin=a*xtb*xtc+b*xta*xtc+c*xta*xtb
    
    return xm, yMaxMin       
    

#%% Método de Euler Cromer (Complexo com RA) - Oscilador Harmónico Duplo



t0 = 0
tf = 15
dt = 0.05

xeq = 1.5
Energia_total = 3

m = 1  # kg
k = 1   # Constante elástica

g = 9.8
x0 = xeq
# Alterando a energia total estamos a alterar a velocidade inicial
v0 = np.sqrt(Energia_total * 2 / m)   # v0 = np.sqrt(Energia_total * 2 / m) 


Nt = int(np.ceil((tf - t0) / dt) + 1)

t = np.linspace(t0, tf, Nt)


x = np.zeros((Nt,)) # np.zeros((Nt,2))  faria com 2 colunas
x[0] = x0  

v = np.zeros((Nt,))
v[0] = v0

a = np.zeros((Nt,))


Ec = np.zeros((Nt,)) # Energia cinética
Ec[0] = (1/2) * m * np.sqrt(v0**2)**2

Ep = np.zeros((Nt,))# Energia potencial 
Ep[0] = (1/2) * k * (x0**2 - xeq**2)**2
Em = np.zeros((Nt,)) # Energia mecanica
Em[0] = Ec[0] + Ep[0]


# Método de Euler Cromer
for i in range(Nt - 1):


    a[i]= ( -2*k * (x[i]**2 - xeq**2)*x[i] )/m # Fx = -2*k * (x**2 - xeq**2)*x --> a = -2 * k * (x**2 - xeq**2)*x/m

    
    v[i+1]=v[i]+a[i]*dt
    
    
    x[i+1]=x[i]+v[i+1]*dt
    
    # Energias
    
    Ec[i+1] = (1/2) * m * np.sqrt(v[i+1]**2)**2
    
    Ep[i+1] = (1/2) * k * (x[i+1]**2 - xeq**2)**2
    
    # Deveria ser constante pois só tem forcas conservativas
    Em[i+1] = Ec[i+1] + Ep[i+1]
    


k = 0
EXTREMOS_COUNT = 4
i = np.where(t == 2)[0][0]
x_ext = []
t_ext = []

while True:
    if(x[i] > x[i-1] and x[i] > x[i+1]):
        tmax, xmax = getMaxMin(t[i-1], t[i], t[i+1], x[i-1], x[i], x[i+1])
        
        t_ext.append(tmax)
        x_ext.append(xmax)
        k += 1
       
    if(x[i] < x[i-1] and x[i] < x[i+1]):
        tmin, xmin = getMaxMin(t[i-1], t[i], t[i+1], x[i-1], x[i], x[i+1])   
        
        t_ext.append(tmin)
        x_ext.append(xmin)
        k += 1
    
    
    if k == EXTREMOS_COUNT:
        break
    
    i += 1

  
      
# Amplitude = (Xmax0 - Xmin0) / 2
amplitude = (np.abs(x_ext[0]) + np.abs(x_ext[1])) / 2
print("Amplitude:",amplitude)

# Periodo = (tMax1 - tMax0)
periodo = np.abs(t_ext[2] - t_ext[0])
#print("Periodo:",periodo) 

frequencia = 1 / periodo
print("Frequencia:",frequencia)    
      


# Plotting

plt.plot(t,x, '-b', linewidth=2)
#plt.plot(t,Em, '-b', linewidth=0.5)



plt.xlabel('t (s)')
plt.ylabel('x (m)')

plt.legend(["Euler Cromer","Energia mecanica"])


#%% Gravar em PNG 


plt.savefig(dir_figs+'ex2-b.png') # Save the graph in png

plt.show() 

# Abrir um separador com o grafico
# Tools > Preferences > IPython console > Graphics > Graphics Backend > Backend: Automatic

