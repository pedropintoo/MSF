# -*- coding: utf-8 -*-
"""
Created on Tue May  9 22:46:13 2023

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
    
#%% Aproximação retangular

# Aproximação retangular:
# I = dx * np.sum(f[0:N]) = dt * np.sum(f[0:N] * vx[0:N])

def integral_retangular(f, v, iInicial, iFinal, dt):
    N = int((iFinal - iInicial) / dt)
    
    if N < 1:
        raise ValueError('Número de subintervalos inválido')
        
    integral = dt * sum(f[0:N]*v[0:N])
    return integral

#%% Aproximação trapezoidal

# Aproximação trapezial:
# I = dx * (f[0] + f[N])/2 + np.sum(f[1:N]) = dt * (f[0]*vx[0] + f[N]*vx[N])/2 + np.sum(f[1:N] * v[1:N])     

def integral_trapezoidal(f, v, iInicial, iFinal, dt):
    N = int((iFinal - iInicial) / dt)
    
    if N < 1:
        raise ValueError('Número de subintervalos inválido')
        
    integral = dt * ((f[0]*v[0] + f[N]*v[N])/2 + sum(f[1:N]*v[1:N]))
    return integral
    
    
#%% Método de Euler-Cromer (Complexo com RA) - 2 dimensão - Trabalho da resistencia do ar + Energia mecanica

PI = np.pi
p_ar = 1.225  # kg/m3

t0 = 0
tf = 2
dt = 0.001

g = 9.8
x0 = 0
y0 = 0
v0 = 140 * 1000/3600  # 140 km/h --> m/s
ang = 7  # 7 graus
v0x = v0 * np.cos(np.radians(ang))
v0y = v0 * np.sin(np.radians(ang))


vt = 100 * 1000/3600  # 100 km/h --> m/s

m = 0.057  # kg


D = g/(vt**2)

Nt = int(np.ceil((tf - t0) / dt) + 1)

t = np.linspace(t0, tf, Nt)

y = np.zeros((Nt,)) # np.zeros((Nt,2))  faria com 2 colunas
y[0] = y0

x = np.zeros((Nt,)) # np.zeros((Nt,2))  faria com 2 colunas
x[0] = x0  

vx = np.zeros((Nt,))
vx[0] = v0x
vy = np.zeros((Nt,))
vy[0] = v0y

ax = np.zeros((Nt,))
ay = np.zeros((Nt,))


Ec = np.zeros((Nt,)) # Energia cinética
Ec[0] = (1/2) * m * np.sqrt(v0x**2 + v0y**2)**2

Ep = np.zeros((Nt,))# Energia potencial
Ep[0] = 0 # partiu do solo

Em = np.zeros((Nt,)) # Energia mecanica
Em[0] = Ec[0] + Ep[0]

FresX = np.zeros((Nt,)) # Forca de resistencia do ar segundo X

FresY = np.zeros((Nt,)) # Forca de resistencia do ar segundo Y




# Método de Euler Cromer
for i in range(Nt - 1):
        
    vv = np.sqrt(vx[i]**2 + vy[i]**2) #intensidade
    
    ax[i]=-D*vv*vx[i]
    ay[i]=-g-D*vv*vy[i]
    
    vx[i+1]=vx[i]+ax[i]*dt
    vy[i+1]=vy[i]+ay[i]*dt
    
    x[i+1]=x[i]+vx[i+1]*dt
    y[i+1]=y[i]+vy[i+1]*dt
    
    
    Ec[i+1] = (1/2) * m * vv**2
    
    Ep[i+1] = m * g * y[i]
    
    # Deveria ser constante pois só tem forcas conservativas
    Em[i+1] = Ec[i+1] + Ep[i+1]
    
    FresX[i] = -D*vv*vx[i+1]*m
    FresY[i] = -D*vv*vy[i+1]*m

    
inters_ground = np.where(y <= 0)[0][1]
alcance = x[inters_ground]
print("Alcance da bola: ",alcance)
    
iInicial = 0
iFinal = t[inters_ground]

# Calcular o trabalho da Forca de resistencia

# Aprox. trapezoidal

TxTra = integral_trapezoidal(FresX, vx, iInicial, iFinal, dt)
TyTra = integral_trapezoidal(FresY, vy, iInicial, iFinal, dt)

TrabalhoTrapezio = TxTra + TyTra

print("Trabalho da resistencia do ar (trapezios) --> ",TrabalhoTrapezio)

# Plotting
plt.plot(x, y, '-r', linewidth=1) # Reta


plt.xlabel('x (m)')
plt.ylabel('y (m)')
plt.grid()

#%% Gravar em PNG 


plt.savefig(dir_figs+'1-a.png') # Save the graph in png

plt.show() 

# Abrir um separador com o grafico
# Tools > Preferences > IPython console > Graphics > Graphics Backend > Backend: Automatic

