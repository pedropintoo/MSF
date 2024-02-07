# -*- coding: utf-8 -*-
"""
Created on Thu May  4 16:58:08 2023

@author: pedro
"""


#%% Importar libraries e diretorio para figuras

import matplotlib.pyplot as plt
import numpy as np
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

#%% Método de Euler (Complexo com RA) - 2 dimensão

PI = np.pi
p_ar = 1.225  # kg/m3

t0 = 0
tf = 1
dt = 0.0001

g = 9.8
x0 = 0
y0 = 0
v0 = 100 * 1000/3600  # 130 km/h --> m/s
ang = 10  # 10 graus
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

    
    
iInicial = 0
iFinal = 0.4

# Calcular o trabalho da Forca de resistencia

# Aprox. retangular
TxRet = integral_retangular(FresX, vx, iInicial, iFinal, dt)
TyRet = integral_retangular(FresY, vy, iInicial, iFinal, dt)

TrabalhoRetangulo = TxRet + TyRet

# Aprox. trapezoidal

TxTra = integral_trapezoidal(FresX, vx, iInicial, iFinal, dt)
TyTra = integral_trapezoidal(FresY, vy, iInicial, iFinal, dt)

TrabalhoTrapezio = TxTra + TyTra



print("Trabalho da resistencia do ar (retangulos) --> ",TrabalhoRetangulo)

print("Trabalho da resistencia do ar (trapezios) --> ",TrabalhoTrapezio)



print()

print("Energia mecanica em t = 0 -->",Em[0])

print("Energia mecanica em t = 0.8 -->",Em[int(0.8*Nt)])


# Plotting
plt.plot(t, Em, '-r', linewidth=1) # Reta


plt.xlabel('t (s)')
plt.ylabel('Em (J)')
plt.grid()



#%% Gravar em PNG 


plt.savefig(dir_figs+'1-c).png') # Save the graph in png

plt.show() 

# Abrir um separador com o grafico
# Tools > Preferences > IPython console > Graphics > Graphics Backend > Backend: Automatic



