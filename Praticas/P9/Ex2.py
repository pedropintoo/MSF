# -*- coding: utf-8 -*-
"""
Created on Thu May  4 17:53:31 2023

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

#%% Função para regressão linear

def lin_reg(x, y, print_res=0):
    
    if len(x) != len(y):
        raise ValueError('ERROR: x and y must be the same length in lin_reg(x, y)!')
    N = len(x)
    if N < 3:
        raise ValueError('ERROR: N must be higher than 2')
    
    # summations
    s_x = np.sum(x)
    s_y = np.sum(y)
    s_xy = np.sum(x * y)
    s_x2 = np.sum(x ** 2)
    s_y2 = np.sum(y ** 2)
    
    # linear regression parameters
    m = (N * s_xy - s_x * s_y) / (N * s_x2 - s_x ** 2)
    b = (s_x2 * s_y - s_x * s_xy) / (N * s_x2 - s_x ** 2)
    
    # squared correlation coefficient
    r2 = ((N * s_xy - s_x * s_y)**2) / ((N * s_x2 - s_x ** 2)*(N * s_y2 - s_y ** 2))

    # errors
    delta_m = np.abs(m) * np.sqrt((1 / r2 - 1) / (N - 2))
    delta_b = delta_m * np.sqrt(s_x2 / N)
    
    if print_res == 1:    # For printing
        print(f"m = {m}")
        print(f"b = {b}")
        print(f"r² = {r2}")
        print(f"delta_m = {delta_m}")
        print(f"delta_b = {delta_b}")
        
    
    return m, b, r2, delta_m, delta_b

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

  
#%% Calculo do trabalho para cada delta utilizando Aproximação retangular e trapezial

    
iInicial = 0
iFinal = 0.4

N = iFinal/dt

TrabalhoExato = np.array([-4.9768522, -4.9768522, -4.9768522, -4.9768522, -4.9768522])

ValoresTrabalhoRetangulo = np.array([])
ValoresTrabalhoTrapezio = np.array([])

ErrosRetangulo = np.array([])
ErrosTrapezio = np.array([])

deltas = np.array([0.1, 0.01, 0.001, 0.0001, 0.00001])

for dt in deltas:
    
    ###########################################################################
    
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
        
    ###########################################################################

    # Calcular o trabalho da Forca de resistencia
    
    # Aprox. retangular
    TxRet = integral_retangular(FresX, vx, iInicial, iFinal, dt)
    TyRet = integral_retangular(FresY, vy, iInicial, iFinal, dt)
    
    TrabalhoRetangulo = TxRet + TyRet
    
    # Aprox. trapezoidal
    
    TxTra = integral_trapezoidal(FresX, vx, iInicial, iFinal, dt)
    TyTra = integral_trapezoidal(FresY, vy, iInicial, iFinal, dt)
    
    TrabalhoTrapezio = TxTra + TyTra
    
    # Guardar valores e erros associados
    
    ValoresTrabalhoRetangulo = np.append(ValoresTrabalhoRetangulo,TrabalhoRetangulo)
    ErrosRetangulo = np.append(ErrosRetangulo, np.abs(TrabalhoRetangulo - TrabalhoExato[0]))

    ValoresTrabalhoTrapezio = np.append(ValoresTrabalhoTrapezio, TrabalhoTrapezio)
    ErrosTrapezio = np.append(ErrosTrapezio, np.abs(TrabalhoTrapezio - TrabalhoExato[0]))





fig, axs = plt.subplots(2, 1,figsize=(5,7)) # sharex=True faria com que o x fosse partilhado


# Aproximações

axs[0].plot(deltas, TrabalhoExato, color='green')

axs[0].plot(deltas, ValoresTrabalhoRetangulo, marker='s', linestyle='none', color='purple')

axs[0].plot(deltas, ValoresTrabalhoTrapezio, marker='x', linestyle='none', color='red')


axs[0].legend({"Aproximação retangular", "Aproximação trapezoidal", "Valor exato"})
axs[0].set_ylim(-5.05 ,-4.75)
axs[0].set_xlabel('Trabalho (W)')
axs[0].set_ylabel('Deltas (s)')
axs[0].grid()

# Erro log-log

log_deltas = np.log(deltas)
log_errosTrapezio = np.log(ErrosTrapezio)
log_errosRetangulo = np.log(ErrosRetangulo)

plt.plot(log_deltas, log_errosTrapezio, marker='s', linestyle='none', color='purple')
plt.plot(log_deltas, log_errosRetangulo, marker='x', linestyle='none', color='red')

m, b, r2, delta_m, delta_b = lin_reg(log_deltas, log_errosTrapezio, 1)
log_fit_trapezio = m*log_deltas + b

print("Aprox. trapezoidal --> erro =  e^",np.round(b,2)," * dt^",np.round(m,2)," (proporcional a dt)",sep="")



m, b, r2, delta_m, delta_b = lin_reg(log_deltas, log_errosRetangulo, 0)
log_fit_retangulo = m*log_deltas + b

print("Aprox. retangular --> erro = e^",np.round(b,2)," * dt^",np.round(m,2)," (proporcional a dt)",sep="")

  

#axs[1].plot(log_deltas, log_fit_trapezio, color='purple')
#axs[1].plot(log_deltas, log_fit_retangulo, color='red')


axs[1].set_xlabel('log(delta)')
axs[1].set_ylabel('log(erro)')
axs[1].grid()




#%% Gravar em PNG 


plt.savefig(dir_figs+'2.png') # Save the graph in png

plt.show() 

# Abrir um separador com o grafico
# Tools > Preferences > IPython console > Graphics > Graphics Backend > Backend: Automatic



