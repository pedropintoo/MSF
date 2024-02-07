# -*- coding: utf-8 -*-
"""
Created on Wed Jun 21 11:23:47 2023

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
    

#%% Função Runge-Kutta 4ª ordem (x,v)

def acelera(t,x,v):
    F_x = - k * x 
    F_res = -b * v
    F_ext = F0 * np.cos(wf * t)
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

#%% Oscilador Harmonico - Função Runge-Kutta 4ª ordem (x,v)



t0 = 0
tf = 400
dt = 0.001
 
m = 1
xeq = 0


k = 1
b = 0.05 
F0 = 7.5
wf = 1.4

g = 9.8

x0 = -3
v0 = 0

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
Ep[0] = 1/2 * k * x0**2
Em = np.zeros((Nt,)) # Energia mecanica
Em[0] = Ec[0] + Ep[0]


# Método de Runge-Kutta de 4ª ordem
for i in range(Nt-1):
    x[i+1],v[i+1] = rk4(t[i], x[i],v[i],acelera,dt)
    
    # Energias

    Ec[i+1] = (1/2) * m * np.sqrt(v[i+1]**2)**2

    Ep[i+1] = 1/2 * k * x[i+1]**2

    Em[i+1] = Ec[i+1] + Ep[i+1]
    


# Resposta


k = 0
i = np.where(t == 350)[0][0]     # Depende de quando começa o regime estacionário
x_ext = []
t_ext = []
i_ext = []    # Indice para os coef. Fourier

while True:
    if(x[i] > x[i-1] and x[i] > x[i+1]):
        tmax, xmax = getMaxMin(t[i-1], t[i], t[i+1], x[i-1], x[i], x[i+1])
        t_ext.append(tmax)
        x_ext.append(xmax)
        
        i_ext.append(i)
        k += 1
       
    if(x[i] < x[i-1] and x[i] < x[i+1]):
        tmin, xmin = getMaxMin(t[i-1], t[i], t[i+1], x[i-1], x[i], x[i+1])
        
        t_ext.append(tmin)
        x_ext.append(xmin)
        
        i_ext.append(i)  # Indice para os coef. Fourier
        k += 1
    
    
    if k == 4:
        break
    
    i += 1

      

        
# Amplitude = (Xmax0 - Xmin0) / 2
amplitude = (np.abs(x_ext[0]) + np.abs(x_ext[1])) / 2
print("Amplitude:",amplitude)

#print("Max:",np.round(np.max(x_ext),3))
#print("Min:",np.round(np.min(x_ext),3))

# Periodo = (tMax1 - tMax0)
# Frequencia = 1/Periodo
#periodo = np.abs(t_ext[2] - t_ext[0])
#frequencia = 1 / periodo
#print("Periodo:",np.round(periodo,3)) 
#print("Frequencia:",np.round(frequencia,3)) 
#frequeciaAngular = 2 * np.pi * frequencia
#print("Frequencia angular:",np.round(frequeciaAngular,3))



# Plotting

# Energia

plt.plot(t, Ec, '-b', linewidth=1)
plt.plot(t, Ep, '-r', linewidth=1)
plt.plot(t, Em, '-g', linewidth=2)

plt.xlabel('t (s)')
plt.ylabel('E (J)')

plt.title('Em = Ec + Ep')

plt.legend(["Ec","Ep","Em"])

plt.grid()
plt.show()


# Posicao

plt.plot(t, x, '-b', linewidth=2)


plt.xlabel('t (s)')
plt.ylabel('x (m)')
plt.title(f'Lei do movimento x0 = {x0} e v0= {v0} ')

plt.grid()

#plt.legend(["Runge-Kutta de 4ª ordem"])


#%% Gravar em PNG 


plt.savefig(dir_figs+'ex2-a).png') # Save the graph in png

plt.show() 

# Abrir um separador com o grafico
# Tools > Preferences > IPython console > Graphics > Graphics Backend > Backend: Automatic

