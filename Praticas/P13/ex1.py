# -*- coding: utf-8 -*-
"""
Created on Thu Jun  1 16:22:49 2023

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
    F_amort = -b * v
    F_ext = F0 * np.cos(wf * t)
    F_x = - k * x * (1 + 2 * alpha * x ** 2)
    return  (F_amort + F_ext + F_x) / m

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

def abfourier(tp,xp,it0,it1,nf):
    # cálculo dos coeficientes de Fourier a_nf e b_nf
    # a_nf = 2/T integral ( xp cos( nf w) ) dt entre tp(it0) e tp(it1)
    # b_nf = 2/T integral ( xp sin( nf w) ) dt entre tp(it0) e tp(it1)
    # integracao numerica pela aproximação trapezoidal
    # input: matrizes tempo tp (abcissas)
    # posição xp (ordenadas)
    # indices inicial it0
    # final it1 (ao fim de um período)
    # nf índice de Fourier
    # output: af_bf e bf_nf
    dt=tp[1]-tp[0]
    per=tp[it1]-tp[it0]
    ome=2*np.pi/per
    s1=xp[it0]*np.cos(nf*ome*tp[it0])
    s2=xp[it1]*np.cos(nf*ome*tp[it1])
    st=xp[it0+1:it1]*np.cos(nf*ome*tp[it0+1:it1])
    soma=np.sum(st)
    q1=xp[it0]*np.sin(nf*ome*tp[it0])
    q2=xp[it1]*np.sin(nf*ome*tp[it1])
    qt=xp[it0+1:it1]*np.sin(nf*ome*tp[it0+1:it1])
    somq=np.sum(qt)
    intega=((s1+s2)/2+soma)*dt
    af=2/per*intega
    integq=((q1+q2)/2+somq)*dt
    bf=2/per*integq
    return af,bf

#%% Oscilador quartico - Função Runge-Kutta 4ª ordem (x,v)



t0 = 0
tf = 300
dt = 0.01
 
k = 1
b = 0.05 
alpha = 0.002
F0 = 7.5
wf = 1

m = 1

g = 9.8


Nt = int(np.ceil((tf - t0) / dt) + 1)

t = np.linspace(t0, tf, Nt)


x = np.zeros((Nt,)) # np.zeros((Nt,2))  faria com 2 colunas
x[0] = 3  

v = np.zeros((Nt,))
v[0] = 0

a = np.zeros((Nt,))




# Método de Runge-Kutta de 4ª ordem
for i in range(Nt-1):
    x[i+1],v[i+1] = rk4(t[i], x[i],v[i],acelera,dt)



# Resposta

k = 0
i = np.where(t == 150)[0][0]
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
    
    
    if k == 3:
        break
    
    i += 1

      

# Alinea b
        
# Amplitude = (Xmax0 - Xmin0) / 2
amplitude = (np.abs(x_ext[0]) + np.abs(x_ext[1])) / 2
print("Amplitude:",amplitude)

# Periodo = (tMax1 - tMax0)
# Frequencia = 1/Periodo
periodo = np.abs(t_ext[2] - t_ext[0])
#frequencia = 1 / periodo
print("Periodo:",periodo) 
#print("Frequencia:",frequencia) 


# Alinea c) - Coeficientes de Fourier

it0 = i_ext[0]
it1 = i_ext[2]

a = []
b = []
x_f = np.zeros(len(t[it0:]))

for nf in range(40):
    
    ai, bi = abfourier(t,x,it0,it1,nf)
    a.append(ai)
    b.append(bi)
    x_f += ai * np.cos(nf * t[it0:]) + bi * np.sin(nf * t[it0:])
    


# Plotting
plt.plot(a, 'ob', linewidth=2)
plt.plot(b, 'or', linewidth=2)

plt.show()

plt.plot(t, x, '-b', linewidth=2)

#plt.plot(t[it0:], x[it0:], '-b', linewidth=2)

#plt.plot(t[it0:], x_f, '-r', linewidth=2)


plt.xlabel('t (s)')
plt.ylabel('y (m)')

plt.legend(["Runge-Kutta de 4ª ordem", "Série de Fourier"])



#%% Gravar em PNG 


plt.savefig(dir_figs+'ex1.png') # Save the graph in png

plt.show() 

# Abrir um separador com o grafico
# Tools > Preferences > IPython console > Graphics > Graphics Backend > Backend: Automatic



