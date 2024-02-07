# -*- coding: utf-8 -*-
"""
Created on Wed Jun 21 15:20:51 2023

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

def acelera(t,x,v):                     # Altera consoante o problema!
    F_x = - k * x - 3 * alpha * x**2            
    return  (F_x) / m

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

#%% Coeficientes de Fourier

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

#%% Oscilador, analise Energia Potencial



x = np.linspace(-5,5,100) # Alterar consoante o problema e o Xeq
k = 1.2
alpha = -0.01
Ep = 1/2 * k * x**2 + alpha * x**3

EP_MENOR = 3

print("xi =",x[np.where(Ep <= EP_MENOR)[0][0]])
print("xf =",x[np.where(Ep <= EP_MENOR)[0][-1]])


plt.plot(x,Ep)
plt.plot(x[np.where(Ep <= EP_MENOR)[0][0]], EP_MENOR, "-or")
plt.plot(x[np.where(Ep <= EP_MENOR)[0][-1]], EP_MENOR, "-or")
plt.ylim(0,10)   # Alterar consoante o problema
plt.xlabel("x [m]")
plt.ylabel("Energia potencial [J]")
plt.grid()


#%% (Calculo ERRO) Método de Runge-Kutta de 4ª ordem - Oscilador não Harmónico forçado

# Nota: Usa outra funcao que já tem o ciclo dentro!

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
    

# Exata

v_exata = vt * np.tanh(g * t / vt)

# Resposta

print("Runge-Kutta - v(2s) = ",v_rungeKutta[t==2])
print("Euler Cromer - v(2s) = ",v_eulerC[t==2])
print("Exata - v(2s) = ",v_exata[t==2])


# Diferença/erro entre Euler Cromer e Runge-Kutta na velocidade

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


  
#%% Oscilador (Simples) - Função Runge-Kutta 4ª ordem (x,v)



t0 = 0
tf = 100
dt = 0.001
 
k = 1
b = 0.05 
alpha = 0.25
F0 = 7.5
wf = 1

x0 = 3 
v0 = 0

m = 1

g = 9.8


Nt = int(np.ceil((tf - t0) / dt) + 1)

t = np.linspace(t0, tf, Nt)


x = np.zeros((Nt,)) # np.zeros((Nt,2))  faria com 2 colunas
x[0] = x0

v = np.zeros((Nt,))
v[0] = v0

a = np.zeros((Nt,))




# Método de Runge-Kutta de 4ª ordem
for i in range(Nt-1):
    x[i+1],v[i+1] = rk4(t[i], x[i],v[i],acelera,dt)

  



plt.plot(t, x, '-b', linewidth=2)

plt.show()

plt.plot(x, v, '-r', linewidth=2)


plt.xlabel('t (s)')
plt.ylabel('y (m)')


plt.show()
    
#%% Oscilador quartico NÃO HARMONICO (Muito Complexo) - Função Runge-Kutta 4ª ordem (x,v)



t0 = 0
tf = 800
dt = 0.001
 
k = 1
alpha = 1
b = 0.05
F0 = 7.5
wf = 1

m = 1

g = 9.8

x0 = 3
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
Ep[0] = 1/2 * k * x0**2 * ( 1 + alpha * x0**2)
Em = np.zeros((Nt,)) # Energia mecanica
Em[0] = Ec[0] + Ep[0]




# Método de Runge-Kutta de 4ª ordem
for i in range(Nt-1):
    x[i+1],v[i+1] = rk4(t[i], x[i],v[i],acelera,dt)
    
    # Energias

    Ec[i+1] = (1/2) * m * np.sqrt(v[i+1]**2)**2

    Ep[i+1] = 1/2 * k * x[i+1]**2 * ( 1 + alpha * x[i+1]**2)

    Em[i+1] = Ec[i+1] + Ep[i+1]


# Resposta


k = 0
EXTREMOS_COUNT = 7
i = np.where(t >= 700)[0][0]
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
    
    
    if k == EXTREMOS_COUNT:
        break
    
    i += 1

      

        
# Amplitude = (Xmax0 - Xmin0) / 2
amplitude = (np.abs(np.max(x_ext)) + np.abs(np.min(x_ext))) / 2
print("Amplitude:",amplitude)

print("Max:",np.max(x_ext),"m")
print("Min:",np.min(x_ext),"m")

# Periodo = (tMax1 - tMax0)
# Frequencia = 1/Periodo
periodo = np.abs(t_ext[-1] - t_ext[0])
frequencia = 1 / periodo
print("Periodo:",periodo,"s") 
#print("Frequencia:",frequencia,"Hz") 
frequeciaAngular = 2 * np.pi * frequencia
#print("Frequencia angular:",frequeciaAngular,"rad/s")


# Coeficientes de Fourier

it0 = i_ext[0]
it1 = i_ext[2]

a = []
b = []
analise_Fourier = np.zeros(20)

for nf in range(20):
    ai, bi = abfourier(t,x,it0,it1,nf)
    a.append(ai)
    b.append(bi)
    analise_Fourier[nf] = np.sqrt(ai**2 + bi**2)
    


# Plotting

# Posicao

#plt.plot(t, x, '-b', linewidth=2)
plt.plot(t[it0:], x[it0:], '-b', linewidth=2)

plt.xlabel('t (s)')
plt.ylabel('x (m)')
plt.title(f'Lei do movimento x0 = {x0} e v0= {v0} ')

plt.grid()

#plt.legend(["Runge-Kutta de 4ª ordem"])

plt.savefig(dir_figs+'ex20-b)-1.png') # Save the graph in png

plt.show()

# Energia

plt.plot(t[it0:], Ec[it0:], '-b', linewidth=2)
plt.plot(t[it0:], Ep[it0:], '-r', linewidth=2)
plt.plot(t[it0:], Em[it0:], '-g', linewidth=4)

plt.xlabel('t (s)')
plt.ylabel('E (J)')

plt.title('Em = Ec + Ep')

plt.legend(["Ec","Ep","Em"])

plt.grid()

plt.savefig(dir_figs+'ex1-e)-1.png') # Save the graph in png

plt.show()

# Grafico de Fase

plt.plot(x[it0:],v[it0:], '-b', linewidth=2)
#plt.plot(x, v, '-b', linewidth=2)


plt.title('Grafico de Fase')


plt.savefig(dir_figs+'ex20-b)-2.png') # Save the graph in png

plt.show()

# Coeficientes

plt.plot(np.abs(a), 'ob', linewidth=2)
plt.plot(np.abs(b), 'or', linewidth=2)
plt.plot(analise_Fourier, '-g', linewidth=2)
plt.ylim(-0.25)

plt.bar(np.arange(len(analise_Fourier)), analise_Fourier, alpha=0.5)


plt.xlabel('n')

plt.title('analise_Fourier = $\sqrt{a_n^2 + b_n^2}$')

plt.legend(["$|an|$","$|bn|$","analise_Fourier","analise_Fourier"])

plt.grid()  
  

#%% Oscilador quartico (Complexo) - Função Runge-Kutta 4ª ordem (x,v)



t0 = 0
tf = 20
dt = 0.001
 
k = 1.2
alpha = -0.01

m = 1.5

g = 9.8

x0 = 3.5
v0 = 2

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
Ep[0] = 1/2 * k * x0**2 + alpha * x0**3
Em = np.zeros((Nt,)) # Energia mecanica
Em[0] = Ec[0] + Ep[0]




# Método de Runge-Kutta de 4ª ordem
for i in range(Nt-1):
    x[i+1],v[i+1] = rk4(t[i], x[i],v[i],acelera,dt)
    
    # Energias

    Ec[i+1] = (1/2) * m * np.sqrt(v[i+1]**2)**2

    Ep[i+1] = 1/2 * k * x[i+1]**2 + alpha * x[i+1]**3

    Em[i+1] = Ec[i+1] + Ep[i+1]


# Resposta

print("A Energia mecanica (neste caso) é constante: Em =",Em[-1],"J")


k = 0
i = np.where(t == 5)[0][0]
x_ext = []
t_ext = []
i_ext = []    # Indice para os coef. Fourier
EXTREMOS_COUNT = 4

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
    
    
    if k == EXTREMOS_COUNT:
        break
    
    i += 1

      

        
# Amplitude = (Xmax0 - Xmin0) / 2
# amplitude = (np.abs(x_ext[0]) + np.abs(x_ext[1])) / 2
# print("Amplitude:",amplitude)

print("Max:",np.max(x_ext),"m")
print("Min:",np.min(x_ext),"m")

# Periodo = (tMax1 - tMax0)
# Frequencia = 1/Periodo
periodo = np.abs(t_ext[2] - t_ext[0])
frequencia = 1 / periodo
print("Periodo:",periodo,"s") 
print("Frequencia:",frequencia,"Hz") 
frequeciaAngular = 2 * np.pi * frequencia
print("Frequencia angular:",frequeciaAngular,"rad/s")


# Coeficientes de Fourier

it0 = i_ext[0]
it1 = i_ext[2]

a = []
b = []
analise_Fourier = np.zeros(20)

for nf in range(20):
    ai, bi = abfourier(t,x,it0,it1,nf)
    a.append(ai)
    b.append(bi)
    analise_Fourier[nf] = np.sqrt(ai**2 + bi**2)
    


# Plotting

# Posicao

plt.plot(t, x, '-b', linewidth=2)
#plt.plot(t[it0:], x[it0:], '-b', linewidth=2)


plt.xlabel('t (s)')
plt.ylabel('x (m)')
plt.title(f'Lei do movimento x0 = {x0} e v0= {v0} ')

plt.grid()

#plt.legend(["Runge-Kutta de 4ª ordem"])

plt.savefig(dir_figs+'ex1-b)-1.png') # Save the graph in png

plt.show()

# Energia

plt.plot(t, Ec, '-b', linewidth=2)
plt.plot(t, Ep, '-r', linewidth=2)
plt.plot(t, Em, '-g', linewidth=4)

plt.xlabel('t (s)')
plt.ylabel('E (J)')

plt.title('Em = Ec + Ep')

plt.legend(["Ec","Ep","Em"])

plt.grid()

plt.savefig(dir_figs+'ex1-b)-2.png') # Save the graph in png

plt.show()


# Coeficientes

plt.plot(np.abs(a), 'ob', linewidth=2)
plt.plot(np.abs(b), 'or', linewidth=2)
plt.plot(analise_Fourier, '-g', linewidth=2)
plt.ylim(-0.25)

plt.bar(np.arange(len(analise_Fourier)), analise_Fourier, alpha=0.5)


plt.xlabel('n')

plt.title('analise_Fourier = $\sqrt{a_n^2 + b_n^2}$')

plt.legend(["$|an|$","$|bn|$","analise_Fourier","analise_Fourier"])

plt.grid()


#%% Gravar em PNG 


plt.savefig(dir_figs+'esteroides_computacionais.png') # Save the graph in png

plt.show() 

# Abrir um separador com o grafico
# Tools > Preferences > IPython console > Graphics > Graphics Backend > Backend: Automatic



