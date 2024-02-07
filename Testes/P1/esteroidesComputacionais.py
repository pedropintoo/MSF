# -*- coding: utf-8 -*-
"""
Created on Sat Mar 18 10:55:41 2023

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

#%% Regressão linear (normal) Método dos mínimos quadrados

tempo = np.array(     [    0,     1,     2,     3,     4,     5,     6,     7,     8,     9 ]) # X
distancia = np.array([ 0.00, 0.735, 1.363, 1.739, 2.805, 3.814, 4.458, 4.955, 5.666, 6.329 ])  # Y

plt.plot(tempo, distancia, marker='s', linestyle='none', color='purple')
plt.xlabel('tempo (min)')
plt.ylabel('distancia (km)')
plt.grid()


# Limite de Plot in X and Y
# plt.xlim(ax, bx) 
# plt.ylim(cy, dy)

plt.title("Ex2 - Ciclista")

# regressão linear
m, b, r2, delta_m, delta_b = lin_reg(tempo,distancia)
normal_fit = m*tempo + b

plt.plot(tempo, normal_fit, '-g', linewidth=1)


#%% Regressão linear (log-log) Método dos mínimos quadrados



temperatura = np.array([200,   300,   400,   500,   600,   700,   800,   900,  1000,  1100 ])   # X
energia = np.array([ 0.6950, 4.363, 15.53, 38.74, 75.08, 125.2, 257.9, 344.1, 557.4, 690.7])    # Y

log_temperatura = np.log(temperatura)
log_energia = np.log(energia)


fig, axs = plt.subplots(2, 1,figsize=(5,7)) # sharex=True faria com que o x fosse partilhado


m, b, r2, delta_m, delta_b = lin_reg(log_temperatura, log_energia, 1)
log_fit = m*log_temperatura + b

axs[1].plot(log_temperatura, log_fit,'g', linewidth=1)

# Converting log_fit --> normal_fit

temperatura = temperatura**m

m, b, r2, delta_m, delta_b = lin_reg(temperatura, energia, 1)
normal_fit = m*temperatura + b

axs[0].plot(temperatura, normal_fit,'g', linewidth=1)


axs[0].plot(temperatura, energia, marker='s', linestyle='none', color='purple')
axs[0].set_xlabel('temperatura^m (J^m)')
axs[0].set_ylabel('energia (s)')
axs[0].grid()


axs[1].plot(log_temperatura, log_energia, marker='s', linestyle='none', color='blue')
axs[1].set_xlabel('log(temperatura)')
axs[1].set_ylabel('log(energia)')
axs[1].grid()


#%% Regressão linear (semilog) Método dos mínimos quadrados



tempo =     np.array([0,      5,     10,    15,    20,    25,    30,     35,     40,     45])     # X
atividade = np.array([9.676 , 6.355, 4.261, 2.729, 1.862, 1.184, 0.7680, 0.4883, 0.3461, 0.2119]) # Y

# Apenas é utilizado o log_atividade !!! (semilog)
log_atividade = np.log(atividade)


fig, axs = plt.subplots(2, 1,figsize=(5,7), sharex=True) # sharex=True faria com que o x fosse partilhado

axs[0].plot(tempo, atividade, marker='s', linestyle='none', color='purple')
axs[0].set_ylabel('activity (mCi)')
axs[0].set_xlabel('time (s)')
axs[0].grid()


axs[1].plot(tempo, log_atividade, marker='s', linestyle='none', color='blue')
axs[1].set_ylabel('log(activity) (mCi)')
axs[1].set_xlabel('time (s)')
axs[1].grid()


m, b, r2, delta_m, delta_b = lin_reg(tempo, log_atividade, 1)
log_fit = m*tempo + b          # Fit for log function
axs[1].plot(tempo, log_fit,'g', linewidth=1)


normal_fit = np.exp(log_fit)   # Convert to normal function
axs[0].plot(tempo, normal_fit,'g', linewidth=1)


axs[0].set_title("Ex4 - Amostra do Isotopo Radioativo")


#%% Método de Euler (Sem RA) - 1 dimensão

t0 = 0
tf = 2
dt = 0.01

g = 9.8
y0 = 0
v0 = 10            # m/s
a0 = -g

Nt = int(np.ceil((tf - t0) / dt) + 1)

y = np.zeros((Nt,)) # np.zeros((Nt,2))  faria com 2 colunas
y[0] = y0
v = np.zeros((Nt,))
v[0] = v0
a = np.zeros((Nt,))

for i in range(Nt - 1):
    # F = -mg - D * |v| * v
    # Atençao que se for acelaração segund x é ax = -D|v|vx
    # Sem RA
    a[i + 1] = a0 # Porque é constante (acelaracao = g)
    v[i + 1] = v[i] + a[i]*dt 
    y[i + 1] = y[i] + v[i]*dt
    
#%% Método de Euler (Com RA) - 1 dimensão

t0 = 0
tf = 2
dt = 0.01

g = 9.8
y0 = 0
v0 = 10            # m/s
vt = 100 * 1000/3600  # 100 km/h --> m/s

D = g/(vt**2)

Nt = int(np.ceil((tf - t0) / dt) + 1)

t = np.linspace(t0, tf, Nt)
# t = np.arange(t0,tf+dt,dt)

y_RA = np.zeros((Nt,)) # np.zeros((Nt,2))  faria com 2 colunas
y_RA[0] = y0
v_RA = np.zeros((Nt,))
v_RA[0] = v0
a_RA = np.zeros((Nt,))


# Método de Euler
for i in range(Nt - 1):
    # F = -mg - D * |v| * v
    # Atençao que se for acelaração segund x é ax = -D|v|vx
    # Com RA
    a_RA[i] = -g-D*v_RA[i]**2                   # D = g/(vt**2)
    v_RA[i + 1] = v_RA[i] + a_RA[i]*dt 
    y_RA[i + 1] = y_RA[i] + v_RA[i]*dt


#%% Método de Euler (Complexo com RA) - 2 dimensão

D = g/(vt**2)

y = np.zeros((Nt,)) # np.zeros((Nt,2))  faria com 2 colunas
y[0] = y0
x = np.zeros((Nt,)) # np.zeros((Nt,2))  faria com 2 colunas
x[0] = 0  # x0

vx = np.zeros((Nt,))
vy = np.zeros((Nt,))

ax = np.zeros((Nt,))
ay = np.zeros((Nt,))


# Método de Euler
for i in range(Nt - 1):
    # F = -mg - D * |v| * v
    
    vv=np.sqrt(vx[i]**2 +vy[i]**2) # Intensidade
    
    ax[i]=-D*vv*vx[i]
    ay[i]=-g-D*vv*vy[i]
    
    vx[i+1]=vx[i]+ax[i]*dt
    vy[i+1]=vy[i]+ay[i]*dt
    
    x[i+1]=x[i]+vx[i]*dt
    y[i+1]=y[i]+vy[i]*dt

#%% Queda livre com RA (método de Euler)


t0 = 0
tf = 2
dt = 0.01

g = 9.8
y0 = 0
v0 = 10            # m/s
a0 = -g
vt = 100 * 1000/3600  # 100 km/h --> m/s

D = g/(vt**2)

Nt = int(np.ceil((tf - t0) / dt) + 1)

t = np.linspace(t0, tf, Nt)
# t = np.arange(t0,tf+dt,dt)


y = np.zeros((Nt,)) # np.zeros((Nt,2))  faria com 2 colunas
y[0] = y0
v = np.zeros((Nt,))
v[0] = v0
a = np.zeros((Nt,))

y_RA = np.zeros((Nt,)) # np.zeros((Nt,2))  faria com 2 colunas
y_RA[0] = y0
v_RA = np.zeros((Nt,))
v_RA[0] = v0
a_RA = np.zeros((Nt,))


# Método de Euler
for i in range(Nt - 1):
    # Sem RA
    a[i + 1] = a0 # Porque é constante (acelaracao = g)
    v[i + 1] = v[i] + a[i]*dt 
    y[i + 1] = y[i] + v[i]*dt
    
    # Com RA
    vv = np.sqrt(v_RA[i]**2 + 0) # Se quisesse fazer um vx e vy teria de trocar o 0 --> vx_RA
    a_RA[i] = -g-D*vv*v_RA[i]  # coloca-se o '-g' porque é em funcao de y
    v_RA[i + 1] = v_RA[i] + a_RA[i]*dt 
    y_RA[i + 1] = y_RA[i] + v_RA[i]*dt
    

# Max position

if(len(y == np.max(y)) >= 0):
    y_max = np.max(y)
    inters_max = np.where(y == np.max(y))[0][0]
    t_inters_max = t[inters_max]
    print("Máximo do objeto:",y_max,"em t =",t_inters_max)

if(len(np.where(y<=0)[0]) > 1):
    inters_ground = np.where(y<=0)[0][1]
    t_inters_ground = t[inters_ground]
    print("Posição inicial em t =",t_inters_ground)



# Plotting

plt.plot(t, y, '-r', linewidth=1) # Reta
plt.plot(t, y_RA, '-k', linewidth=1) # Reta




plt.xlabel('t (s)')
plt.ylabel('y (m)')
plt.grid()

plt.legend(['posicao sem RA','posicao com RA'])

plt.title("Ex1 - Queda livre")

#%% Método de Euler + Força de Magnus (Complexo com RA) - 2 dimensão

PI = np.pi
p_ar = 1.225  # kg/m3

t0 = 0
tf = 2
dt = 0.001

g = 9.8
x0 = -10
y0 = 1
v0 = 130 * 1000/3600  # 130 km/h --> m/s
ang = 10  # 10 graus
v0x = v0 * np.cos(np.radians(ang))
v0y = v0 * np.sin(np.radians(ang))


vt = 100 * 1000/3600  # 100 km/h --> m/s

m = 0.057  # kg
d = 0.067  # m   # Diametro
r = d / 2  # m   # Raio
A = PI * r**2


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


# Método de Euler
for i in range(Nt - 1):
    # 𝐹_𝑀𝑎𝑔𝑛𝑢s = 1/2 𝐴 𝜌_𝑎𝑟 𝑟 𝜔⃗ × 𝑣 
    # ang_w é o ang da rotacao
    
    ang_w = (0,0,100)
    F_magnus = (1/2) * A * p_ar * r * (np.cross(ang_w,(vx[i],vy[i],0)))
    ax_magnus = F_magnus[0]/m
    ay_magnus = F_magnus[1]/m
    
    vv=np.sqrt(vx[i]**2 +vy[i]**2) # Intensidade
    
    ax[i]=-D*vv*vx[i] + ax_magnus
    ay[i]=-g-D*vv*vy[i] + ay_magnus
    
    vx[i+1]=vx[i]+ax[i]*dt
    vy[i+1]=vy[i]+ay[i]*dt
    
    x[i+1]=x[i]+vx[i]*dt
    y[i+1]=y[i]+vy[i]*dt


# Max position

if(len(y == np.max(y)) >= 0):
    y_max = np.max(y)
    inters_max = np.where(y == np.max(y))[0][0]
    t_inters_max = t[inters_max]
    print("Altura máxima da bola: y = ",y_max,"em t =", )

if(len(np.where(y<=0)[0]) > 1):
    inters_ground = np.where(y<=0)[0][0]
    alcance = x[inters_ground]
    print("Alcance da bola: ",alcance)





#%% Cálculo simbólico Sympy

y,v,a,t,vt,g = sy.symbols('y,v,a,t,vt,g') # defenir nomes de variaveis

# Fazer as fórmulas
y = (((vt)**2)/g) * sy.log(sy.cosh(g*t/vt))

# Substituir os valores que conhecemos
y = y.subs([(g,9.8),(vt,6.8)])

v = sy.diff(y,t)
a = sy.diff(v,t)

# d)
a_exato = g - (g*v*np.abs(v))/(vt**2)
# Substituir os valores que conhecemos
a_exato = a_exato.subs([(g,9.8),(vt,6.8)])

# e) sem resistência do ar a = 9.8

y_semRA, v_semRA, a_semRA, t, v0_semRA, y0_semRA = sy.symbols('y_semRA,v_semRA,a_semRA,t,v0_semRA,y0_semRA') # definir nomes de variaveis

a_semRA = 9.8   # aceleração da gravidade

v_semRA = sy.integrate(a_semRA, (t, 0, t)) + v0_semRA

y_semRA = sy.integrate(v_semRA, (t, 0, t)) + y0_semRA


# Substituir os valores iniciais nas equações finais de v_semRA e y_semRA
v_semRA = v_semRA.subs(v0_semRA, 0)
y_semRA = y_semRA.subs([(y0_semRA, 0), (v0_semRA, 0)])

# Imprimir as equações de posição e velocidade
print("y_semRA = ", y_semRA)
print("v_semRA = ", v_semRA)

ts = np.linspace(0,4,100) # pontos para a funcao


# Criar uma função em de x
v_plot = sy.lambdify(t,v,"numpy")
a_plot = sy.lambdify(t,a,"numpy")
a_exato_plot = sy.lambdify(t,a_exato,"numpy")
y_semRA_plot = sy.lambdify(t,y_semRA,"numpy")
v_semRA_plot = sy.lambdify(t,v_semRA,"numpy")
a_semRA_plot = sy.lambdify(t,a_semRA,"numpy")


# Calcular o tempo que demora a atingir o solo

tfinal = sy.nsolve(y-20,t,0)
tfinal_semRA = sy.nsolve(y_semRA-20,t,0)
print("Interseção com ochão a 20 metros de altura")
print("Com Resistencia do ar: (",tfinal,", ",20,")",sep='')
print("Sem Resistencia do ar: (",tfinal_semRA,", ",20,")",sep='')

#%% Representar vetores centrados num ponto


# (x0,y0) = ponto inicial do vetor
# v1 = (x,y) = comprimentos x e y do vetor

PI = np.pi
F = 5
x0 = 0
y0 = 0


conv = PI/180 # converter graus em radianos G*conv = R
theta = np.array([PI/2, 60*conv, -7*PI/6, 310*conv])

vetores = np.array([F*np.cos(theta), F*np.sin(theta)])


fig, ax = plt.subplots()
ax.plot(x0,y0, 'o', markersize = 5)

for i in range(len(theta)):
    ax.arrow(x0, y0, vetores[0,i], vetores[1,i], color='r', width=0.05)
ax.set_aspect('equal')


# Confirmar com uma circunferencia de raio 5

# Creating equally spaced 100 data in range 0 to 2*pi
theta = np.linspace(0, 2 * np.pi, 100)

# Setting radius
radius = 5

# Generating x and y data
x = radius * np.cos(theta)
y = radius * np.sin(theta)

# Plotting
ax.plot(x, y)

plt.show()

#%% Vetores posicao e velocidade


# (x0,y0) = ponto inicial do vetor
# v1 = (x,y) = comprimentos x e y do vetor

x0,y0 = 0, 0

time = np.arange(1, 5, 1)
r = np.array([2*time, time**2])
v = np.array([2*np.ones(len(time)), 2*time])


fig, axs = plt.subplots(1, 2)
axs[0].plot(x0,y0, 'o', markersize = 5)

for i in range(len(time)):
    
    print(r[0,i])
    axs[0].arrow(x0, y0, r[0,i], r[1,i], color='k', width=0.1)
    axs[0].plot(r[0,i],r[1,i], 'ro', markersize = 5)
    axs[1].plot(r[0,i],r[1,i], 'ro', markersize = 5)
    axs[1].arrow(r[0,i], r[1,i], v[0,i], v[1,i], color='k', width=0.1)


axs[0].set_xlim(0,15)
axs[0].set_ylim(0,25)
axs[1].set_xlim(0,15)
axs[1].set_ylim(0,25)
    
axs[0].set_aspect('equal')



plt.show()


#%% Gravar em PNG 


plt.savefig(dir_figs+'esteroides_computacionais.png') # Save the graph in png

plt.show() 

# Abrir um separador com o grafico
# Tools > Preferences > IPython console > Graphics > Graphics Backend > Backend: Automatic

