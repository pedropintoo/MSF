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
    
    
#%% Método de Euler (Complexo com RA) - 2 dimensão Resistencia do ar


t0 = 0
tf = 1.5
dt = 0.1

g = 9.8
x0 = 0
y0 = 0
v0 = 100 * 1000 / 3600  # 130 km/h --> m/s
ang = 10  # 10 graus
v0x = v0 * np.cos(np.radians(ang))
v0y = v0 * np.sin(np.radians(ang))


vt = 100 * 1000 / 3600 # m/s



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
   
    
    vv=np.sqrt(vx[i]**2 +vy[i]**2) # Intensidade
    
    ax[i]=-D*vv*vx[i]
    ay[i]=-g-D*vv*vy[i]
    
    vx[i+1]=vx[i]+ax[i]*dt
    vy[i+1]=vy[i]+ay[i]*dt
    
    x[i+1]=x[i]+vx[i]*dt
    y[i+1]=y[i]+vy[i]*dt


inters_ground = np.where(y == y.max())[0][0]
alcance = x[inters_ground]
print("Altura máxima da bola: ",alcance)

tempoal = t[inters_ground]
print("No instante: ",tempoal)

print("\n")


inters_ground = np.where(y <= 0)[0][1]
alcance = x[inters_ground]
print("Alcance da bola: ",alcance)

tempoal = t[inters_ground]
print("No instante: ",tempoal)



# Plotting



plt.plot(x, y, '-r', linewidth=1) # Reta
plt.xlabel('x (m)')
plt.ylabel('y (m)')
plt.grid()




#%% Método de Euler + Força de Magnus + Resistencia do Ar (Complexo com RA) - 3 dimensão


PI = np.pi
p_ar = 1.225  # kg/m3

t0 = 0
tf = 10
dt = 0.0001

g = 9.8

x0 = 0
y0 = 0
z0 = 0

ang = 16

v0 = 100 * 1000 / 3600 # 100 km/h --> m/s
v0x = v0 * np.cos(np.deg2rad(ang))
v0y = v0 * np.sin(np.deg2rad(ang))
v0z = 0


vt = 100 * 1000/3600  # 100 km/h --> m/s

m = 0.45  # kg
r = 0.11  # m   # Raio
A = PI * r**2


D = g/(vt**2)

Nt = int(np.ceil((tf - t0) / dt) + 1)

t = np.linspace(t0, tf, Nt)

y = np.zeros((Nt,)) # np.zeros((Nt,2))  faria com 2 colunas
y[0] = y0

x = np.zeros((Nt,)) # np.zeros((Nt,2))  faria com 2 colunas
x[0] = x0  

z = np.zeros((Nt,)) # np.zeros((Nt,2))  faria com 2 colunas
z[0] = z0  

vx = np.zeros((Nt,))
vx[0] = v0x
vy = np.zeros((Nt,))
vy[0] = v0y
vz = np.zeros((Nt,))
vz[0] = v0z

ax = np.zeros((Nt,))
ay = np.zeros((Nt,))
az = np.zeros((Nt,))


# Método de Euler
for i in range(Nt - 1):
    # 𝐹_𝑀𝑎𝑔𝑛𝑢s = 1/2 𝐴 𝜌_𝑎𝑟 𝑟 𝜔⃗ × 𝑣 
    # ang_w é o ang da rotacao
    
    
    # F_magnus = (1/2) * A * p_ar * r * (np.cross(ang_w,(vx[i],vy[i],vz[i])))
    
    # Neste caso !!!
    ang_w = (0,0,-10)
    
    F_magnus = (1/2) * A * p_ar * r * (np.cross(ang_w,(vx[i],vy[i],vz[i])))
    
    ax_magnus = F_magnus[0]/m
    ay_magnus = F_magnus[1]/m  # 0
    az_magnus = F_magnus[2]/m
    
    vv=np.sqrt(vx[i]**2 + vy[i]**2 + vz[i]**2) # Intensidade
    
    ax[i]=-D*vv*vx[i] + ax_magnus
    ay[i]=-g-D*vv*vy[i] + ay_magnus
    az[i]=-D*vv*vz[i] + az_magnus

    
    vx[i+1]=vx[i]+ax[i]*dt
    vy[i+1]=vy[i]+ay[i]*dt
    vz[i+1]=vz[i]+az[i]*dt
    
    x[i+1]=x[i]+vx[i]*dt
    y[i+1]=y[i]+vy[i]*dt
    z[i+1]=z[i]+vz[i]*dt
    
    if x[i] > 20 :
        print("x = ",x[i])
        print("y = ",y[i])
        print("z = ",z[i])
        if (-3.75 < z[i] < 3.75 and 0 < y[i] < 2.4): 
            print("A bola entrou!")
        else:
            print("A bola nao entrou!")
        break;
        



# Plotting

# fig, axs = plt.subplots(3, 1,figsize=(5,7)) # sharex=True faria com que o x fosse partilhado


# axs[0].plot(t, x, '-r', linewidth=1) # Reta
# axs[0].set_xlabel('t (s)')
# axs[0].set_ylabel('x (m)')
# axs[0].grid()

# axs[1].plot(t, y, '-r', linewidth=1) # Reta
# axs[1].set_xlabel('t (s)')
# axs[1].set_ylabel('y (m)')
# axs[1].grid()

# axs[2].plot(t, z, '-r', linewidth=1) # Reta
# axs[2].set_xlabel('t (s)')
# axs[2].set_ylabel('z (m)')
# axs[2].grid()


plt.figure(figsize=(8,8))
ax = plt.axes(projection='3d')
ax.plot3D(x[x>=0],-z[x>=0],y[x>=0], 'r')
goalx = [0,0,0,0]
goaly = [0,2.4,2.4,0]
goalz = [-3.66,-3.66,3.66,3.66]
ax.plot3D(goalx,goalz,goaly, 'k')
ax.set_xlim3d(0, 15)
ax.set_ylim3d(-25, 25)
ax.set_zlim3d(0, 5)
ax.set_box_aspect((2,6,2))
ax.set_xlabel('x')
ax.set_ylabel('z')
ax.set_zlabel('y')




#%% Método de Euler Cromer (Complexo com RA) - Oscilador Harmónico Simples


t0 = 0
tf = 60
dt = 0.001

g = 9.8
x0 = 4
v0 = 0  

A = 4 # Amplitude máxima

m = 1  # kg
k = 1   # Constante elástica
w = np.sqrt(k/m)  # velocidade angular

Nt = int(np.ceil((tf - t0) / dt) + 1)

t = np.linspace(t0, tf, Nt)


x = np.zeros((Nt,)) # np.zeros((Nt,2))  faria com 2 colunas
x[0] = x0  

v = np.zeros((Nt,))
v[0] = v0

a = np.zeros((Nt,))

v_exata = - A * w * np.sin(w * t)

# Método de Euler Cromer
for i in range(Nt - 1):


    a[i]=-k/m * x[i] # Fx = -k * x(t) --> a = -k/m * x(t)

    
    v[i+1]=v[i]+a[i]*dt
    
    
    x[i+1]=x[i]+v[i+1]*dt
    


# Plotting

plt.plot(t, v, '-r', linewidth=1) # Reta
plt.plot(t, v_exata, '-b', linewidth=1) # Reta



plt.xlabel('t (s)')
plt.ylabel('v (m/s^2)')
plt.grid()





#%% Método de Euler Cromer (Complexo com RA) - Evolução temporal de ciclista com empurao

par = 1.225  # densidade do ar

t0 = 0
tf = 200
dt = 0.001

g = 9.8
x0 = 0
v0 = 0.5  # empurao de 1m/s
a0 = 0  

pot = 0.48 *  735.4975 # potencia (W)

m = 72  # kg
Catr = 0.01    # coeficiente de resistência 𝜇 de um piso liso de alcatrão
Cres = 0.9  # coeficiente de resistência do ar
A = 0.50 # area frontal

D = (Cres/2) * A * par

Nt = int(np.ceil((tf - t0) / dt) + 1)

t = np.linspace(t0, tf, Nt)


x = np.zeros((Nt,)) # np.zeros((Nt,2))  faria com 2 colunas
x[0] = x0  

v = np.zeros((Nt,))
v[0] = v0

a = np.zeros((Nt,))
a[0] = a0


# Método de Euler Cromer
for i in range(Nt - 1):
    
    v_norma = np.sqrt(v[i]**2 - 0) # Norma da velocidade

    # D = D = (Cres/2) * A * par
    # Px = m * g * sin(ang)  # Peso segundo Ox (com inclinação de ang)
    # N = m * g * cos(ang)   # Forca normal (com inclinação de ang)
    
    # F = Fcic + Px + Fres + Frol + Nx (segundo x) 
    # a = 1/m * (Potencia/velocidade - N * Catr - Px - D * velocidade * v_norma)

    a[i] = 1 / m * (pot/v[i] - m * g * np.cos( np.deg2rad(0)) * Catr - D * v[i] * v_norma )  # Px = 0 neste caso e N = m*g porque ang = 0
 
    
    v[i+1]=v[i]+a[i]*dt
    
    
    x[i+1]=x[i]+v[i+1]*dt
   
    
v_terminal = v[-1]
print("a) Velocidade terminal: ",v[-1]," (m/s)")



tempo_v = np.where(v>=v_terminal*0.9)[0][0]
t_tempo_v = t[tempo_v]
print("b) Depois de ",t_tempo_v," segundos")

tempo_v = np.where(x>=2000)[0][0]
t_tempo_v = t[tempo_v]
print("c) Percorre 2000 metros depois de",t_tempo_v," segundos")


# P = F * V

V = 30 * 1000 / 3600
P = 298
F = P / V
print("a) ",F)

V = 40 * 1000 / 3600
P = 298
F = P / V
print("b) ",F)

   
    

# Plotting

plt.plot(t, v, '-r', linewidth=1) # Reta



plt.xlabel('t (s)')
plt.ylabel('v (m/s)')
plt.grid()



#%% Método de Euler Cromer (Complexo com RA) - (considerando inclinação) Evolução temporal de ciclista com empurao

par = 1.225  # densidade do ar

t0 = 0
tf = 400
dt = 0.001

g = 9.8
x0 = 0
v0 = 0.5  # empurao de 1m/s
a0 = 0  

ang = 4 # graus

pot = 0.48 *  735.4975 # potencia (W)

m = 72  # kg
Catr = 0.01    # coeficiente de resistência 𝜇 de um piso liso de alcatrão
Cres = 0.9  # coeficiente de resistência do ar
A = 0.50 # area frontal

D = (Cres/2) * A * par



Nt = int(np.ceil((tf - t0) / dt) + 1)

t = np.linspace(t0, tf, Nt)


x = np.zeros((Nt,)) # np.zeros((Nt,2))  faria com 2 colunas
x[0] = x0  

v = np.zeros((Nt,))
v[0] = v0

a = np.zeros((Nt,))
a[0] = a0


# Método de Euler Cromer
for i in range(Nt - 1):
    
    v_norma = np.sqrt(v[i]**2 - 0) # Norma da velocidade

    # D = D = (Cres/2) * A * par
    # Px = m * g * sin(ang)  # Peso segundo Ox (com inclinação de ang)
    # N = m * g * cos(ang)   # Forca normal (com inclinação de ang)
    
    # Fx = Femp - Nx
    
    # a = 1/m * (Potencia/velocidade - N * Catr - Px - D * velocidade * v_norma)
    
    Px = m * g * np.sin(np.deg2rad(ang)) # Peso (com inclinação de 5 graus)
    N = m * g * np.cos( np.deg2rad(ang))
    
    a[i] = 1 / m * (pot/v[i] - N * Catr - Px - D * v[i] * v_norma ) 
 
    
    v[i+1]=v[i]+a[i]*dt
    
    
    x[i+1]=x[i]+v[i+1]*dt
    
    if(x[i] > 1500):
        ang = -1    # Desce com inclinacao de 1 grau
        # a = 1/m * (Potencia/velocidade - N * Catr + Px - D * velocidade * v_norma)
    
    
    
# Plotting

plt.plot(t, v, '-r', linewidth=1) # Reta



plt.xlabel('t (s)')
plt.ylabel('v (m/s)')
plt.grid()

tempo_v = np.where(x>=2000)[0][0]
t_tempo_v = t[tempo_v]
print("b) Percorre 2000 metros depois de",t_tempo_v," segundos")


v_terminal = v[-1]
print("c) Velocidade terminal: ",v[-1]," (m/s)")


#%% Método de Euler-Cromer (Complexo com RA) - 2 dimensão - Trabalho da resistencia do ar + Energia mecanica

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

   
inters_ground = np.where(y <= 0)[0][1]
alcance = x[inters_ground]
print("Alcance da bola: ",alcance," em t =",t[inters_ground],"s") 
    
iInicial = 0
iFinal = t[inters_ground]

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


plt.savefig(dir_figs+'esteroides_computacionais.png') # Save the graph in png

plt.show() 

# Abrir um separador com o grafico
# Tools > Preferences > IPython console > Graphics > Graphics Backend > Backend: Automatic



