# -*- coding: utf-8 -*-
"""
Created on Thu Mar  2 16:57:28 2023

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


#%% EX2: Volante de badmington

t0 = 0
tf = 4
dt = 0.01

g = 9.8

vt = 6.8

Nt = int(np.ceil((tf - t0) / dt) + 1)

t = np.linspace(t0, tf, Nt)
# t = np.arange(t0,tf+dt,dt)

# Volante de badmington

a_v = np.zeros((Nt,))
v_v = np.zeros((Nt,))



# Formulas para posição e acelaração exatas
y_exata_v = (((vt)**2)/g) * np.log(np.cosh(g*t/vt))





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


#%% Ploting

fig, axs = plt.subplots(3, 1,figsize=(5,7),sharex=True) # sharex=True faria com que o x fosse partilhado
t = np.linspace(t0, tf, Nt) # pontos para as funcao

axs[0].plot(t, y_exata_v, '-g', linewidth=1)
axs[0].plot(t,y_semRA_plot(t),'-r')
axs[0].plot(tfinal,20,'go')
axs[0].plot(tfinal_semRA,y_semRA_plot(tfinal_semRA),'ro')
axs[0].set_xlabel('t (s)')
axs[0].set_ylabel('y (m)')
axs[0].grid()
axs[0].legend(['Sem Resistência do ar','Com Resistência do ar'])


axs[1].plot(ts,v_plot(ts),'-g')
axs[1].plot(tfinal,v_plot(float(tfinal)),'go')
axs[1].plot(ts,v_semRA_plot(ts),'-r')
axs[1].plot(tfinal_semRA,v_semRA_plot(tfinal_semRA),'ro')
axs[1].set_xlabel('t (s)')
axs[1].set_ylabel('v(t)')
axs[1].grid()


axs[2].plot(ts,a_plot(ts),'-g')
axs[2].plot(tfinal,a_plot(float(tfinal)),'go')
a_semRA_plot_vec = np.vectorize(a_semRA_plot) # Como a_semRA_plot é uma constante precisamos de vetorizar
axs[2].plot(ts, a_semRA_plot_vec(ts), '-r')
axs[2].plot(tfinal_semRA,a_semRA_plot(tfinal_semRA),'ro')
axs[2].set_xlabel('t (s)')
axs[2].set_ylabel('a(t)')
axs[2].grid()

#axs[2].plot(ts,a_exato_plot(ts),'-r') # Função a_exato para verificar que é igual alínea d)


axs[0].set_title("Ex2 - Volante de badmington")

#%%


#------------  Gravar em PNG  ------------#


plt.savefig(dir_figs+'Ex2_volante_de_badmington.png') # Save the graph in png

plt.show() 

# Abrir um separador com o grafico
# Tools > Preferences > IPython console > Graphics > Graphics Backend > Backend: Automatic




