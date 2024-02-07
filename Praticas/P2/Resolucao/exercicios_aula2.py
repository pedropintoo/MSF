# -*- coding: utf-8 -*-
"""

MSF - Aula prática 2
Nuno Monteiro

"""

#%% Importar libraries e diretório para figuras

import matplotlib.pyplot as plt
import numpy as np
import os

dir_figs = './figures/'
if not os.path.exists(dir_figs): 
    os.mkdir(dir_figs)

#%% EX1: Difração fenda única

# L and X datasets (in cm)
L = np.array([222.0, 207.5, 194.0, 171.5, 153.0, 133.0, 113.0, 92.0])
X = np.array([2.3, 2.2, 2.0, 1.8, 1.6, 1.4, 1.2, 1.0])


# "quick" method for plotting
# plt.figure()
# plt.plot(L,X,'r.')
# plt.ylabel('X (cm)')
# plt.xlabel('L (cm)')
# plt.grid()
# plt.ylim(0.8,2.5)
# plt.xlim(80,240)
# plt.savefig(dir_figs+'fenda_unica_obs.png')
# plt.show()

# plotting
fig, ax = plt.subplots() # create figure and plot axis
ax.plot(L, X, 'r.')
ax.set_xlabel('L (cm)')
ax.set_ylabel('X (cm)')
ax.set_title('L vs. X')
ax.grid()
fig.savefig(dir_figs+'fenda_unica_obs.png')
plt.show()


# linear regression function
def lin_reg(x, y, print_res):
    
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
    dm = np.abs(m) * np.sqrt((1 / r2 - 1) / (N - 2))
    db = dm * np.sqrt(s_x2 / N)
    
    if print_res == 1:
        print("s_x, s_y, s_xy, s_x2, s_y2, m, b, r2, dm, db:")
        print([s_x, s_y, s_xy, s_x2, s_y2, m, b, r2, dm, db])
    
    return m, b, r2, dm, db

# get and print results
m, b, r2, dm, db = lin_reg(x=L, y=X, print_res=1)

X_fit = L * m + b
fig, ax = plt.subplots()
ax.plot(L, X, 'r.')
ax.plot(L,X_fit,"b")
ax.set_xlabel('L (cm)')
ax.set_ylabel('X (cm)')
ax.set_title('L vs. X')
ax.grid()
fig.savefig(dir_figs+'fenda_unica_fit.png')
plt.show()

# e) Encontre o valor de X, quando L = 165.0 cm. Use a reta determinada pela regressão linear.

print(f'Valor de X(L = 165.0 cm) = {165 * m + b:.1f} cm')

# f) Afaste da reta encontrada um dos valores medidos de y (X).
# Compare o coeficiente de determinação com o valor anterior. Faça
# um gráfico com os novos pontos experimentais e a nova reta.

L = np.array([222.0, 207.5, 194.0, 171.5, 153.0, 133.0, 113.0, 92.0])
X = np.array([2.3, 2.2, 5.0, 1.8, 1.6, 1.4, 1.2, 1.0])

m, b, r2, dm, db = lin_reg(x=L, y=X, print_res=0)
X_fit = L * m + b

fig, ax = plt.subplots()
ax.plot(L, X, 'r.')
ax.plot(L,X_fit,"b")
ax.set_xlabel('L (cm)')
ax.set_ylabel('X (cm)')
ax.set_title('L vs. X')
ax.grid()
fig.savefig(dir_figs+'fenda_unica_fit_2.png')
plt.show()

print(f"new R^2 = {r2:.2f}")

#%% EX2: ciclista

# data
t = np.linspace(0,9,10)
# t = np.arange(0,10)
d = np.array([0.00, 0.735, 1.363, 1.739, 2.805, 3.814, 4.458, 4.955, 5.666, 6.329])

# plotting
fig, axs = plt.subplots(2, 1, sharex=True,figsize=(5,7))
ax = axs[0]
ax.plot(t, d, marker='s', linestyle='none', color='purple')
ax.set_ylabel('d (km)')
ax.grid()

m, b, r2, dm, db = lin_reg(t, d, 0)

# using polyfit to get lin reg parameters
malt,balt = np.polyfit(t, d, 1)
print(m, malt, b, balt)

dfit = m * t + b

# plotting
ax = axs[1]
ax.plot(t, d, marker='s', linestyle='none', color='purple')
ax.plot(t, dfit, linestyle='-', color='green')

ax.set_xlabel('t (s)')
ax.set_ylabel('d (km)')
ax.grid()
# place a text box in upper left in axes coords
textstr = f'$d = {m:.3f} \cdot t + ({b:.3f})$ \n'+f'$R^2={r2:.3f}$'
props = dict(boxstyle='round', facecolor='lightgreen', alpha=.5)
ax.text(0.05, 0.95, textstr, transform=ax.transAxes,
        verticalalignment='top', bbox=props)

fig.suptitle('d vs. t',fontsize=20, y=.93)
plt.subplots_adjust(hspace=0.05)
fig.savefig(dir_figs+'/ciclista.png')
plt.show()

# É uma relação linear bem aproximada? 
# O ciclista conseguiu manter a mesma velocidade uniforme 
# durante o percurso?

# É uma relação linear bem aproximada?
print(f'R^2 = {r2:.3f}')

# Qual a velocidade média do ciclista?
print(f'v = {m:.3f} +/- {dm:.3f} m/s')
print(f'erro relativo: {(dm/m)*100:.2f}%')

# polyfit

mp, bp = np.polyfit(t, d, 1)
dfit = mp*t + bp
fig, ax = plt.subplots(figsize=(6,6))
ax.plot(t, d, marker='s', linestyle='none', color='purple')
ax.plot(t, dfit, linestyle='-', color='green')

ax.set_xlabel('t (s)')
ax.set_ylabel('d (km)')
ax.set_title('d vs. t (polyfit)')
ax.grid()
# place a text box in upper left in axes coords
textstr = f'$d = {mp:.3f} \cdot t + ({bp:.3f})$ \n'+f'$R^2={r2:.3f}$'
props = dict(boxstyle='round', facecolor='lightgreen', alpha=.5)
ax.text(0.05, 0.95, textstr, transform=ax.transAxes, fontsize=12,
        verticalalignment='top', bbox=props)
# ax.set(aspect=1)
fig.savefig(dir_figs+'/ciclista_fit_polyfit.png')
plt.show()

# Apresente a velocidade em km/hora
print(f"v = {m*60:.2f} km/h")

#%% EX3: potência corpo negro

# a) Apresente estas medições num gráfico. A analisar o gráfico, a relação entre a energia emitida e a temperatura é linear?
# b) Apresente as medições num gráfico log-log. Qual a dependência entre as quantidade energia emitida e a temperatura?

T = np.array([200, 300, 400, 500, 600, 700, 800, 900, 1000, 1100])
E = np.array([0.6950, 4.363, 15.53, 38.74, 75.08, 125.2, 257.9, 344.1, 557.4, 690.7])

logT = np.log(T)
logE = np.log(E)
m, b, r2, dm, db = lin_reg(logT, logE, print_res=1)
logEfit = m * logT + b

fig, axs = plt.subplots(2,1)
axs[0].plot(T,E,'ko')
axs[1].plot(logT,logE,'k^')
axs[1].plot(logT,logEfit,'b-')
plt.savefig(dir_figs+'emissao_corpo_negro_1.png')
plt.show()

Efit = np.exp(logEfit)
fig, ax = plt.subplots()
ax.plot(T,E,'k.-')
ax.plot(T,Efit,'b^-')
plt.savefig(dir_figs+'emissao_corpo_negro_2.png')
plt.show()

