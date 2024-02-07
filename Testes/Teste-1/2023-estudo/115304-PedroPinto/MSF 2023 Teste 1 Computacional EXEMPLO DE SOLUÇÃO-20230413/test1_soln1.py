#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 20 10:32:27 2023

@author: gareth
"""
import numpy as np
import matplotlib.pyplot as plt


def linreg_MQ(x,y):
# regressão linear usando o método dos minimos quadrados
# retorna melhor valores de m e b para uma reta y = m*x+b
# com as incertezas dm e db, e o coeficiente de determinação r^2
    x = np.array(x) #convertir inputs para array numpy
    y = np.array(y)
    
    N = len(x)
    if N < 3:
        raise ValueError("ERROR linreg_MQ: N deve ser maior do que 2")

    else:
        sumxy = sum(x*y)
        sumx = sum(x)
        sumy = sum(y)
        sumx2 = sum(x**2.0)
        sumy2 = sum(y**2.0)
        
        # para verificar os cálculos intermédios
        #print(N,sumxy,sumx,sumy,sumx2,sumy2)
        denom = N*sumx2 - sumx**2.0
        m = (N*sumxy - sumx*sumy)/denom
        b= (sumx2*sumy - sumx*sumxy)/denom
        r2 = (N*sumxy-sumx*sumy)**2.0/(denom*(N*sumy2-sumy**2.0))
        dm = abs(m)*np.sqrt((1.0/r2-1.0)/(N-2))
        db = dm*np.sqrt(sumx2/N)
        return [m,b,dm,db,r2]


#dados
t = np.linspace(0,8*48,9)
act =np.array([10.03,  7.06,  4.88,  3.38,  2.26,  1.66,   1.14,  0.79,  0.58])

#--1a--
plt.figure(1)
plt.plot(t,act,'.') # plto dos dados originais
m,b,dm,db,r2 = linreg_MQ(t,act) #regressão linear
plt.plot(t,m*t+b) # adicionar fit ao plot

plt.xlabel('t (h)')
plt.ylabel('atividade (mBq)')
plt.savefig('fig_1a.png')
print("linear:")
print(f"coeficiente de determinação: r^2 = {r2:.4f}")
#resposta:
# a relação entre a atividade e o tempo  não ´w linear, pois visualmente vemos
# que os dados não seguem muito bem uma reta, que é confirmada pelo valor da
# coeficiente de determinação, que não é perto a 1.

#--1b--
plt.figure(2)
logact = np.log(act) #logaritmo natural da atividade
plt.plot(t,logact,'.') # plot dos dados originais
m,b,dm,db,r2 = linreg_MQ(t,logact) #regressão linear
plt.plot(t,m*t+b) # adicionar fit ao plot

plt.xlabel('t (h)')
plt.ylabel('log atividade (log mBq)')
plt.savefig('fig_1b.png')
print("log-linear:")
print(f"declive: m = {m:.5f} +/- {dm:.5f} log(mBq)/h")
print(f"coeficiente de determinação: r^2 = {r2:.4f}")

#--1c--
thalf = -np.log(2)/m
dthalf = thalf*dm/abs(m) #erro relativo em t1/2 = erro relativo em m

print(f"tempo meia-vida: t_1/2 = {thalf:.1f} +/- {dthalf:.1f} horas")





