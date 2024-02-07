import numpy as np
import matplotlib.pyplot as plt


# parámetros e condições iniciais
tf=1.5
t0 = 0.0

x0 = 0
y0 = 2.0

v0 = 15
thet = 30*np.pi/180
v0x = v0*np.cos(thet)
v0y = v0*np.sin(thet)

g=9.80
vt = 20
dres=g/vt**2

for dt in [0.1,0.01,0.001,0.0001]:			# experimentar várias passos de tempo

    #criar variáveis
    n=int((tf-t0)/dt+0.1)	# +0.1 para garantir não arredondar para baixo
    t=np.zeros(n+1)		# n+1 elementos; último índice n
    x=np.zeros(n+1)
    vx=np.zeros(n+1)
    y=np.zeros(n+1)
    vy=np.zeros(n+1)
    ax=np.zeros(n+1)
    ay=np.zeros(n+1)
    
    #valores iniciais
    vx[0]=v0x
    vy[0]=v0y
    t[0]=t0
    x[0]=x0
    y[0]=y0
    
    #método de Euler
    for i in range(n):
        t[i+1]=t[i]+dt
        vv=np.sqrt(vx[i]**2+vy[i]**2)
        ax[i]=-dres*vv*vx[i]
        ay[i]=-g-dres*vv*vy[i]
        vx[i+1]=vx[i]+ax[i]*dt  
        vy[i+1]=vy[i]+ay[i]*dt             
        x[i+1]=x[i]+vx[i]*dt  
        y[i+1]=y[i]+vy[i]*dt 
    
    

    #valores só em cima das 3m
    yup = y[y>3.0]
    xup = x[y>3.0]
    tup = t[y>3.0]

    #vemos qual o dt necessário para convergir
    print("dt",dt, " x",xup[-1]," t",tup[-1])


plt.plot(x,y)
plt.xlabel('x (m)')
plt.ylabel('y (m)')
    
plt.plot(xup[-1],yup[-1],'o') #markar o ponto onde passa a altura do cesto

plt.savefig('fig_2.png')

print(f"a bola passa a altura do cesto (3m) a descer a {xup[-1]:.2f}m em {tup[-1]:.2f}s")

#print(max(y))