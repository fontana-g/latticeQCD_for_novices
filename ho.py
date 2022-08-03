#!/Applications/anaconda3/envs/finelli/bin/python

import numpy as np
from numpy.random import uniform 
import matplotlib.pyplot as plt
from math import *

N=20
N_cor = 20
N_cf = 1000
a = 0.5
epsilon = 1.4
m=1
x = np.zeros((N,),np.float64)
G = np.zeros((N_cf,N),np.float64)


#action with only the terms concerning a certain x_j (harmonic oscillator)
def S(j,x,m,a):
    jp = (j+1)%len(x)    # next site
    jm = (j-1)%len(x)    # previous site
    return a*((x[j])**2)/2 + x[j]*(x[j]-x[jp]-x[jm])*(m/a)

#Metropolis sweep for x (x is an array of length N)
#meanwhile it calculates the acceptance rate for this x
def update(x,m,a,epsilon):
    accept=0
    
    for i in range(len(x)):
        x_old=x[i]
        S_old=S(i,x,m,a)
        
        noise=uniform(-epsilon,epsilon)
        x[i]+=noise
        S_new=S(i,x,m,a)
        dS= S_new-S_old
        if dS > 0 and exp(-dS)<uniform(0,1):
            x[i]=x_old
        else:
            accept += 1
            
    return x, accept/len(x)

def computeG(x,n): #compute G_n
    g=0
    N=len(x)
    for i in range(N):
        g+=x[i]*x[(i+n)%N]
    return g/N

def MCaverage(x,G,N_corr,N_cf):
    N=len(x)
    for j in range(10*N_corr):
        update(x,m,a,epsilon)
    for alpha in range(N_cf):
        for i in range(N_corr):
            update(x,m,a,epsilon)
        for n in range(N):
            G[alpha][n] = computeG(x,n)
    #in the end we have a matrix with N_cf rows and len(x) columns
    #we have to compute averages summing over alphas
    G_avgs=[]
    for n in range(N):
        avg=0
        for alpha in range(N_cf):
            avg+=G[alpha][n]
        G_avgs.append(avg/N_cf)
    print(f"averages: {G_avgs}")
    return G_avgs
#returns an array of averages of length N [Gavg0, ..., Gavg(N-1)]
def DeltaE(G_avgs):
    dE=np.zeros((N),np.float64)
    dE = np.log(np.divide(G_avgs[:-1],G_avgs[1:]))
    return dE/a

DE=np.zeros((N),np.float64)
DE=DeltaE(MCaverage(x,G,N_cor,N_cf))

xaxis=range(N-1)
yaxis=DE
plt.plot(xaxis,yaxis)

plt.show()