#!/Applications/anaconda3/envs/finelli/bin/python


# generare il pool di Su(3) ranom matrices (le M) -> farlo in un file separato e load here

#generare il lattice e per ogni punto mettere quattro matrici di SU3 inizializz a identità

#generatore del pool di matrici random di SU3
import numpy as np
from numpy.random import uniform 
import math
import cmath
from numpy import *
import matplotlib.pyplot as plt
#generarne 50 o 100
eps=0.24
Nmat=100 #quante matrici voglio generare

#generare matrice hermitiana H
def generate_hermitean():
    H=np.zeros((3,3),np.complex_)
    
   # for i in range(3):
    #    for j in range (i,3):
    #        if (i==j):
      #          diag=uniform(-1,1)
     #           H[i][i]= diag + 1j*diag
      #      else:
       #         nodiag=uniform(-1,1)
        #        H[i][j]=nodiag + 1j*nodiag
         #       H[j][i]=nodiag + 1j*nodiag
    for i in range(3):
        for j in range (0,3):
            H[i][j]= uniform(-1,1) + 1j*uniform(-1,1)
    H = (1/2)*(H+np.matrix.getH(H))
    return H

#taylor expansion e matrice unitaria come esponenziale di una hermitiana
def unitarization(H):
    M=np.zeros((3,3),np.complex_)
   
    for n in range(90):
        M+=(((1j*eps)**n)/(np.math.factorial(n))*np.linalg.matrix_power(H,n))
    #print(M)
    detM=np.linalg.det(M)**(1/3)
    M=M/detM
    #print(np.linalg.det(M))
    
    return M

def generate_pool(N_mat):
    MatList=np.zeros((N_mat*2,3,3),np.complex_)
    H=np.zeros((3,3),np.complex_)
    
    y = []
    x= range(Nmat*2)    
    
    for i in range(N_mat):
        H=generate_hermitean()
        MatList[i]=unitarization(H)
        #print(np.linalg.det(MatList[i]).real)
        MatList[Nmat+i]=np.matrix.getH(unitarization(H))
        #print(MatList[Nmat+i])
        #print(np.matmul(MatList[i],MatList[i+1]))
        #x[i] = i
        y.append(np.real(np.linalg.det(MatList[i])))
        #print(y[i])
        y.append(np.real(np.linalg.det(MatList[Nmat+i])))
        
        #debug per vedere se il det = 1 
    plt.plot(x,y,'ro')
    plt.axis([0,210,0,3])
    plt.show()
    print(np.average(y))
    print(np.std(y))
    return MatList


generate_pool(Nmat)




    #altro metodo
    #M=M+1j*eps*H
    #first=M[:,0]
    #second=M[:,1]
    #third=M[:,2]
    
    #fnorm=np.linalg.norm(first)
    #first=first/fnorm
    #print(np.linalg.norm(first)) #debug ok è 1
    #for i in range(3):
        #second[i]=second[i]-(np.dot(first,np.conjugate(second))*first[i])
    #snorm=np.linalg.norm(second)
    #second=second/snorm
    #print(np.linalg.norm(second)) #debug ok è 1
    #print(np.dot(np.conjugate(first),second)) #è ok????
    #exp=np.zeros((3,3),np.complex_)
     #third=np.cross(np.conjugate(first),second)
    #print(np.dot(third,second))
    #M[:,0]=first
    #M[:,1]=second
    #M[:,2]=third
    #q,r = np.linalg.qr(M)
    #print(np.linalg.det(M))