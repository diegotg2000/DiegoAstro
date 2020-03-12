import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import entropy

def entropia(A):
    return entropy(A,base=2)

def D(A,B):
    return entropia((A+B)/2)-entropia(A)/2-entropia(B)/2         

def H(A):
    return entropia(A)/np.log2(len(A))

def Dstar(n):
    return -0.5*((n+1)/n*np.log2(n+1)+np.log2(n)-2*np.log2(2*n))

def complejidad(A,B):
    n=len(A)
    return D(A,B)*H(A)/Dstar(n)


BETAS=np.arange(10,51,5)/10

plt.figure(figsize=(5,4))
for j in range(10,110,10):
    Com=np.zeros(len(BETAS))
    for i in range(len(BETAS)):
        #se lee los datos random
        A=np.loadtxt('random'+str(j)+'_'+str(BETAS[i])+'.BSKIndex', usecols=0)
        B=np.loadtxt('random'+str(j)+'_'+str(BETAS[i])+'.BSKIndex', usecols=1)
        C=np.concatenate((A,B),axis=0)
    
        nodos_ran, conexiones_ran = np.unique(C, return_counts=True)
    
        #se lee los datos simulados
        A=np.loadtxt('norandom'+str(j)+'_'+str(BETAS[i])+'.BSKIndex', usecols=0)
        B=np.loadtxt('norandom'+str(j)+'_'+str(BETAS[i])+'.BSKIndex', usecols=1)
        C=np.concatenate((A,B),axis=0)
    
        nodos_dat, conexiones_dat = np.unique(C, return_counts=True)
    
        #se hace el histograma del numero de conexiones
    
        max_con=max([max(conexiones_dat),max(conexiones_ran)])
    
        y_ran,x_ran=np.histogram(conexiones_ran, bins=max_con, range=(1,max_con+1))
    
        y_dat,x_dat=np.histogram(conexiones_dat, bins=max_con, range=(1,max_con+1))
    
        #se analiza random
        puntos_ran=np.loadtxt('random'+str(j)+'.dat',usecols=0)
        total_ran=len(puntos_ran)
        nonzero_ran=len(nodos_ran)
        p0_ran=np.array([(total_ran-nonzero_ran)])
        P=np.concatenate((p0_ran,y_ran),axis=0)
        y_ran=P/total_ran
        #se analiza datos
        puntos_dat=np.loadtxt('norandom'+str(j)+'.dat',usecols=0)
        total_dat=len(puntos_dat)
        nonzero_dat=len(nodos_dat)
        p0_dat=np.array([(total_dat-nonzero_dat)])
        P=np.concatenate((p0_dat,y_dat),axis=0)
        y_dat=P/total_dat

    
        if len(x_dat)==len(x_ran):   
            Com[i]=complejidad(y_dat,y_ran)
        else:
            print('ojo!')
            
    np.savetxt('p_{}'.format(j), Com)


    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    