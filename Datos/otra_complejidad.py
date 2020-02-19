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

NOMBRES=[]
BETAS=np.arange(10,51)/10
for i in range(len(BETAS)):
    NOMBRES.append(['random'+str(BETAS[i])+'.BSKIndex','norandom'+str(BETAS[i])+'.BSKIndex']) 
    
Srandom=np.zeros(len(BETAS))
Snorandom=np.zeros(len(BETAS))
Crandom=np.zeros(len(BETAS))
Cnorandom=np.zeros(len(BETAS))

for i in range(len(NOMBRES)):
    #se lee los datos random
    A=np.loadtxt(NOMBRES[i][0], usecols=0)
    B=np.loadtxt(NOMBRES[i][0], usecols=1)
    C=np.concatenate((A,B),axis=0)
    
    nodos_ran, conexiones_ran = np.unique(C, return_counts=True)
    
    #se lee los datos simulados
    A=np.loadtxt(NOMBRES[i][1], usecols=0)
    B=np.loadtxt(NOMBRES[i][1], usecols=1)
    C=np.concatenate((A,B),axis=0)
    
    nodos_dat, conexiones_dat = np.unique(C, return_counts=True)
    
    #se hace el histograma del numero de conexiones
    
    max_con=max([max(conexiones_dat),max(conexiones_ran)])
    
    y_ran,x_ran=np.histogram(conexiones_ran, bins=max_con, range=(1,max_con+1))
    
    y_dat,x_dat=np.histogram(conexiones_dat, bins=max_con, range=(1,max_con+1))
    
    #se analiza random
    puntos_ran=np.loadtxt('random.dat',usecols=0)
    total_ran=len(puntos_ran)
    nonzero_ran=len(nodos_ran)
    p0_ran=np.array([(total_ran-nonzero_ran)])
    P=np.concatenate((p0_ran,y_ran),axis=0)
    y_ran=P/total_ran
    #se analiza datos
    puntos_dat=np.loadtxt('norandom.dat',usecols=0)
    total_dat=len(puntos_dat)
    nonzero_dat=len(nodos_dat)
    p0_dat=np.array([(total_dat-nonzero_dat)])
    P=np.concatenate((p0_dat,y_dat),axis=0)
    y_dat=P/total_dat
    #calculo la entropia
    Srandom[i]=entropia(y_ran)
    Snorandom[i]=entropia(y_dat)
    
    if len(x_dat)==len(x_ran):   
        Crandom[i]=complejidad(y_ran,y_dat)
        Cnorandom[i]=complejidad(y_dat, y_ran)
    else:
        print('ojo!')
                
plt.figure(figsize=(5,4))
plt.plot(BETAS, Srandom, label='Random',marker='o')
plt.plot(BETAS, Snorandom, label='NoRandom',marker='o')
plt.xlabel('\u03b2')
plt.ylabel('S')
plt.legend()
plt.savefig('SvBeta_otra_def2.png')

plt.figure(figsize=(5,4))
plt.plot(BETAS, Crandom, label='Random',marker='o')
plt.plot(BETAS, Cnorandom, label='NoRandom',marker='o')
plt.xlabel('\u03b2')
plt.ylabel('C')
plt.legend()
plt.savefig('CvBeta_otra_def2.png')

plt.figure(figsize=(5,4))
plt.plot(Srandom, Crandom, label='Random',marker='o')
plt.plot(Snorandom, Cnorandom, label='NoRandom',marker='o')
plt.xlabel('S')
plt.ylabel('C')
plt.legend()
plt.savefig('SvC_otra_def2.png')
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    