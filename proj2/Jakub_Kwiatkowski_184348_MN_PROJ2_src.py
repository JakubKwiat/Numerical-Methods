#!/usr/bin/env python
# coding: utf-8

# In[35]:


import math
import time
import matplotlib.pyplot as plt
from matplotlib.pyplot import figure


# In[45]:


#A
# index = 184348
a1=8
a2=-1
a3=-1
f=4
N=948
b=[math.sin(n*(f+1)) for n in range(N)]
def stworz_macierz_A(N,a1,a2,a3):
    macierz=[[0 for i in range(N)] for i in range(N)]
    for i in range(N):
        macierz[i][i]=a1
        if i<N-1:
            macierz[i][i+1]=a2
        if i>0:
            macierz [i][i-1]=a2
        if i>1:
            macierz [i][i-2]=a3
        if i<N-2:
            macierz [i][i+2]=a3
    return macierz
A=stworz_macierz_A(N,a1,a2,a3)


# In[6]:


def norma_euklidesowa_wektora_residuum(wektor):
    N=len(wektor)
    norma=0
    for i in range(N):
        norma = norma + wektor[i]**2
    return math.sqrt(norma)


# In[7]:


def wektor_residuum(A,x,b):
    N=len(x)
    res=[0 for i in range(N)]
    for i in range(N):
        for j in range(N):
            res[i]=res[i]+A[i][j]*x[j]
    res = [(res[i]-b[i]) for i in range(N)]
    return res


# In[27]:


#B - metoda Jacobiego
bariera=10**-9
def metoda_Jacobiego(A,b,bariera):
    start=time.time()
    N=len(A)
    licznik_iteracji=0
    x=[1 for i in range(N)]
    x_poprzednie = [1 for i in range(N)]
    res = [1 for i in range(N)]
    while norma_euklidesowa_wektora_residuum(res) > bariera:
        for i in range(N):
            sumaL=0
            sumaU=0
            for j in range(i):
                sumaL=sumaL+A[i][j]*x_poprzednie[j]
            for j in range(i+1,N):
                sumaU=sumaU+A[i][j]*x_poprzednie[j]
            x[i]=(b[i]-sumaL-sumaU)/A[i][i]
        x_poprzednie=[x[i] for i in range(N)]
        res=wektor_residuum(A,x,b)
        licznik_iteracji+=1
    koniec=time.time()
    czas_trwania=koniec-start
    return x,licznik_iteracji,czas_trwania
Jacobi_wyniki,liczba_iteracji_Jacobi,czas_trwania_Jacobi=metoda_Jacobiego(A,b,bariera)
print("Metoda Jacobiego")
print("Liczba iteracji:",liczba_iteracji_Jacobi)
print("Czas trwania:",czas_trwania_Jacobi)


# In[28]:


#B - metoda Gaussa-Seidla
def metoda_Gaussa_Seidla(A,b,bariera):
    start=time.time()
    N=len(A)
    licznik_iteracji=0
    x=[1 for i in range(N)]
    x_poprzednie = [1 for i in range(N)]
    res = [1 for i in range(N)]
    while norma_euklidesowa_wektora_residuum(res) > bariera:
        for i in range(N):
            sumaL=0
            sumaU=0
            for j in range(i):
                sumaL=sumaL+A[i][j]*x[j]
            for j in range(i+1,N):
                sumaU=sumaU+A[i][j]*x_poprzednie[j]
            x[i]=(b[i]-sumaL-sumaU)/A[i][i]
        x_poprzednie=[x[i] for i in range(N)]
        res=wektor_residuum(A,x,b)
        licznik_iteracji+=1
    koniec=time.time()
    czas_trwania=koniec-start
    return x,licznik_iteracji,czas_trwania
Gauss_Seidl_wyniki,liczba_iteracji_Gauss_Seidl,czas_trwania_Gauss_Seidl=metoda_Gaussa_Seidla(A,b,bariera)
print("Metoda Gaussa-Seidla")
print("Liczba iteracji:",liczba_iteracji_Gauss_Seidl)
print("Czas trwania:",czas_trwania_Gauss_Seidl)


# Metoda Gaussa-Seidla jest około 1,5 razy szybsza niż metoda Jacobiego jeśli chodzi o czas trwania oraz potrzebuje około 1,5 razy mniej iteracji.

# In[10]:


#C
a1=3
a2=-1
a3=-1
N=948
A=stworz_macierz_A(N,a1,a2,a3)
Gauss_Seidl_wyniki,liczba_iteracji_Gauss_Seidl,czas_trwania_Gauss_Seidl=metoda_Gaussa_Seidla(A,b,bariera)
print("Metoda Gaussa-Seidla")
print("Liczba iteracji:",liczba_iteracji_Gauss_Seidl)
print("Czas trwania:",czas_trwania_Gauss_Seidl)


# In[11]:


Jacobi_wyniki,liczba_iteracji_Jacobi,czas_trwania_Jacobi=metoda_Jacobiego(A,b,bariera)
print("Metoda Jacobiego")
print("Liczba iteracji:",liczba_iteracji_Jacobi)
print("Czas trwania:",czas_trwania_Jacobi)


# Przy następujących parametrach podczas liczenia normy euklidesowej wektora residuum pojawia się błąd przy obu metodach. Wnioskiem z tego jest, że przy tych parametrach metody iteracyjne nie zbiegają się

# In[40]:


#D
a1=8
a2=-1
a3=-1
f=4
N=948
b = [math.sin(i * (f + 1)) for i in range(N)]
A=stworz_macierz_A(N,a1,a2,a3)
def tworzenie_macierzy_LU(A):
    N=len(A)
    U=[[0 for i in range(N)] for i in range(N)]
    for i in range (N):
        for j in range(N):
            U[i][j]=A[i][j]
    L=[[0 for i in range(N)] for i in range(N)]
    for i in range(N):
        L[i][i]=1
    for k in range(N-1):
        for j in range(k+1,N):
            L[j][k]=U[j][k]/U[k][k]
            for i in range(k,N):
                U[j][i]=U[j][i]-L[j][k]*U[k][i]
    return L,U
def faktoryzacja_LU(A,b):
    start = time.time()
    N=len(A)
    x = [0 for i in range(N)]
    y = [0 for i in range(N)]
    L, U = tworzenie_macierzy_LU(A)
    for i in range(N):
        sumaL = 0
        for j in range(i):
            sumaL = sumaL + L[i][j]*y[j]
        y[i] = (b[i] - sumaL)/L[i][i]
    for i in range(N-1, -1, -1):
        sumaU = 0
        for j in range(i+1, N):
            sumaU = sumaU + U[i][j] * x[j]
        x[i] = (y[i] - sumaU)/U[i][i]
    res = wektor_residuum(A, x, b)
    norma=norma_euklidesowa_wektora_residuum(res)
    koniec = time.time()
    czas_trwania=koniec-start
    return x, y, czas_trwania,norma
x,y,czas_trwania_faktoryzacji_LU,norma_residuum_faktoryzacji_LU=faktoryzacja_LU(A,b)
print("Faktoryzacja LU")
print("Norma residuum",norma_residuum_faktoryzacji_LU)


# Norma residuum w przypadku metody bezpośredniej - metody faktoryzacji LU jest rzędu 10^(-15) czyli bardzo bliska zeru, co oznacza dużą dokładność obliczeń

# In[32]:


#E
czas_trwania_faktoryzacja_LU_tab=[]
czas_trwania_Gauss_Seidl_tab=[]
czas_trwania_Jacobi_tab=[]
a1=8
a2=-1
a3=-1
f=4
N=100
b=[math.sin(i * (f + 1)) for i in range(N)]
A=stworz_macierz_A(N,a1,a2,a3)
x,y,czas_trwania_faktoryzacji_LU,norma_residuum_faktoryzacji_LU=faktoryzacja_LU(A,b)
czas_trwania_faktoryzacja_LU_tab.append(czas_trwania_faktoryzacji_LU)
Gauss_Seidl_wyniki,liczba_iteracji_Gauss_Seidl,czas_trwania_Gauss_Seidl=metoda_Gaussa_Seidla(A,b,bariera)
czas_trwania_Gauss_Seidl_tab.append(czas_trwania_Gauss_Seidl)
Jacobi_wyniki,liczba_iteracji_Jacobi,czas_trwania_Jacobi=metoda_Jacobiego(A,b,bariera)
czas_trwania_Jacobi_tab.append(czas_trwania_Jacobi)
N=500
b=[math.sin(i * (f + 1)) for i in range(N)]
A=stworz_macierz_A(N,a1,a2,a3)
x,y,czas_trwania_faktoryzacji_LU,norma_residuum_faktoryzacji_LU=faktoryzacja_LU(A,b)
czas_trwania_faktoryzacja_LU_tab.append(czas_trwania_faktoryzacji_LU)
Gauss_Seidl_wyniki,liczba_iteracji_Gauss_Seidl,czas_trwania_Gauss_Seidl=metoda_Gaussa_Seidla(A,b,bariera)
czas_trwania_Gauss_Seidl_tab.append(czas_trwania_Gauss_Seidl)
Jacobi_wyniki,liczba_iteracji_Jacobi,czas_trwania_Jacobi=metoda_Jacobiego(A,b,bariera)
czas_trwania_Jacobi_tab.append(czas_trwania_Jacobi)
N=1000
b=[math.sin(i * (f + 1)) for i in range(N)]
A=stworz_macierz_A(N,a1,a2,a3)
x,y,czas_trwania_faktoryzacji_LU,norma_residuum_faktoryzacji_LU=faktoryzacja_LU(A,b)
czas_trwania_faktoryzacja_LU_tab.append(czas_trwania_faktoryzacji_LU)
Gauss_Seidl_wyniki,liczba_iteracji_Gauss_Seidl,czas_trwania_Gauss_Seidl=metoda_Gaussa_Seidla(A,b,bariera)
czas_trwania_Gauss_Seidl_tab.append(czas_trwania_Gauss_Seidl)
Jacobi_wyniki,liczba_iteracji_Jacobi,czas_trwania_Jacobi=metoda_Jacobiego(A,b,bariera)
czas_trwania_Jacobi_tab.append(czas_trwania_Jacobi)
N=2000
b=[math.sin(i * (f + 1)) for i in range(N)]
A=stworz_macierz_A(N,a1,a2,a3)
x,y,czas_trwania_faktoryzacji_LU,norma_residuum_faktoryzacji_LU=faktoryzacja_LU(A,b)
czas_trwania_faktoryzacja_LU_tab.append(czas_trwania_faktoryzacji_LU)
Gauss_Seidl_wyniki,liczba_iteracji_Gauss_Seidl,czas_trwania_Gauss_Seidl=metoda_Gaussa_Seidla(A,b,bariera)
czas_trwania_Gauss_Seidl_tab.append(czas_trwania_Gauss_Seidl)
Jacobi_wyniki,liczba_iteracji_Jacobi,czas_trwania_Jacobi=metoda_Jacobiego(A,b,bariera)
czas_trwania_Jacobi_tab.append(czas_trwania_Jacobi)
N=3000
b=[math.sin(i * (f + 1)) for i in range(N)]
A=stworz_macierz_A(N,a1,a2,a3)
x,y,czas_trwania_faktoryzacji_LU,norma_residuum_faktoryzacji_LU=faktoryzacja_LU(A,b)
czas_trwania_faktoryzacja_LU_tab.append(czas_trwania_faktoryzacji_LU)
Gauss_Seidl_wyniki,liczba_iteracji_Gauss_Seidl,czas_trwania_Gauss_Seidl=metoda_Gaussa_Seidla(A,b,bariera)
czas_trwania_Gauss_Seidl_tab.append(czas_trwania_Gauss_Seidl)
Jacobi_wyniki,liczba_iteracji_Jacobi,czas_trwania_Jacobi=metoda_Jacobiego(A,b,bariera)
czas_trwania_Jacobi_tab.append(czas_trwania_Jacobi)


# In[44]:


#E
argumenty=[100,500,1000,2000,3000]
plt.figure(dpi=200)
plt.plot(argumenty,czas_trwania_Jacobi_tab,label='Metoda Jacobiego')
plt.plot(argumenty,czas_trwania_Gauss_Seidl_tab,label='Metoda Gaussa-Seidla')
plt.plot(argumenty,czas_trwania_faktoryzacja_LU_tab,label='Faktoryzacja LU')
plt.xlabel('Liczba N')
plt.ylabel('Czas trwania[s]')
plt.title('Zależność czasu trwania poszczególnych algorytmów w zależnośći od N')
plt.legend()
plt.show()

