#!/usr/bin/env python
# coding: utf-8

# In[139]:


import csv
import matplotlib.pyplot as plt
import math


# In[140]:


def odczyt(sciezka):
    with open(sciezka,newline='')as plikcsv:
        czytanie=csv.reader(plikcsv,delimiter=',')
        dystanse=[]
        wysokosci=[]
        for wiersz in czytanie:
            if(wiersz[0]!='Dystans (m)' and wiersz[1]!='WysokoĹ›Ä‡ (m)'):
                dystanse.append(float(wiersz[0]))
                wysokosci.append(float(wiersz[1]))  
    return dystanse,wysokosci


# In[141]:


dystanse_ME,wysokosci_ME=odczyt('dane\MountEverest.csv')
dystanse_SG,wysokosci_SG=odczyt('dane\SpacerniakGdansk.csv')
dystanse_WKK,wysokosci_WKK=odczyt('dane\WielkiKanionKolorado.csv')


# In[142]:


def wybór_punktow(lista_x,lista_y,ile_punktow):
    nowa_lista_x=[]
    nowa_lista_y=[]
    odstep=round(len(lista_x)/ile_punktow-1)
    for i in range(ile_punktow):
        nowa_lista_x.append(lista_x[i*odstep])
        nowa_lista_y.append(lista_y[i*odstep])
        indeks=i*odstep
    return nowa_lista_x,nowa_lista_y,indeks
def obliczanie_fi(wartosci_x,j,x):
    wartosc=1
    for i in range(len(wartosci_x)):
        if i!=j:
            wartosc*=(x-wartosci_x[i])/(wartosci_x[j]-wartosci_x[i])
    return wartosc
def Interpolacja_Lagrange(wartosci_x,wartosci_y,x):
    wartosc=0
    for j in range(len(wartosci_x)):
        wartosc+=wartosci_y[j]*obliczanie_fi(wartosci_x,j,x)
    return wartosc
def Lagrange(wartosci_x,wartosci_y,ile_punktow,nazwa):
    interpolowane_wartosci=[]
    wartosci_x_pkt,wartosci_y_pkt,indeks_ostatni=wybór_punktow(wartosci_x,wartosci_y,ile_punktow)
    for i,x in enumerate(wartosci_x):
        if i>indeks_ostatni:
            break
        interpolowane_wartosci.append(Interpolacja_Lagrange(wartosci_x_pkt,wartosci_y_pkt,x))
    plt.figure(dpi=100)
    plt.plot(wartosci_x[:indeks_ostatni],wartosci_y[:indeks_ostatni])
    plt.plot(wartosci_x_pkt,wartosci_y_pkt,'o')
    plt.plot(wartosci_x[:indeks_ostatni+1],interpolowane_wartosci)
    plt.xlabel("dystans [m]")
    plt.ylabel("wysokość [m]")
    plt.title('Metoda Lagranga dla: ' + str(ile_punktow) + ' punktów - ' + str(nazwa))
    plt.show()


# In[125]:


Lagrange(dystanse_ME,wysokosci_ME,6,'Mount Everest')
Lagrange(dystanse_ME,wysokosci_ME,9,'Mount Everest')
Lagrange(dystanse_ME,wysokosci_ME,15,'Mount Everest')


# In[143]:


Lagrange(dystanse_SG,wysokosci_SG,6,'Spacerniak Gdańsk')
Lagrange(dystanse_SG,wysokosci_SG,9,'Spacerniak Gdańsk')
Lagrange(dystanse_SG,wysokosci_SG,15,'Spacerniak Gdańsk')


# In[145]:


Lagrange(dystanse_WKK,wysokosci_WKK,6,'Wielki Kanion Kolorado')
Lagrange(dystanse_WKK,wysokosci_WKK,9,'Wielki Kanion Kolorado')
Lagrange(dystanse_WKK,wysokosci_WKK,15,'Wielki Kanion Kolorado')


# In[146]:


def pivoting(i,U,L,P,N):
        pivot = abs(U[i][i])
        pivot_ind = i
        for j in range(i+1, len(U)):
            if abs(U[j][i] > pivot):
                pivot = U[j][i]
                pivot_ind = j
        if pivot_ind != i:
            for j in range(0, len(U)):
                if j >= i:
                    temp = U[i][j]
                    U[i][j] = U[pivot_ind][j]
                    U[pivot_ind][j] = temp
                else:
                    temp = L[i][j]
                    L[i][j] = L[pivot_ind][j]
                    L[pivot_ind][j] = temp
                temp = P[i][j] 
                P[i][j] = P[pivot_ind][j]
                P[pivot_ind][j] = temp
def faktoryzacja(b,U,L,P,N,A):
        for i in range(N):
            for j in range(N):
                U[i][j] = A[i][j]
        for i in range(N):
            for j in range(N):
                if j == i:
                    L[i][j] = 1
                    P[i][j] = 1
                else:
                    L[i][j] = 0
                    P[i][j] = 0
        for k in range(N - 1):
            pivoting(k,U,L,P,N)
            for j in range(k + 1,N):
                L[j][k] = U[j][k] / U[k][k]
                for p in range(k, N):
                    U[j][p] = U[j][p] - (L[j][k] * U[k][p])
        y = [0 for i in range(N)]
        x = [0 for i in range(N)]
        for i in range(N):
            alfa =b[i]
            for j in range(i):
                alfa = alfa - (L[i][j] * y[j])
            y[i] = alfa / L[i][i]
        for i in range(N-1, -1, -1):
            alfa = y[i]
            for j in range(i + 1,N):
                alfa = alfa - (U[i][j] * x[j])
            x[i] = alfa /U[i][i]
        return x
def Interpolacja_funkcjami_sklejanymi(wartosci_x,wartosci_y,ile_punktow,nazwa):
    wartosci_x_pkt,wartosci_y_pkt,indeks_ostatni=wybór_punktow(wartosci_x,wartosci_y,ile_punktow)
    N=4*(len(wartosci_x_pkt)-1)
    U=[[0 for i in range(N)] for i in range(N)]
    L=[[0 for i in range(N)] for i in range(N)]
    P=[[0 for i in range(N)] for i in range(N)]
    A=[[0 for i in range(N)] for i in range(N)]
    b=[0 for i in range(N)]
    for i in range(N):
        A[0][0]=1
        b[0]=wartosci_y_pkt[0]
        h=wartosci_x_pkt[1]-wartosci_x_pkt[0]
        A[1][0]=1
        A[1][1]=h
        A[1][2]=h**2
        A[1][3]=h**3
        b[1]=wartosci_y_pkt[1]
        A[2][2]=1
        b[2]=0
        h=wartosci_x_pkt[len(wartosci_x_pkt)-1]-wartosci_x_pkt[len(wartosci_x_pkt)-2]
        A[3][4*(len(wartosci_x_pkt)-2)+2]=2
        A[3][4*(len(wartosci_x_pkt)-2)+3]=6*h
        b[3]=0
        for i in range(1,len(wartosci_x_pkt)-1):
            h=wartosci_x_pkt[i]-wartosci_x_pkt[i-1]
            A[4*i][4*i]=1
            b[4*i]=wartosci_y_pkt[i]
            A[4*i+1][4*i] = 1
            A[4*i+1][4*i+1] = h
            A[4*i+1][4*i+2] = h**2
            A[4*i+1][4*i+3] = h**3
            b[4*i+1] = wartosci_y_pkt[i+1]
            A[4*i+2][4*(i-1)+1] = 1
            A[4 * i + 2][4 * (i - 1) + 2] = 2 * h
            A[4 * i + 2][4 * (i - 1) + 3] = 3 * (h**2)
            A[4 * i + 2][4 * i + 1] = -1
            b[4 * i + 2] = 0
            A[4 * i + 3][4 * (i - 1) + 2] = 2
            A[4 * i + 3][4 * (i - 1) + 3] = 6*h
            A[4 * i + 3][4 * i + 2] = -2
            b[4 * i + 3] = 0
        wynik=faktoryzacja(b,U,L,P,N,A)
        Sklejanie=[]
        for i in range(indeks_ostatni):
            przedzial = int(i/indeks_ostatni * (len(wartosci_x_pkt)-1))
            Sklejanie.append(wynik[przedzial*4] + wynik[przedzial*4 + 1] * (wartosci_x[i] - wartosci_x_pkt[przedzial]) + wynik[przedzial*4 + 2] * (wartosci_x[i] - wartosci_x_pkt[przedzial]) ** 2 + wynik[przedzial*4 + 3] * (wartosci_x[i] - wartosci_x_pkt[przedzial]) ** 3)
    plt.figure(dpi=100)
    plt.plot(wartosci_x[:indeks_ostatni], Sklejanie[:indeks_ostatni])
    plt.plot(wartosci_x[:indeks_ostatni], wartosci_y[:indeks_ostatni])
    plt.plot(wartosci_x_pkt, wartosci_y_pkt,'o')
    plt.title('Interpolacja funkcjami sklejanymi dla: ' + str(ile_punktow) + ' punktów - ' + str(nazwa))
    plt.xlabel('dystans [m]')
    plt.ylabel('wysokość [m]')
    plt.show()


# In[ ]:


Interpolacja_funkcjami_sklejanymi(dystanse_ME,wysokosci_ME,10,'Mount Everest')
Interpolacja_funkcjami_sklejanymi(dystanse_ME,wysokosci_ME,18,'Mount Everest')
Interpolacja_funkcjami_sklejanymi(dystanse_ME,wysokosci_ME,40,'Mount Everest')


# In[137]:


Interpolacja_funkcjami_sklejanymi(dystanse_SG,wysokosci_SG,10,'Spacerniak Gdańsk')
Interpolacja_funkcjami_sklejanymi(dystanse_SG,wysokosci_SG,18,'Spacerniak Gdańsk')
Interpolacja_funkcjami_sklejanymi(dystanse_SG,wysokosci_SG,40,'Spacerniak Gdańsk')


# In[138]:


Interpolacja_funkcjami_sklejanymi(dystanse_WKK,wysokosci_WKK,10,'Wielki Kanion Kolorado')
Interpolacja_funkcjami_sklejanymi(dystanse_WKK,wysokosci_WKK,18,'Wielki Kanion Kolorado')
Interpolacja_funkcjami_sklejanymi(dystanse_WKK,wysokosci_WKK,40,'Wielki Kanion Kolorado')


# In[ ]:




