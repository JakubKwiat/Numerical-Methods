#!/usr/bin/env python
# coding: utf-8

# In[29]:


import csv
import matplotlib.pyplot as plt
dane=[]
daty=[]
macd_tab=[]
signal_tab=[]
def odczyt(nazwa_pliku,daty,wartosci):
    with open(nazwa_pliku,newline='')as plikcsv:
        czytanie=csv.reader(plikcsv,delimiter=',')
        for wiersz in czytanie:
            daty.append(wiersz[0])
            wartosci.append(float(wiersz[4]))
    return wartosci,daty

def EMA(liczba_okresow,dzien,wartosci):
    alfa=2/(liczba_okresow+1)
    tab=[]
    for j in range(liczba_okresow+1):
        tab.append(wartosci[dzien-liczba_okresow+j-1])
    tab.reverse()
    p_0=float(tab[0])
    licznik=p_0
    mianownik=1
    for i in range(liczba_okresow):
        licznik += ((1-alfa)**(i+1))*tab[i+1]
        mianownik += (1-alfa)**(i+1)
    return (licznik/mianownik)

def MACD_obliczanie(tab,wartosci):
    for i in range(1000):
        if i >= 26:
            tab.append(EMA(12,i,wartosci)-EMA(26,i,wartosci))
        else:
            tab.append(0)
    return tab

def SIGNAL_obliczanie(signal_wartosci,macd_wartosci):
    for i in range(1000):
        if i>=35:
            signal_wartosci.append(EMA(9,i,macd_wartosci))
        else:
            signal_wartosci.append(0)
    return signal_wartosci
def Handlowanie_akcjami(signal_wartosci,macd_wartosci,wartosci,pieniadze):
    akcje=0
    for i in range(1,1000):
        if macd_wartosci[i-1]<=signal_wartosci[i-1] and macd_wartosci[i]>signal_wartosci[i] and pieniadze!=0:
            akcje=float(pieniadze/wartosci[i])
            pieniadze=0
        if macd_wartosci[i-1]>=signal_wartosci[i-1] and macd_wartosci[i]<signal_wartosci[i] and akcje!=0:
            pieniadze=float(akcje*wartosci[i])
            akcje=0
    if pieniadze==0:
        pieniadze=float((akcje*wartosci[999]/1000)*100-100)
    else:
        pieniadze=(pieniadze/1000)*100-100
    return pieniadze
def Handlowanie_akcjami_opoznione(signal_wartosci,macd_wartosci,wartosci,pieniadze):
    akcje=0
    for i in range(1,1000):
        if macd_wartosci[i-4]<=signal_wartosci[i] and macd_wartosci[i]>signal_wartosci[i] and pieniadze!=0:
            akcje=float(pieniadze/wartosci[i])
            pieniadze=0
        if macd_wartosci[i-4]>=signal_wartosci[i] and macd_wartosci[i]<signal_wartosci[i] and akcje!=0:
            pieniadze=float(akcje*wartosci[i])
            akcje=0
    if pieniadze==0:
        pieniadze=float((akcje*wartosci[999]/1000)*100-100)
    else:
        pieniadze=(pieniadze/1000)*100-100
    return pieniadze
dane,daty=odczyt('KO.csv',daty,dane)
miesiace=[]
for j in range(1000):
    if j%30==0:
        miesiace.append(daty[j])
macd_tab=MACD_obliczanie(macd_tab,dane)
signal_tab=SIGNAL_obliczanie(signal_tab,macd_tab)


# In[28]:


plt.figure(dpi=200)
plt.plot(daty,signal_tab,color='green',lw=0.2)
plt.plot(daty,macd_tab,color='red',lw=0.2)
plt.legend(['SIGNAL','MACD'])
plt.title('Wykres Wskaźnika MACD')
plt.xlabel('Czas')
plt.ylabel('Wartosc')
plt.xticks(miesiace,rotation=270)
plt.savefig('Wskaznik MACD')
plt.show()


# Wskaźnik MACD dla notowań spółki THE COCA-COLA Company w okresie od 3 kwietnia 2017 roku do 23 marca 2021 roku.

# In[27]:


plt.figure(dpi=200)
plt.plot(daty,signal_tab,color='green',lw=0.2)
plt.plot(daty,macd_tab,color='red',lw=0.2)
plt.plot(dane,color='blue',lw=0.2)
plt.legend(['SIGNAL','MACD','Cena Zamknięcia'])
plt.title('Wykres Wskaźnika MACD i Ceny Zamknięcia')
plt.xlabel('Czas')
plt.ylabel('Wartosc')
plt.xticks(miesiace,rotation=270)
plt.savefig('Wskaznik MACD I CLOSING PRICE')
plt.show()


# Notowanie akcji spółki THE COCA-COLA Company zestawione z wyznaczonym na ich podstawie wskaźnikiem MACD w okresie od 3 kwietnia 2017 roku do 23 marca 2021 roku. Na wykresie można zauważyć, że wskaźnik MACD odwzorowuje rzeczywiste zachowanie notowań akcji, lecz robi to z pewnym opóźnieniem, dlatego jest on zaliczany do wskaźników opóźnionych, które generują sygnały na podstawie już zaistniałych zmian.

# In[30]:


kapital=1000
zysk=Handlowanie_akcjami(signal_tab,macd_tab,dane,kapital)
print("Nasz zysk po 1000 próbkach wynosi:", zysk,"procent")


# Wyniki symulacji zakupu i sprzedaży akcji, przy zastosowaniu algorytmu dokonującego natychmiastowej sprzedaży lub zakupu akcji po każdym przecięciu się linii MACD i SIGNAL. Zysk przy użyciu tego algorytmu wyniósł niewiele ponad 4%.Jest to bardzo niski wynik i raczej mało opłacalny.

# In[31]:


kapital=1000
zysk=Handlowanie_akcjami_opoznione(signal_tab,macd_tab,dane,kapital)
print("Nasz zysk po 1000 próbkach wynosi:", zysk,"procent")


# Wyniki symulacji zakupu i sprzedaży akcji, przy zastosowaniu algorytmu dokonującego sprzedaży lub zakupu akcji z kilkudniowym
# opóźnieniem po przecięciu się lniii MACD i SIGNAL. Zysk przy użyciu tego algorytmu okazał się dużo większy niż w poprzednim przypadku i wynosi ponad 26% co przy odpowiednio duzym kapitale może powodować całkiem spore zyski z inwestycji.

# In[ ]:




