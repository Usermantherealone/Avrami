# -*- coding: utf-8 -*-
"""
Created on Fri Jan 26 08:58:32 2024

@author: Nick-Alexander Merker
"""

# ========================================================================
# Import der nötigen Bibliotheken und Funktionen
# ========================================================================

import numpy as np
import scipy as sc
from scipy.signal import argrelmin
from scipy.signal import argrelmax
import pandas as pd
import os
import math

# ========================================================================
# Diese Funktion übersetzt die Zahl aus dem Dateinamen 
# in einen Füllstoffgehalt
# ========================================================================

def comparator(key):
    x = key.split('-')
    x = x[-1]
    
    if x == '1':
        j = 0.999
        return j
    
    elif x == '2':
        j = 0.995
        return j
    
    elif x == '3':
        j = 0.99
        return j
    
    elif x == '4':
        j = 0.97
        return j
    
    elif x == '5':
        j = 0.95
        return j
    
    elif x == '10':
        j = 0.90
        return j
    
    elif x == '15':
        j = 0.85
        return j
    
    elif x == '20':
        j = 0.80
        return j
    
    elif x == '25':
        j = 0.75
        return j
    
    elif x == '30':
        j = 0.70
        return j

# ========================================================================
# Hier beginnt das eigentliche Skript zur Anwendung der Avrami-Gleichung
# ========================================================================

# Hier werden der Input- und Output-Pfad angegeben
pfad: str = r"F:\07_Python\Avrami-Input"
pfadOUT: str = r"F:\07_Python\Avrami-Output"

#Listet die Dateien am Dateipfad in einer list auf
dateien = os.listdir(pfad)

# Hier werden einige Listen [] definiert, in die später wichtige Variablen 
# gespeichert werden
Avrami = pd.DataFrame()
Name = []
Kurve = []
Füllstoffanteil = []
Temperatur = []
n = []
K = []
R = []
v = []
induction = []
taulist = []
fitpoints = []
Bericht = []

#Schleife, die über die Dateien im Dateiordner iteriert
for datei in dateien:  
    
    #Importiert die DSC-Datei mit Pandas.
    df = pd.read_fwf(os.sep.join([pfad, datei]), widths=[15, 15, 15, 15, 15], header=None, skiprows=4, skipfooter=2, decimal=',', encoding='ANSI')
    
    #Teilt die Liste an der Stelle mit Inhalt 'Curve Name:'
    df_list = np.split(df, df[df[0]=='Kurvenname:'].index)   
    df_list.pop(0)
    
    # Schleife, die über die Experimentsegmente iteriert
    for i in df_list:
        # An dieser Stelle werden die wichtigen Information über das 
        # Experiment in der Variable i gespeichert
        Name = i.iat[1, 0] + i.iat[1, 1]
        Name = Name.split('-')
        Name.pop(0)
        sep = '-'
        Kurve.append(sep.join(Name))
        
        new_data = pd.DataFrame()
        new_data[0] = i.iloc[5:, 1]
        new_data[1] = pd.to_numeric(i.iloc[5:, 4])
        new_data.reset_index(drop=True, inplace=True)
        
        # Hier wird die Funktion vom Anfang aufgerufen und die Zahl in einen 
        # Füllstoffgehalt übersetzt
        key = Name[2]
        j = comparator(key)
        for g, h in enumerate(new_data.iloc[:, 1]):
            # An dieser Stelle wird der Wärmefluss durch den Massenanteil der 
            # Probe geteilt um sie für den Füllstoffgehalt zu korrigieren 
            new_data.iloc[g, 1] = h/j
               
        Füllstoffanteil.append(100-(j*100))
        Temperatur.append(i.iat[6, 3])
        temperature = str(i.iat[6, 3])
        new_data.dropna(inplace=True)
        i = new_data

        # Hier werden nochmal Listen definiert, die nur in dieser 
        # Schleife vorkommen sollen        
        k = []
        tau = []
        Fehler = []
        m = []
        fit_datax = []
        fit_datay = []

        # argrelmin() findet das lokale Minimum der DSC-Kurve und 
        # speichert es in der Variable k ab
        try:        
            k = argrelmin(i.iloc[:, 1].to_numpy(), order=10)
            if k[0][0] > 0:
                i = i[k[0][0]:].reset_index(drop=True)
                induction.append(k[0][0])
            
                # argrelmax() findet das lokale Maximum der DSC-Kurve und 
                # speichert es in der Variable tau ab. tau*3 wird als das 
                # zweite Integrationslimit festgelegt
                try:
                    tau = argrelmax(i.iloc[:, 1].to_numpy(), order=10)
                    if tau[0][0] > 0:
                        tau3 = tau[0][0]*3
                        taulist.append(tau[0][0])
                    
                    else:
                        print('error2')
                        Fehler.append('Kein Maximum')
                        taulist.append(0)
                        
                    if tau3 > len(i):
                        print('error')
                        Fehler.append('Experiment zu kurz')
                except:
                    print('error3')
                    Fehler.append('Kein Maximum')
                    taulist.append(0)       
                
                # An dieser Stelle wird die Avrami-Gleichung angewendet
                try:   
                    #Rechnet die Summe der Werte in Wärmefluss aus
                    htotal = i.iloc[:tau3+1, 1].sum()
                    
                    #Schleife über die Wärmefluss Zeilen
                    for x, y in enumerate(i.iloc[:tau3, 1]):
                        
                        #Rechnet die Avrami-Gleichung aus
                        deltah = i.iloc[:x+1, 1].sum()
                        try:
                            z = math.log(math.log(1/(1-(deltah/htotal))))
                            i.loc[x, 'Avrami'] = z
                        except: print('Der Fehler tritt auf bei:', deltah/htotal, i.iloc[:x+1, 4].sum(), htotal)
                
                        #Wenn z zwischen 0.03 und 0.20 liegt, wird dieser Wert für den Fit kopiert
                        if (z >= np.log(0.03) and z <= np.log(0.20) and x > 0):
                            fit_datax.append(np.log(x))
                            fit_datay.append(z)
                except:
                    print('Rechnung fehlgeschlagen')
                    Fehler.append('Rechnung fehlgeschlagen')
            else:
                print('Kein Minimum')
                Fehler.append('Kein Minimum')
                induction.append(0)
        except:
            print('Kein Minimum1')
            Fehler.append('Kein Minimum')
            induction.append(0)
         
        fitpoints.append(len(fit_datay))   
        
        # sc.stats.linregress() führt die lineare Regression durch und 
        # speichert die Ergebnisse in m
        try: 
            m = sc.stats.linregress(np.array(fit_datax), np.array(fit_datay))
            fit = []
            for s in fit_datax:
                fit.append(m[0]*s + m[1])
            
            # Die Ergebnisse aus m werden ggf. umgerechnet und für die n- und 
            # K-Parameter sowie das Pearson-R gespeichert
            n.append(m[0])
            K.append(np.exp(m[1]))
            R.append(m[2])
            v.append(np.exp(m[1])**-m[0])
        except:
            n.append(0)
            K.append(0)
            R.append(0)
            v.append(0)
        
        # Während des vorherigen Teils wurden auftretende Fehler abgefangen 
        # und gespeichert. Hier werden die Fehler zusammengeführt in Bericht gespeichert
        Fehlerstr = ' '
        Fehlerstr = Fehlerstr.join(Fehler)        
        Bericht.append(Fehlerstr)
        if Fehlerstr == '':
            print('okay')
        else:
            print(Fehlerstr)           
   
# ========================================================================
# Zusammenfassung der Daten zu einer Tabelle und Export in csv-Format
# ========================================================================
    
Temperatur = pd.DataFrame(Temperatur, columns= [['Temperatur'], ['°C'], [datei[:-4]]])
n = pd.DataFrame(n, columns= [['n'], ['[a. u.]'], [datei[:-4]]])
K = pd.DataFrame(K, columns= [['K'], ['[s^-n]'], [datei[:-4]]])
R = pd.DataFrame(R, columns= [['Pearson R'], [' '], [datei[:-4]]])
v = pd.DataFrame(v, columns= [['v'], ['1/s'], [' ']])
induction = pd.DataFrame(induction, columns= [['Induktionszeit'], ['s'], [' ']])
taulist = pd.DataFrame(taulist, columns= [['Halbwertszeit'], ['s'], [' ']])
fitpoints = pd.DataFrame(fitpoints, columns= [['Datenpunkte für Fit'], [' '], [' ']])
Bericht = pd.DataFrame(Bericht, columns= [['Bericht'], [' '], [' ']])
Kurve = pd.DataFrame(Kurve, columns= [['Kurve'], [' '], [' ']])
Füllstoffanteil = pd.DataFrame(Füllstoffanteil, columns= [['Füllstoffanteil'], ['wt%'], [' ']])
Avrami = pd.concat([Avrami, Kurve, Füllstoffanteil, Temperatur, n, K, R, v, induction, taulist, fitpoints, Bericht], axis=1)                
    
Avrami.to_csv(path_or_buf=os.path.join(pfadOUT, 'Auswertung-Avrami-Füllstoffe.csv'), sep=';', encoding='latin-1', index=False)