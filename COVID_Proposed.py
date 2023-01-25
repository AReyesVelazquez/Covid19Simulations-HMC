#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 25 18:52:02 2022

@author: abrahamreyes
"""

import numpy             as np
import random            as rnd
import matplotlib.pyplot as plt
import time

#States and reactions

class State:
    def __init__(self,statename):
        self.statename  = str(statename)
        self.population = [init_pob[sts.index(statename)]]
        
class MonoReac:
    def __init__(self,A,B):
        self.a = sts.index(A)
        self.b = sts.index(B)
    def ford_ef(self):
        if init_pob[self.a] >= 1:
            init_pob[self.a] -= 1
            init_pob[self.b] += 1
    def back_ef(self):
        if init_pob[self.b] >= 1:
            init_pob[self.a] += 1
            init_pob[self.b] -= 1
    def ford_wgt(self):
        return init_pob[self.a]
    def back_wgt(self):
        if init_pob[self.a] > 0:
            return (init_pob[self.a]+1)
        else:
            return 0

class BinReac:
    def __init__(self,A,B,C):
        self.a = sts.index(A)
        self.b = sts.index(B)
        self.c = sts.index(C)
    def ford_ef(self):
        if init_pob[self.a] >= 1:
            init_pob[self.a] -= 1
            init_pob[self.c] += 1
    def back_ef(self):
        if init_pob[self.c] >= 1:
            init_pob[self.a] += 1
            init_pob[self.c] -= 1
    def ford_wgt(self):
        return init_pob[self.a]*init_pob[self.b]
    def back_wgt(self):
        if init_pob[self.a] > 0:
            return (init_pob[self.a]+1)*init_pob[self.b]
        else:
            return 0
        
def Simulation():
    # Variant triggers data
    week = 0
    soc_wgt = soc_rates[0]
    alph = 0
    delt = 0
    omic = 0
    # Times
    t        = 0
    t_reac   = 0.3
    t_stop   = 730
    timex    = [0]
    N        = np.sum(init_pob)
    # Data arrays
    st_num    = len(st_indx)
    mono_num  = len(mono_indx)
    bin_num   = len(bin_indx)
    States    = []
    Rates     = []
    Reactions = []
    Daychange = np.zeros(t_stop)
    for i in range(st_num):
        States.append(State(sts[i]))
    for i in range(mono_num):
        Rates.append(mono_prob[i]/mono_time[i])
        Reactions.append(MonoReac(mono_a[i], mono_b[i]))
    for i in range(bin_num):
        Rates.append(bin_prob[i]/(N*bin_time[i]))
        Reactions.append(BinReac(bin_a[i], bin_b[i], bin_c[i]))
    reac_num = len(Reactions)
    O = np.zeros(reac_num)
    P = np.zeros(reac_num)
    o = np.zeros(reac_num)
    p = np.zeros(reac_num)
    # Monte-Carlo Loop
    while t < t_stop:
        t_day = np.floor(t).astype(int)
        # Week counter
        if t_day % 7 == 0 and soc_indx[t_day] == 0:
            soc_wgt = soc_rates[week]
            soc_indx[t_day] = 1
            week += 1
            print(t_day)
            print(init_pob)
        # Daily births and vaccinations
        if Daychange[t_day] == 0:
            # Births 
            init_pob[0] += 28
            # Vaccinations start after t=329
            if t_day >= 329:
                vac = 3*(t_day-329)
                if init_pob[0] > vac:
                    init_pob[0] -= vac
                    init_pob[7] += vac
            Daychange[t_day] = 1
        # Variants triggers
        if t_day >= 344 and alph == 0:
            bin_prob[0]  = 3.627
            bin_prob[1]  = 0.063*3.627
            bin_prob[2]  = 0.15*3.627
            mono_prob[0] = 3.0
            alph = 1
        if t_day >= 430 and delt == 0:
            bin_prob[0]  = 5.08
            bin_prob[1]  = 0.12*5.08
            bin_prob[2]  = 0.15*5.08
            mono_prob[0] = 4.0
            delt = 1
        if t_day >= 700 and omic == 0:
            bin_prob[0]  = 16.205
            bin_prob[1]  = 0.3*16.205
            bin_prob[2]  = 0.81*16.205
            mono_prob[0] = 3.0
            omic = 1
        N = np.sum(init_pob)
        Rates[10] = (soc_wgt*bin_prob[0])/(N*bin_time[0])
        Rates[11] = (soc_wgt*bin_prob[1])/(N*bin_time[1])
        Rates[12] = (soc_wgt*bin_prob[2])/(N*bin_time[2])
        for i in range(reac_num):
            O[i] = Reactions[i].ford_wgt()
            P[i] = O[i]*Rates[i]
            o[i] = -Reactions[i].back_wgt()
            p[i] = o[i]*Rates[i]
        Q = np.sum(P)
        q = np.sum(p)
        if Q == 0:
            print('broken')
            print(init_pob)
            break
        tau = np.log(1/rnd.uniform(0,1))/Q
        Tau = np.log(1/rnd.uniform(0,1))/q
        if tau <= t_reac:
            uwu = rnd.uniform(0,1)*(Q/reac_num)
            m   = 0
            n   = 0
            while m < uwu:
                m += P[n]/reac_num
                n += 1
            nu = n-1
            Reactions[nu].ford_ef()
            t += tau
        if tau > t_reac:
            t = t
        for i in range(st_num):
            States[i].population.append(init_pob[i])
        timex.append(t)
        if tau <= t_reac and np.abs(Tau) <= tau:
            ewe = np.abs(rnd.uniform(0,1)*(q/reac_num))
            mm  = 0
            nn  = 0
            while mm < ewe:
                mm += np.abs(p[nn])/reac_num
                nn += 1
            nuu = nn-1
            Reactions[nuu].back_ef()
            t += Tau
        if tau <= t_reac and np.abs(Tau) > tau:
            t = t
        for i in range(st_num):
            States[i].population.append(init_pob[i])
        timex.append(t)

    print(init_pob)
    print("--- %s seconds ---" % (time.time()-start_time))
    
    fig, (ax1) = plt.subplots(1)
    fig.suptitle('IvT (Proposed Algorithm)')
    plt.xlabel('Days')
    plt.ylabel('Population')
    ax1.plot(timex,States[1].population,label=States[1].statename,c='c')
    ax1.plot(timex,States[2].population,label=States[2].statename,c='m')
    ax1.legend(loc='upper left')
    
    fig, (ax2) = plt.subplots(1)
    fig.suptitle('HvT (Proposed Algorithm)')
    plt.xlabel('Days')
    plt.ylabel('Population')
    ax2.plot(timex,States[3].population,label=States[2].statename,c='navy')
    ax2.plot(timex,States[4].population,label=States[4].statename,c='green')
    ax2.legend(loc='upper left')

    plt.show()


#Simulation variables
#States of the network
st_indx  = [   1   ,  2  ,  3  ,  4  ,  5  ,  6  ,  7  ,  8  ]
sts      = [  'S'  , 'E' , 'I' , 'H' ,'ICU', 'R' , 'D' , 'V' ]
init_pob = [999999 ,  0  ,  1  ,  0  ,  0  ,  0  ,  0  ,  0  ]

#Mono-particle transitions S_a -> S_b
mono_indx = [ 1   ,  2  ,  3  ,  4  ,  5  ,  6  ,  7  ,  8  ,  9  ,  10 ]
mono_a    = [ 'E' , 'I' , 'I' , 'I' , 'H' , 'H' ,'ICU','ICU', 'R' , 'V' ]
mono_b    = [ 'I' , 'H' ,'ICU', 'R' , 'R' , 'D' , 'R' , 'D' , 'S' , 'S' ]
mono_prob = [ 1.0 ,0.138,0.061,0.801,0.962,0.038,0.584,0.416, 1.0 , 1.0 ]
mono_time = [ 5.2 ,14.00,14.00,14.00,12.00,12.00,8.000,8.000,180.0,180.0]

#Bi-particle transitions S_a + S_b -> S_c + S_b
bin_indx = [  1  ,  2  ,  3  ]
bin_a    = [ 'S' , 'V' , 'R' ]
bin_b    = [ 'I' , 'I' , 'I' ]
bin_c    = [ 'E' , 'E' , 'E' ]
bin_prob = [2.790,0.05*2.790,0.15*2.790]
bin_time = [14.00,14.00,14.00]

start_time = time.time()

soc_rates = [1.0,1.0,1.0,1.0,1.001,1.03,1.06,1.05,0.94,0.71,
            0.65,0.64,0.63,0.66,0.69,0.74,0.76,0.77,0.78,0.82,
            0.84,0.85,0.85,0.85,0.82,0.84,0.83,0.83,0.84,0.83,
            0.83,0.83,0.83,0.81,0.83,0.84,0.84,0.84,0.83,0.83,
            0.83,0.84,0.82,0.82,0.80,0.77,0.82,0.82,0.82,0.73,
            0.77,0.78,0.77,0.77,0.77,0.78,0.75,0.80,0.84,0.85,
            0.86,0.87,0.87,0.87,0.88,0.88,0.88,0.90,0.91,0.90,
            0.91,0.89,0.92,0.92,0.90,0.90,0.88,0.90,0.91,0.91,
            0.92,0.91,0.91,0.90,0.90,0.88,0.91,0.90,0.91,0.92,
            0.91,0.92,0.91,0.92,0.91,0.91,0.89,0.86,0.92,0.93,
            0.92,0.78,0.81,0.84,0.82]

soc_indx = np.zeros(730)


Simulation()