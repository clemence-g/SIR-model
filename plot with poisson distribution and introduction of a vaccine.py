# -*- coding: utf-8 -*-
"""
Created on Sun Dec 18 12:17:03 2022

@author: User
"""


import numpy as np
import matplotlib.pyplot as plt

def poisson(S, I, R, N, p1p2, beta, gamma, mu):
    '''
    Returns a vector with random values following the 
    Poisson distribution with lambda determined by the 
    differential equations.

    '''
    # Defining the reduced term beta_star with the effects of vaccine
    # p1 is the proportion of the vaccinated population
    # p2 is the efficacity of the vaccine
    beta_star = (1-p1p2)*beta
    
    beta_term = (beta_star*I*S)/N
    gamma_term = gamma*I
    mu_term = mu*I
    
    # Setting security limits
    if beta_term < 0:
        beta_term = 0
        
    if gamma_term < 0:
        gamma_term = 0
        
    if mu_term < 0:
        mu_term = 0
        
    # Defining the Poisson distribution 
    
    i_0 = np.random.poisson(beta_term, None)
    i_1 = np.random.poisson(gamma_term, None)
    i_2 = np.random.poisson(mu_term, None)
           
    return np.array([i_0, i_1, i_2])

def evolution(poisson, p1p2, beta, gamma, mu, t0, tf, h, N0, I0):
    '''
    Describes the evolution of the disease
    by using the Poisson distribution.
    
    '''
    # Initial values
    s = N0-I0 
    i = I0
    r = 0
    n = N0
    
    S1 = []
    I1 = []
    R1 = []
    N1 = []
    T1 = []
    
    t = t0
    
    while t < tf:
       
        s_last = s
        i_last = i
        r_last = r
        n_last = n
        
        # Calling the poisson function with following set arguments
        pois = poisson(s_last,i_last,r_last,n_last,p1p2,beta,gamma,mu) 
        
        s += -pois[0]
        i += pois[0] -pois[1] -pois[2]
        r += pois[1]
        n += -pois[2]
        
       
        S1.append(s)
        I1.append(i)
        R1.append(r)
        N1.append(n)
        T1.append(t)
        
        t += h
        
    return np.array([S1, I1, R1, N1, T1])

def graph(p1p2, beta, gamma, mu, k):
    '''
    Creates our plot and represents the evolution
    of the disease through time.

    '''
    # Calling the evolution function with following set arguments
    sol = evolution(poisson, p1p2, beta, gamma, mu, t0=0, tf=100, h=1, N0=1000, I0=10)
    
    if beta==0.8:
         contagion = " et très contagieuse"
    if beta == 0.1:
         contagion = " et peu contagieuse"
    if beta != 0.8 and beta !=0.1:
        contagion =""
    if mu == 0:
         mortalite = "non-mortelle"
    if mu == 0.02:
         mortalite= "mortelle"    
    if mu == 0.1:
         mortalite = "très mortelle"
    if mu == 0.5:
         mortalite = "fatale"
    if mu == 0.01:
         mortalite = "peu mortelle"
    if gamma == 0.5:
         recuperation = " avec récupération rapide"
    if gamma != 0.5:
        recuperation = ""
        
    plt.figure(k) # In order to have different plots for the different conditions : k = 0, 1, 2, etc.
    plt.plot(sol[4],sol[0], color="blue", label="S")
    plt.plot(sol[4],sol[1], color= "red", label="I")
    plt.plot(sol[4],sol[2], color="green", label="R")
    plt.plot(sol[4],sol[3], color="black", label="N")
    plt.legend()
    plt.xlabel("Temps") 
    plt.ylabel("Population")
    plt.title("beta = {}, gamma = {}, mu = {}, p1p2={}".format(beta, gamma, mu, p1p2), fontsize=10) #sous-titre
    plt.suptitle("Modèle SIR pour une maladie {}{}{}".format(mortalite, contagion, recuperation), fontsize= 15)
    plt.grid()
    plt.show()
           
graph(p1p2=0.5, beta=0.5, gamma=0.1, mu=0.02, k=0) 
graph(p1p2=0.7, beta=0.5, gamma=0.1, mu=0.02, k=1)
graph(p1p2=0.8, beta=0.5, gamma=0.1, mu=0.02, k=2)
graph(p1p2=0.9, beta=0.5, gamma=0.1, mu=0.02, k=3)

# First type of disease    
graph(p1p2=0.7, beta=0.5, gamma=0.1, mu=0.02, k=5)

# Second type      
graph(p1p2=0, beta=0.8, gamma=0.1, mu=0.5, k=3)

# Third type      
graph(p1p2=0, beta=0.8, gamma=0.1, mu=0.01, k=4)

# Fourth type      
graph(p1p2=0, beta=0.1, gamma=0.1, mu=0.01, k=5)

# Fifth type      
graph(p1p2=0, beta=0.1, gamma=0.5, mu=0.01, k=6)
