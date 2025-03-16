# -*- coding: utf-8 -*-
"""
Created on Mon Nov 28 09:07:37 2022

@author: denis & clem
"""
import numpy as np
import matplotlib.pyplot as plt

def derivatives(S, I, R, N, beta, gamma, mu):
    '''
    Returns a vector with the derivatives 
    of S, I, R and N.
    
    '''
    
    dS_dt= -beta*I*S/N
    dI_dt= (beta*I*S/N) - (gamma*I) - (mu*I)
    dR_dt= gamma*I
    dN_dt= -mu*I
    
    return np.array([dS_dt, dI_dt, dR_dt, dN_dt])


def euler(derivatives, beta, gamma, mu, t0, tf, h, N0, I0):
    '''
    Using Euler's method with a step of h 
    for approximating the evolution over time.
    
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
    
    while t < tf :
        
        s_last = s
        i_last = i
        r_last = r
        n_last = n
        
        s += h * derivatives(s_last,i_last,r_last,n_last, beta,gamma,mu)[0]
        i += h * derivatives(s_last,i_last,r_last,n_last, beta,gamma,mu)[1]
        r += h * derivatives(s_last,i_last,r_last,n_last, beta,gamma,mu)[2]
        n += h * derivatives(s_last,i_last,r_last,n_last, beta,gamma,mu)[3]
        
        S1.append(s)
        I1.append(i)
        R1.append(r)
        N1.append(n)
        T1.append(t)
        
        t += h
        
    
    return np.array([S1, I1, R1, N1, T1])

def graph(beta, gamma, mu, k):
    '''
    Creates our plot and represents the evolution
    of the disease through time.
    
    '''
    # Calling the euler function with following set agruments 
    sol = euler(derivatives, beta, gamma, mu, t0=0, tf=100, h=1, N0=100, I0=1)
    
    if beta == 0.8:
         contagion = " and very contagious"
    if beta == 0.1:
         contagion = " and not very contagious"
    if beta != 0.8 and beta != 0.1:
        contagion = ""
    if mu == 0:
         mortalite = "non-lethal"
    if mu == 0.02:
         mortalite= "lethal"    
    if mu == 0.1:
         mortalite = "very lethal"
    if mu == 0.5:
         mortalite = "fatal"
    if mu == 0.01:
         mortalite = "not very lethal"
    if gamma == 0.5:
         recuperation = "with fast recovery"
    if gamma != 0.5:
        recuperation = ""
        
    plt.figure(k) # In order to have different plots for the different conditions : k = 0, 1, 2, etc.  
    plt.plot(sol[4],sol[0], color="blue", label="S = healthy population")
    plt.plot(sol[4],sol[1], color= "red", label="I = infected population")
    plt.plot(sol[4],sol[2], color="green", label="R = recovered population")
    plt.plot(sol[4],sol[3], color="black", label="N = total population")
    plt.legend()
    plt.xlabel("Time") 
    plt.ylabel("Population")
    plt.title("beta = {}, gamma = {}, mu = {}".format(beta, gamma, mu), fontsize=10) #subtitle
    plt.suptitle("SIR model for a {}{} disease {}".format(mortalite, contagion ,recuperation), fontsize= 15) #title
    plt.grid()
    plt.show()



# First point             
graph(beta=0.3, gamma=0.1, mu=0, k=0)

# Second point     
graph(beta=0.3, gamma=0.1, mu=0.02, k=1)

# First type of disease    
graph(beta=0.8, gamma=0.1, mu=0.1, k=2)

# Second type      
graph(beta=0.8, gamma=0.1, mu=0.5, k=3)

# Third type      
graph(beta=0.8, gamma=0.1, mu=0.01, k=4)

# Fourth type      
graph(beta=0.1, gamma=0.1, mu=0.01, k=5)

# Fifth type      
graph(beta=0.1, gamma=0.5, mu=0.01, k=6)
