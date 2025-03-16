# -*- coding: utf-8 -*-
"""
Created on Sat Dec 17 11:58:58 2022

@author: User
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation
from IPython.display import HTML

answer = int(input("Choose a condition: \n 0 = no wall \n 1 = wall without door \n 2 = wall with a door : "))

beta = float(input("Choose the value for beta between 0.1 and 0.8 \nThe greater the beta, the greater the contagion : "))
gamma = float(input("Choose the value for gamma between 0.1 and 0.5 \nThe greater the gamma, the faster the recovery : "))
mu = float(input("Choose the value for mu between 0 and 0.5 \nThe greater the mu, the greater the mortality : "))

# Matrix dimensions
x_max = 100
y_max = 100

# Initial conditions
I0 = float(10)    #initial number of infected people in an infected cell
N0 = float(100)   #initial number of people in a cell
S0 = float(N0-I0) #initial number of healthy people in an infected cell
R0 = float(0)     #initial number of healed people in a cell 

# 4 matrices at t=0
I_0 = np.zeros((x_max,y_max))
S_0 = np.full((x_max,y_max), N0)
R_0 = np.zeros((x_max,y_max))
N_0 = np.full((x_max,y_max), N0)

# Initial values for the first cell
I_0[0,0] = I0
S_0[0,0] = S0
    
def ddx(i_last, x, y):
    '''
    Defines the second partial derivative of I with respect to x.
    '''
    
    if x < 99 and x > 0:
        d2x = i_last[x+1,y] + i_last[x-1,y] - 2*i_last[x,y] # We take the elements of the latest version of the i matrix
    if x == 99:
        d2x = i_last[x-1,y] - 2*i_last[x,y]
    if x == 0:
        d2x = i_last[x+1,y] - 2*i_last[x,y]
        
    return d2x

def ddy(i_last, x, y):
    '''
    Defines the second partial derivative of I with respect to y.
    '''    
    
    if y < 99 and y > 0:
        d2y = i_last[x,y+1] + i_last[x,y-1] - 2*i_last[x,y]
    elif y == 99:
        d2y = i_last[x,y-1] - 2*i_last[x,y]
    elif y == 0:
        d2y = i_last[x,y+1] - 2*i_last[x,y]
    
    return d2y


def derivatives(S, I, R, N, beta, gamma, mu):
    '''
    Returns a list of 4 matrices that evaluate the derivatives at a given time t.
    '''
    devS = np.zeros((x_max,y_max))
    devI = np.zeros((x_max,y_max))
    devR = np.zeros((x_max,y_max))
    devN = np.zeros((x_max,y_max))
    
    for xx in range (0,100): 
        for yy in range(0,100):
            
            # Derivative terms
            
            beta_term = (beta*S[xx,yy]/N[xx,yy])*(I[xx,yy]+ddx(I,xx,yy)+ddy(I,xx,yy))
            gamma_term = gamma*I[xx,yy]
            mu_term = mu*I[xx,yy]
            
            # Setting security limits
            
            if beta_term > S[xx,yy]: 
                beta_term = S[xx,yy]
            
            # Differential equations 
            
            dS_dt = -beta_term
            dI_dt = beta_term - gamma_term - mu_term
            dR_dt = gamma_term
            dN_dt = -mu_term
            
            devS[xx,yy] = dS_dt
            devI[xx,yy] = dI_dt
            devR[xx,yy] = dR_dt
            devN[xx,yy] = dN_dt
            
    return [devS, devI, devR, devN]


def euler(derivatives, beta, gamma, mu, t0, tf, h):
    '''
    Returns a list, in which each element 
    is a list of matrices that vary over time.
    '''
    # Initial matrices
    s = S_0
    i = I_0
    r = R_0
    n = N_0
    
    # Setting the lists of the different types of matrices
    S1 = [s]
    I1 = [i]
    R1 = [r]
    N1 = [n]
    
    t = t0
    
    while t < tf:
        s_last = s
        i_last = i
        r_last = r
        n_last = n
        
        # Shows "time" passing as matrices are created
        if t == 100:
            print("It takes a bit time to animate all this...")
        elif t == 250:
            print("Half way through.")
        elif t == 400:
            print("Just a little more...")
        elif t == 450:
            print("Almost there !")
        elif t == 490:
            print("Are you ready ?")
        else:
            print("t =", t)
        
        s = np.add(s, h*derivatives(s_last,i_last,r_last,n_last, beta,gamma,mu)[0])
        i = np.add(i, h*derivatives(s_last,i_last,r_last,n_last, beta,gamma,mu)[1])
        r = np.add(r, h*derivatives(s_last,i_last,r_last,n_last, beta,gamma,mu)[2])
        n = np.add(n, h*derivatives(s_last,i_last,r_last,n_last, beta,gamma,mu)[3])
        
        # Defining a wall with or without a door
        
        if answer == 1:             # Having just a wall
            s[50,] = np.zeros(100)
            i[50,] = np.zeros(100)
            
        elif answer == 2:             # Having a wall with a door
            s[50,:69] = np.zeros(69)
            s[50,71:] = np.zeros(29)
            i[50,:69] = np.zeros(69)
            i[50,71:] = np.zeros(29)
                 
        S1.append(s)
        I1.append(i)
        R1.append(r) 
        N1.append(n)
        
        t += h

    return [S1,I1,R1,N1]

# Calling the euler function with step = 1 and timespan between 0 and 500.
sol = euler(derivatives, beta, gamma, mu, t0=0, tf=500, h=1)

solS = sol[0]
solI = sol[1]
solR = sol[2]
solN = sol[3]
solS[0][99,99]

# Setting the animation #
fig, ax = plt.subplots() # Creating a plot
mat=ax.matshow(solI[0]) # Starting animation with initial I matrix

def update(data): # Function that update the matrix
    mat.set_data(data)
    return mat

plt.colorbar(mat) # Color bar for the matrix
plt.title("beta = {}, gamma = {}, mu = {}".format(beta, gamma, mu))
# Starting the animation for I
anim=animation.FuncAnimation(fig, update, solI, interval=100, blit=False)
HTML(anim.to_html5_video())

FFwriter = animation.FFMpegWriter(fps=10)
anim.save(filename= 'animation-sans-mur.mp4', writer=FFwriter)
