from scipy.integrate import odeint
import matplotlib.pyplot as plt
import numpy as np

d1 = 0.639/5
d2 = 0.693/35

k1 = 2.5*d1
k2 = 1000*d2

Km = 1
R = 1
n = 1

state = []

mRNA0 = 0
Prot0 = 0
mRNA01 = 0
Prot01 = 0

state.append([mRNA0, Prot0, mRNA01, Prot01])

def constitutive_promoter(state, t, k1, k2, Km, R, n):
    
    mRNA = state[0]
    Prot = state[1]
    mRNA1 = state[2]
    Prot1 = state[3]
    
    dmRNA_dt = k1*Km**n/(Km**n +R**n) - d1*mRNA
    dProt_dt = k2*mRNA - d2*Prot
    
    dmRNA_dt1 = k1 - d1*mRNA1
    dProt_dt1 = k2*mRNA1 - d2*Prot1

    return dmRNA_dt, dProt_dt, dmRNA_dt1, dProt_dt1
    
    
t  = np.linspace(0, 500., 1000)

for i in state:

# Solve the ODEs
    solution = odeint(constitutive_promoter, i, t, args = (k1, k2, Km, R, n))
    X = solution[:, 0]
    Y = solution[:, 1]
    Z = solution[:, 2]
    W = solution[:, 3]
    
    plt.figure()
    plt.plot(t, Y, 'k', label='Prot')
    plt.plot(t, W, 'b', label='Prot2')
    plt.xlabel('Time')
    plt.ylabel('Concentration')
    plt.legend(loc=0)
    
    plt.figure()
    plt.plot(t, X, 'y', label='mRNA')
    plt.xlabel('Time')
    plt.ylabel('Concentration')
    plt.plot(t, Z, 'r', label='mRNA2')
    plt.legend(loc=0)
    
