#! /usr/env/python
#help: Runs a user chosen number of iterations for pitting two phages against each other. Both phages get to scan through a f(MOI), the winner keeps its f(MOI) while the loser gets to alter it slightly. We also store the payoff values and plot figures to get an idea of which is the best f(MOI) in a random search. 

#tags: moi search random boxplot winner 
import numpy as np
from scipy.integrate import odeint
import math
import random as rand


# Parameters and initial conditions
eta=20
beta=100
gamma=1
delta=1
bact0=1e-3
phage0=1e-7

F=np.linspace(0.00,1.00, 101)
F1=[]   
for f in F:
  F1.append(f)


# ODE system,3 intermediates
def fd(x,t):
      B0=x[0]
      B11=x[1]
      B12=x[2]
      B13=x[3]
      B21=x[4]
      B22=x[5]
      B23=x[6]
      L1=x[7]
      L2=x[8]
      P1=x[9]
      P2=x[10]
      
      Bt=sum(x[0:9])
      
      #f21,f22,f23 are called from outside
      dB0 = g*B0*(1-Bt) - a*B0*(P1+P2)
      dB11 = a*B0*P1 - a*B11*P1 - 3*d*B11
      dB12 = a*B11*P1 - a*B12*P1 - 3*d*B12
      dB13 = a*B12*P1 - 3*d*B13
      dB21 = a*B0*P2 - a*B21*P2 - 3*d*B21
      dB22 = a*B21*P2 - a*B22*P2 - 3*d*B22
      dB23 = a*B22*P2 - 3*d*B23
      dL1 = g*L1*(1-Bt) + 3*d*(f11*B11+f12*B12+f13*B13)
      dL2 = g*L2*(1-Bt) + 3*d*(f21*B21+f22*B22+f23*B23)
      dP1 = b*3*d*((1-f11)*B11+(1-f12)*B12+(1-f13)*B13) - a*P1*Bt
      dP2 = b*3*d*((1-f21)*B21+(1-f22)*B22+(1-f23)*B23) - a*P2*Bt
      return [dB0,dB11,dB12,dB13,dB21,dB22,dB23,dL1,dL2,dP1,dP2]
#Time
t=np.linspace(0,100,101000)
T=[]
for time in t:
  T.append(time)
  
user_set_iterations=int(raw_input("Enter integral number of game iterations:"))
#user_set_iterations=10000

moiwin=open("moi_winners_v_time","w") 
'''
File to store output. 
'''
a=eta 
b=beta
g=gamma
d=delta
B0_0=bact0
p0=phage0
  #Common initial conditions
L1_0=0
L2_0=L1_0
P1_0=p0
P2_0=P1_0

y0=[B0_0,0,0,0,0,0,0,L1_0,L2_0,P1_0,P2_0]
iteration=1
winner=[]
loser=[]
draw=[]
f1=0.1 #The best fixed strategy


#random start of f1i and f2i
f11=rand.choice(F1)
f12=rand.choice(F1)
f13=rand.choice(F1)
f21=f11
f22=f12
f23=f13
winners=[]
while iteration < user_set_iterations:
  
    soln=odeint(fd,y0,t)
    B0=soln[:,0]
    B11=soln[:,1]
    B12=soln[:,2]
    B13=soln[:,3]
    B21=soln[:,4]
    B22=soln[:,5]
    B23=soln[:,6]
    L1=soln[:,7]
    L2=soln[:,8]
    PF1 = (L1[-1]-L2[-1])/(L1[-1]+L2[-1])
    toggle=0 # This temporary variable toggles between 1 and 2, just to keep track of the previous winner so we don't append the winning strategy to that player's quota. 
    if L1[-1]>L2[-1]: # P1 wins
        moiwin.write(str(f11)+"\t"+str(f12)+"\t"+str(f13)+"\n")
        winners.append([f11,f12,f13])
        toggle=1
        f21=min(1.00,max(0.00,rand.choice(F1[max(0,F1.index(f11)-3):min(100,F1.index(f11)+3)])))
        f22=min(1.00,max(0.00,rand.choice(F1[max(0,F1.index(f12)-3):min(100,F1.index(f12)+3)])))
        f23=min(1.00,max(0.00,rand.choice(F1[max(0,F1.index(f13)-3):min(100,F1.index(f13)+3)])))
    elif L1[-1]<L2[-1]: # P2 wins
        moiwin.write(str(f21)+"\t"+str(f22)+"\t"+str(f23)+"\n")
        winners.append([f21,f22,f23])
        toggle=2
        f11=min(1.00,max(0.00,rand.choice(F1[max(0,F1.index(f21)-3):min(100,F1.index(f21)+3)])))
        f12=min(1.00,max(0.00,rand.choice(F1[max(0,F1.index(f22)-3):min(100,F1.index(f22)+3)])))
        f13=min(1.00,max(0.00,rand.choice(F1[max(0,F1.index(f23)-3):min(100,F1.index(f23)+3)])))
        
        
    elif abs(L1[-1]-L2[-1])<0.0000000001: #draw condition. 
        #print f11,f12,f13,f21,f22,f23
        winners.append([f21,f22,f23])
        moiwin.write(str(f21)+"\t"+str(f22)+"\t"+str(f23)+"\n") #record P2 victory and mutate it.
        f21=min(1.00,max(0.00,rand.choice(F1[max(0,F1.index(f11)-3):min(100,F1.index(f11)+3)])))
        f22=min(1.00,max(0.00,rand.choice(F1[max(0,F1.index(f12)-3):min(100,F1.index(f12)+3)])))
        f23=min(1.00,max(0.00,rand.choice(F1[max(0,F1.index(f13)-3):min(100,F1.index(f13)+3)])))

    iteration=iteration+1

moiwin.close()
