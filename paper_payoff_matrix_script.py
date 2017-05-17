#help: Integrates ODE systems for the default model and saves the output in a file which can be processed in MATLAB using "payoff_maker.m" to generate the payoff matrix and find its minimax point.

import numpy as np
from scipy.integrate import odeint

#Parameter set
a=20
g=1
b=100
d=1

# Time
t=np.linspace(0,50,51000)

#ODE
# One intermediate
def f(x,t):
      B0=x[0]
      B1=x[1]
      B2=x[2]
      L1=x[3]
      L2=x[4]
      P1=x[5]
      P2=x[6]
            
      Bt=B0+B1+B2+L1+L2

      dB0 = g*B0*(1-Bt) - a*B0*(P1+P2)
      dB1 = a*B0*P1 - d*B1
      dB2 = a*B0*P2 - d*B2
      dL1 = g*L1*(1-Bt) + f1*d*B1
      dL2 = g*L2*(1-Bt) + f2*d*B2
      dP1 = b*(1-f1)*d*B1 - a*P1*Bt
      dP2 = b*(1-f2)*d*B2 - a*P2*Bt
      return [dB0,dB1,dB2,dL1,dL2,dP1,dP2]
  
# Initial conditions
B0_0=0.001
B1_0=0
B2_0=0
L1_0=0
L2_0=L1_0
P0=1e-7

P1_0=P0
P2_0=P0
y0=[B0_0,B1_0,B2_0,L1_0,L2_0,P1_0,P2_0]

outf=open("pydel","w")

#Sweep over all calues of f1 and f2 from 0.01-1.00
f1=0.01
while f1<1:
    f2=0.01
    while f2<1:
        soln=odeint(f,y0,t)
        B0=soln[:,0]
        B1=soln[:,1]
        B2=soln[:,2]
        L1=soln[:,3]
        L2=soln[:,4]
        PF=(L1[-1]-L2[-1])/(L1[-1]+L2[-1])
        outf.write(str(f1*100)+"\t"+str(f2*100)+"\t"+str(PF)+"\n")
        f2=f2+0.01000
    f1=f1+0.01000
outf.close()
