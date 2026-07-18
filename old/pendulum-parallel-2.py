from numpy import *
from matplotlib.pyplot import *

dt=0.01
g=9.8
L=1
N=20
t=0

thetas=2*logspace(0,-2,N)
omegas=zeros(N)
periods=zeros(N)

for i in range(N):
  print ("trl",i,100);

while 1:
   # Aspel update
   t=t+dt
   thetas=thetas+omegas*dt
   omegas=omegas-g/L*(thetas)*dt

   for i in range(N):
       # fancy colors!
       print ("C",0.4*sin(i)+0.6,0.4*sin(i+2)+0.6,0.4*sin(i+4)+0.6)
       x=sin(thetas[i])
       y=-cos(thetas[i])
       z=i*0.1        
       print ("l3",0,0,z,x,y,z)
       print ("ct3",i,x,y,z,0.04)
   print ("T -0.1 0.9\nSimulation B\n");
   print ("F\n")
