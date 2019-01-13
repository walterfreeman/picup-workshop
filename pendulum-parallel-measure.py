from numpy import *
from matplotlib.pyplot import *

dt=0.0001
g=9.8
L=1
N=100
t=0
steps=0
animinterval=100

thetas=3.1*logspace(0,-5,N)
omegas=zeros(N)
periods=zeros(N)

for i in range(N):
  print ("trl",i,100);

while 1:
   steps = steps + 1
   # leapfrog update
   t=t+dt
   omegas=omegas-g/L*sin(thetas)*dt
   thetas=thetas+omegas*dt

   # animate one every animinterval frames
   if (steps % animinterval == 0):
       for i in range(N):
           # fancy colors!
           print ("C",0.4*sin(i)+0.6,0.4*sin(i+2)+0.6,0.4*sin(i+4)+0.6)
           print ("l3",0,0,i*0.1,sin(thetas[i]),-cos(thetas[i]),i*0.1)
           print ("ct3",i,sin(thetas[i]),-cos(thetas[i]),i*0.1,0.04)
       print ("F\n")

   # check for periods
   for i in range(N):
       if (periods[i]==0 and omegas[i]>0):
#           periods[i]=2*(t + omegas[i]/(g/L*sin(thetas[i])))
           periods[i]=2*t
           print ("!",-thetas[i],abs(periods[i]-2*pi*sqrt(L/g)))
