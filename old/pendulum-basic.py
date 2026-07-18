from numpy import *

dt=0.01
g=9.8
L=1
theta=1
omega=0

while True: 
    # Euler-Cromer-Aspel update
    theta=theta+omega*dt
    omega=omega-g/L*sin(theta)*dt
    
    # Animation
    x=sin(theta)
    y=-cos(theta)
    print ("l3",0,0,0,x,y,0)
    print ("ct3",0,x,y,0,0.05)
    print ("F")
