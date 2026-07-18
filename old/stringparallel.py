from numpy import *
from time import *
from sys import *

stiffness = 10
L = 1
density = 1
N = 50 
M = 20
tension = 10
dt=1e-5
step = 0
skip = 10
zstep = 1/M
lasttime = 0
framedelay = 60 # milliseconds
laststeps = 0

ktotal = stiffness / L
r0 = L / N
Lstr = L + tension/ktotal
k = ktotal * N
m = density * L / N

pos    = zeros((M,N+1,3))
vel    = zeros((M,N+1,3))
fleft  = zeros((M,N+1,3))
unitvecs = zeros((M,N,3))
ampl=logspace(0,-2,M)* 0.5


pos[:,:,0]= tile(linspace(0,Lstr,N+1),(M,1))
pos[:,:,1] = ampl[:,newaxis] * sin(pos[:,:,0]*pi/Lstr)

while 1:
    pos += vel * dt/2
    
    offsets = pos[:,:-1]-pos[:,1:]   # vector from i to i+1
    radii = linalg.norm(offsets,axis=2)   # distance from i to i+1
    unitvecs = offsets / radii[:,:,newaxis] # unit vector from i to i+1
    fleft[:,1:]=unitvecs*k*(radii[:,:,newaxis] - r0)

    vel[:,1:-1] += (fleft[:,1:-1] - fleft[:,2:])/m * dt

    pos += vel * dt/2
    
    if time()-lasttime > framedelay*0.001:
        temptime=time()
        for j in range(0,M):
            for i in range(0,N+1):
                print ("C %.2f %.2f %.2f" % (1+pos[j,i,1]/ampl[j],0.5,1-pos[j,i,1]/ampl[j]))
                print ("c3 %.3f %.3f %.3f %.3f" % (pos[j,i,0],pos[j,i,1],pos[j,i,2]+j*zstep,r0*.8))
                if (i < N):
                    print ("l3 %.3f %.3f %.3f %.3f %.3f %.3f" % (pos[j,i,0],pos[j,i,1],pos[j,i,2]+j*zstep,pos[j,i+1,0],pos[j,i+1,1],pos[j,i+1,2]+j*zstep))
        print ("T -0.5 0.9\nSPS %.2f" % ((step-laststeps) / (time() - lasttime)))
        print ("T -0.5 0.85\nSPF %d" % ((step-laststeps) ))
        print ("T -0.5 0.80\n%f ms to draw\n" % ( 1000*(  time() - temptime)  ) )
        print ("F")
        stdout.flush()
        lasttime=time()
        laststeps=step

    step += 1
