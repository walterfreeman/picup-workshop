#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

// #define ONE_PERIOD_ONLY // uncomment this to do only one period, then print stats and exit. 
void get_forces(double *x, double *y, double *Fx, double *Fy, int N, double k, double r0) // pointers here are array parameters
{
  int i;
  double r;
  Fx[0]=Fy[0]=Fx[N]=Fy[N]=0; // have to set the end forces properly to avoid possible uninitialized memory shenanigans
  double U=0,E=0,T=0;
 
  for (i=1;i<N;i++)
  {
    // left force
    r = hypot(x[i]-x[i-1],y[i]-y[i-1]);
    Fx[i] = -(x[i]-x[i-1]) * k * (r-r0)/r;
    Fy[i] = -(y[i]-y[i-1]) * k * (r-r0)/r;

    // right force
    r = hypot(x[i]-x[i+1],y[i]-y[i+1]);
    Fx[i] += -(x[i]-x[i+1]) * k * (r-r0)/r;
    Fy[i] += -(y[i]-y[i+1]) * k * (r-r0)/r;
  }
}

 void evolve_leapfrog(double *x, double *y, double *vx, double *vy, int N, double k, double m, double r0, double dt) 
{
  int i;


  double Fx[N+1],Fy[N+1]; // this could be made faster by mallocing this once, but we leave it this way for students
                          // to avoid having to deal with malloc(). In any case memory allocation is faster than a bunch
                          // of square root calls in hypot().
  for (i=1;i<N;i++)
  {
    x[i] += vx[i]*dt/2;
    y[i] += vy[i]*dt/2;
  }
  get_forces(x,y,Fx,Fy,N,k,r0);
  for (i=1;i<N;i++)
  {
    vx[i] += Fx[i]/m*dt;
    vy[i] += Fy[i]/m*dt;
  }
  for (i=1;i<N;i++)
  {
    x[i] += vx[i]*dt/2;
    y[i] += vy[i]*dt/2;
  }
}


void evolve_euler(double *x, double *y, double *vx, double *vy, int N, double k, double m, double r0, double dt) 
{
  int i;


  double Fx[N+1],Fy[N+1]; // this could be made faster by mallocing this once, but we leave it this way for students
                          // to avoid having to deal with malloc(). In any case memory allocation is faster than a bunch
                          // of square root calls in hypot().
  get_forces(x,y,Fx,Fy,N,k,r0);
  for (i=1;i<N;i++)
  {
    x[i] += vx[i]*dt;
    y[i] += vy[i]*dt;
    vx[i] += Fx[i]/m*dt;
    vy[i] += Fy[i]/m*dt;
  }
}

// this function is around to go from Euler-Cromer to leapfrog, if we want second-order precision
void evolve_velocity_half(double *x, double *y, double *vx, double *vy, int N, double k, double m, double r0, double dt)
{
  int i;

  double Fx[N+1],Fy[N+1]; // this could be made faster by mallocing this once, but we leave it this way for students
  
  get_forces(x,y,Fx,Fy,N,k,r0);
  
  for (i=1;i<N;i++)
  {
    vx[i] += Fx[i]/m*dt/2;
    vy[i] += Fy[i]/m*dt/2;
  }
}

// Students might not be familiar with pass-by-reference as a trick for returning multiple values yet. 
// Ideally they should be coding this anyway, and there are a number of workarounds, in particular 
// just not using a function for this.
void get_energy(double *x, double *y, double *vx, double *vy, int N, double k, double m, double r0, double *E, double *T, double *U)
{
  *E=*T=*U=0;
  int i;
  double r;
  for (i=0;i<N;i++)
  {
    *T+=0.5*m*(vx[i]*vx[i] + vy[i]*vy[i]);
    r = hypot(x[i]-x[i+1],y[i]-y[i+1]);
    *U+=0.5*k*(r-r0)*(r-r0);
  }
  *E=*T+*U;
}

// does what it says on the tin
void evolve_euler_cromer(double *x, double *y, double *vx, double *vy, int N, double k, double m, double r0, double dt)
{
  int i;
  double Fx[N+1],Fy[N+1]; // this could be made faster by mallocing this once, but we leave it this way for students
  for (i=1;i<N;i++)
  {
    x[i] += vx[i]*dt;
    y[i] += vy[i]*dt;
  }
  
  get_forces(x,y,Fx,Fy,N,k,r0);
  
  for (i=1;i<N;i++)
  {
    vx[i] += Fx[i]/m*dt;
    vy[i] += Fy[i]/m*dt;
  }
}

// does what it says on the tin
void evolve_rk2(double *x, double *y, double *vx, double *vy, int N, double k, double m, double r0, double dt)
{
  int i;

  double Fx[N+1],Fy[N+1]; // this could be made faster by mallocing this once, but we leave it this way for students
  double xh[N+1],yh[N+1],vxh[N+1],vyh[N+1];
  vxh[0]=vyh[0]=vxh[N]=vyh[N]=0;
  get_forces(x,y,Fx,Fy,N,k,r0);

  for (i=0;i<=N;i++)
  {
    xh[i] = x[i] + vx[i]*dt/2;
    yh[i] = y[i] + vy[i]*dt/2;
    vxh[i] = vx[i] + Fx[i]/m*dt/2;
    vyh[i] = vy[i] + Fy[i]/m*dt/2;
  }
  
  get_forces(xh,yh,Fx,Fy,N,k,r0);
  
  for (i=0;i<=N;i++) // need two for loops -- can't interleave halfstep/fullstep updates (students have trouble with this sometimes!)
  {
    x[i] = x[i] + vx[i]*dt;
    y[i] = y[i] + vy[i]*dt;
    vx[i] = vx[i] + Fx[i]/m*dt;
    vy[i] = vy[i] + Fy[i]/m*dt;
  }
  
}

// function to encapsulate determining whether we need to shovel another frame to the animator. delay is the delay in msec.
// I usually aim to hide library calls, like clock() and CLOCKS_PER_SEC, from students in the beginning, since it's not really
// relevant to their development of computational skills, which are what I really care about.
int istime(int delay)
{
  static int nextdraw=0;
  if (clock() > nextdraw)
  {
    nextdraw = clock() + delay * CLOCKS_PER_SEC/1000.0;
    return 1;
  }
  return 0;
}

int main(int argc, char **argv)
{
  int i,N=80; //  number of links, not number of nodes!! Careful for the lurking fencepost errors
  int modenumber=3; // put in some defaults just in case 
  double t, dt=2e-6, amplitude=0.1;
  double stiffness=10, density=1, length=1; // unstretched properties of original string
  double k, m, r0; // properties of single string
  double tension=1,Ls;
  int frame=0, frameskip;

  if (argc < 9) // if they've not given me the parameters I need, don't just segfault -- tell the user what to do, then let them try again
  {
    printf("!Usage: <this> <N> <modenumber> <dt> <stiffness> <density> <length> <tension> <color>\n"); 
    exit(0);
  }
  N=atoi(argv[1]);
  modenumber=atoi(argv[2]);
  dt=atof(argv[3]);
  stiffness=atof(argv[4]);
  density=atof(argv[5]);
  length=atof(argv[6]);
  tension=atof(argv[7]);
  int col=atoi(argv[8]);
  //printf("!#x l y l lx \"Amplitude\" ly \"Deviation from SAA period\" r 0.05 cm %d\n",col);
  
  
  double x[N+1], y[N+1], vx[N+1], vy[N+1], E, T, U;

  // compute microscopic properties from macroscopic ones 

  r0=length/N;
  m=density*length/N;
  k=stiffness*N/length;

  // figure out stretched length

  Ls=length + tension * length / stiffness;

  // make predictions based on what our freshman mechanics class taught us

  double density_stretched = density * length / Ls;
  double wavespeed = sqrt(tension/density_stretched);
  double period_predict = 2 * Ls / wavespeed / modenumber;
  double vym_last=0;

  int monitor_node = N/modenumber/2; // this is the node that we'll be watching to see when a period has elapsed.
  int nperiods=0;  

  for (amplitude=1e-5; amplitude<0.5; amplitude*=1.5)
  {
  for (i=0;i<=N;i++) // remember, we have N+1 of these
  {
    x[i] = Ls*i/N;
    y[i] = amplitude*sin(modenumber * M_PI * x[i] / Ls);
    vx[i]=0;
    vy[i]=0;
  }

  // now, loop over time forever...
  for (t=0;1;t+=dt)
  {
    vym_last=vy[monitor_node];
    evolve_leapfrog(x,y,vx,vy,N,k,m,r0,dt);
    // "if we were going up, but now we're going down, then a period is complete"
    // this is crude and will fail if it has enough "wobble", but it's sufficient for this project. 
    if (vym_last > 0 && vy[monitor_node] < 0)
    {
      t+=dt;
      // now do extrapolation since we overshot
      t += vy[monitor_node] / (vym_last - vy[monitor_node]) * dt;
      nperiods++;
//      printf("!%d\t%d\t%.2e\t%.2e\t%.4e\t%.1e\t%.1e\t%.4f\t%.4f\t%.4e\t%.4e\n",N,
//modenumber,stiffness,length,amplitude,density,tension,t,period_predict*nperiods,1-t/period_predict/nperiods,nperiods/t);
      printf("!%e %e\n",amplitude,fabs(1-t/period_predict));
      break;
    }
    frame++; 
  //  if (frame % frameskip == 0) // wait 15ms between frames sent to anim; this gives us about 60fps.
    if (istime(15))
    {
      printf("C 1 0 0\nc %f %f 0.02\nC 1 1 1\n",x[monitor_node],y[monitor_node]); // draw a big red blob around the node we're watching
      for (i=0;i<=N;i++)
      {
        printf("C %e 0.5 %e\n",0.5+y[i]/amplitude,0.5-y[i]/amplitude); // use red/blue shading; this will make vibrations visible even if amp<<1
	printf("c %f %f %f\n",x[i],y[i],length/N/3); // draw circles, with radius scaled to separation
	if (i<N) printf("l %f %f %f %f\n",x[i],y[i],x[i+1],y[i+1]); // the if call ensures we don't drive off the array
      }
      printf("T -0.5 -0.7\ntime = %.4e amp=%.2e\n",t,amplitude); 
      get_energy(x,y,vx,vy,N,k,m,r0,&E,&T,&U);
      printf("T -0.5 -0.63\nenergy = %e + %e = %e\n",T,U,E);
      printf("F\n"); // flush frame
    }
  }
  }
  printf("Q\n");
} 
