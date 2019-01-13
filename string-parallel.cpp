#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <omp.h>
#include <sys/time.h>
// #define ONE_PERIOD_ONLY // uncomment this to do only one period, then print stats and exit. 
long int myclock(void)
{
  struct timeval tv;
  gettimeofday(&tv,NULL);
//  return (1000000 * tv.tv_sec + tv.tv_usec);
  return clock();
}

double hypot3(double x, double y, double z)
{
//  return hypot(x,y);
//  return sqrt(x*x+y*y);
  return sqrt(x*x+y*y+z*z);
}

void get_forces(double *x, double *y, double *z, double *Fx, double *Fy, double *Fz, int N, double k, double r0) // pointers here are array parameters
{
  int i;
  double r;
  Fx[0]=Fy[0]=Fx[N]=Fy[N]=0; // have to set the end forces properly to avoid possible uninitialized memory shenanigans
  double U=0,E=0,T=0;
 
  for (i=1;i<N;i++)
  {
    // left force
    r = hypot3(x[i]-x[i-1],y[i]-y[i-1],z[i]-z[i-1]);
    Fx[i] = -(x[i]-x[i-1]) * k * (r-r0)/r;
    Fy[i] = -(y[i]-y[i-1]) * k * (r-r0)/r;
    Fz[i] = -(z[i]-z[i-1]) * k * (r-r0)/r;

    // right force
    r = hypot3(x[i]-x[i+1],y[i]-y[i+1],z[i]-z[i+1]);
    Fx[i] += -(x[i]-x[i+1]) * k * (r-r0)/r;
    Fy[i] += -(y[i]-y[i+1]) * k * (r-r0)/r;
    Fz[i] += -(z[i]-z[i+1]) * k * (r-r0)/r;
  }
}

 void evolve_leapfrog(double *x, double *y, double *z, double *vx, double *vy, double *vz, int N, double k, double m, double r0, double dt) 
{
  int i;


  double Fx[N+1],Fy[N+1],Fz[N+1]; // this could be made faster by mallocing this once, but we leave it this way for students
                          // to avoid having to deal with malloc(). In any case memory allocation is faster than a bunch
                          // of square root calls in hypot().
  for (i=1;i<N;i++)
  {
    x[i] += vx[i]*dt/2;
    y[i] += vy[i]*dt/2;
    z[i] += vz[i]*dt/2;
  }
  get_forces(x,y,z,Fx,Fy,Fz,N,k,r0);
  for (i=1;i<N;i++)
  {
    vx[i] += Fx[i]/m*dt;
    vy[i] += Fy[i]/m*dt;
    vz[i] += Fz[i]/m*dt;
  }
  for (i=1;i<N;i++)
  {
    x[i] += vx[i]*dt/2;
    y[i] += vy[i]*dt/2;
    z[i] += vz[i]*dt/2;
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

// function to encapsulate determining whether we need to shovel another frame to the animator. delay is the delay in msec.
// I usually aim to hide library calls, like myclock() and CLOCKS_PER_SEC, from students in the beginning, since it's not really
// relevant to their development of computational skills, which are what I really care about.
int istime(int delay)
{
  static int nextdraw=0;
  if (myclock() > nextdraw)
  {
    nextdraw = myclock() + delay * CLOCKS_PER_SEC/1000.0;
    return 1;
  }
  return 0;
}

int main(int argc, char **argv)
{
  int i,N=80; //  number of links, not number of nodes!! Careful for the lurking fencepost errors
  int modenumber=3; // put in some defaults just in case 
  double t, dt=2e-6;
  double stiffness=10, density=1, length=1; // unstretched properties of original string
  double k, m, r0; // properties of single string
  double tension=1,Ls;
  int frame=0, frameskip;
  int lastdrawframe=0;
  long int lastdrawtime=0;
  if (argc < 12) // if they've not given me the parameters I need, don't just segfault -- tell the user what to do, then let them try again
  {
    printf("!Usage: <this> <N> <modenumber> <dt> <stiffness> <density> <length> <tension> <N_par> <amp_first> <amp_last> <cheat>\n"); 
    exit(0);
  }
  N=atoi(argv[1]);
  modenumber=atoi(argv[2]);
  dt=atof(argv[3]);
  stiffness=atof(argv[4]);
  density=atof(argv[5]);
  length=atof(argv[6]);
  tension=atof(argv[7]);
  int Np=atoi(argv[8]);
  double amp_first=atof(argv[9]);
  double amp_last=atof(argv[10]);
  int cheat=atoi(argv[11]);
  double amp_step=pow(amp_last/amp_first,1.0/(Np-1));
  double amplitude[Np];
  for (int p=0; p<Np; p++) amplitude[p]=amp_first * pow(amp_step,p);
  double pstep=2.0/Np;

  
  double x[(N+1)*Np], y[(N+1)*Np], z[(N+1)*Np], vx[(N+1)*Np], vy[(N+1)*Np], vz[(N+1)*Np], E[Np], T[Np], U[Np];

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
  double vym_last[Np]; for (int p=0;p<Np;p++) vym_last[p]=0;

  int nperiods[Np]; for (int p=0;p<Np;p++) nperiods[p]=0;

  for (int p=0;p<Np;p++)
  {
    for (i=0;i<=N;i++) // remember, we have N+1 of these
    {
      z[i+(N+1)*p] = 0; 
      x[i+(N+1)*p] = Ls*i/N - Ls/2;
      if (cheat==0) y[i+(N+1)*p] = amplitude[p]*sin(modenumber * M_PI * (x[i+(N+1)*p]+Ls/2) / Ls);
      if (cheat==1) y[i+(N+1)*p] = 1e-8*sin(modenumber * M_PI * (x[i+(N+1)*p]+Ls/2) / Ls);
      vx[i+(N+1)*p]=0;
      vy[i+(N+1)*p]=0;
      vz[i+(N+1)*p]=0;
    }
  }

  double velampnow[Np], velamplast[Np];
  for (int p=0; p<Np; p++) velampnow[p]=velamplast[p]=0;
  // now, loop over time forever...
  for (t=0;1;t+=dt)
  {
    for (int p=0; p<Np; p++)
    {
      
      velamplast[p]=velampnow[p];
      velampnow[p]=0;
      for (i=0; i<N; i++)
      {
	velampnow[p] += vy[i+(N+1)*p]*sin(modenumber * M_PI * (x[i+(N+1)*p]+Ls/2) / Ls); // more rigorous way to examine period and not get distracted by coupling to other modes
      }

    // "if we were going up, but now we're going down, then a period is complete"
    if (velamplast[p] > 0 && velampnow[p] < 0)
    {
      double tf=t + velampnow[p] / (velamplast[p] - velampnow[p]) * dt;
      // now do extrapolation since we overshot
      nperiods[p]++;
      if (nperiods[p]==1) printf("!%e %e\n",amplitude[p],fabs(1-tf/period_predict));
    }
    
    }
    //  if (frame % frameskip == 0) // wait 30ms between frames sent to anim; this gives us about 30fps.
    if (istime(30))
    {
      for (int p=0; p<Np; p++)
      {
	for (i=0;i<=N;i++)
	{
          double mult;
          if (cheat==0) mult=1;
          if (cheat==1) mult=1e8 * amplitude[p]; 
	  printf("C %a %a %a\n",0.5+y[i+(N+1)*p]*mult/amplitude[p],0.5,0.5-y[i+(N+1)*p]*mult/amplitude[p]); // use red/blue shading; this will make vibrations visible even if amp<<1
	  if (N < 100) printf("c3 %a %a %a %a\n",x[i+(N+1)*p],-(float)p*pstep+Np/2*pstep,y[i+(N+1)*p]*mult,length/N/2); // draw circles, with radius scaled to separation
	  if (i<N) printf("l3 %a %a %a %a %a %a\n",x[i+(N+1)*p],-(float)(p)*pstep+Np/2*pstep,y[i+(N+1)*p]*mult,x[i+1+(N+1)*p],-(float)p*pstep+Np/2*pstep,y[i+1+(N+1)*p]*mult); // the if call ensures we don't drive off the array
	}
      }
      printf("T -0.9 0.9\nSteps per frame: %d Steps per second: %.3f\n",frame-lastdrawframe,(float)(frame-lastdrawframe)*CLOCKS_PER_SEC/(myclock()-lastdrawtime));
      printf("F\n");fflush(stdout); // flush frame
      lastdrawframe=frame;
      lastdrawtime=myclock();
    }
    for (int p=0; p<Np; p++)
      evolve_leapfrog(x+(N+1)*p,y+(N+1)*p,z+(N+1)*p,vx+(N+1)*p,vy+(N+1)*p,vz+(N+1)*p,N,k,m,r0,dt);  // pointer hacking to point at the right parallel spring

    frame++;
  }
}
