#include <stdio.h>
#include <math.h>
#include <stdlib.h>
int nr,nth;
int N=120;
float get_bessel_zero(int m, int n)
{ 
  float zeroes[30]={2.4048,3.8317,5.1356,6.3802,7.5883,8.7715,5.5201,7.0156,8.4172,9.7610,11.0647,12.3386,8.6537,10.1735,11.6198,13.0152,14.3725,15.7002,11.7915,13.3237,14.7960,16.2235,17.6160,18.9801,14.9309,16.4706,17.9598,19.4094,20.8269,22.2178};
  int index=(n-1)*6+m;
  if (m>5 || n>5)  {   
    fprintf(stderr,"Higher Bessel function order requested than in my hacked-together table of the zeroes of the Bessel functions; you'll have to fix that. Exiting...\n");
    exit(-1);
  }  return zeroes[index];
}
float init_bessel(float x, float y,int m, int n){  float r=hypot(x,y);
  if (r>1) return 0;
  float zero=get_bessel_zero(m,n);
  float phi=atan2(x,y);
  return jn(m,zero*r)*0.3*cos(m*phi);
  printf("!Initializing with Bessel function: angular number %d, radial number %d, zero %f\n",m,n,zero);
}
int main(int argc, char **argv)
{
  if (argc < 7) 
  {
    printf("!Usage: <this> [b <angular mode number> <radial mode number>] or [t <center> <radius> <amplitude>] N dt frameskip\n");
    exit(0);
  }
  char mode=argv[1][0];
  int m,n,frameskip;
  double center,rad,str,dt,damping=0;
  if (mode == 'b')
  {
    m=atoi(argv[2]);
    n=atoi(argv[3]);
    N=atoi(argv[4]);
    dt=atof(argv[5]);
    frameskip=atoi(argv[6]);
  }
  if (mode == 't')  
  {
    center=atof(argv[2]);
    rad=atof(argv[3]);
    str=atof(argv[4]);
    N=atoi(argv[5]);
    dt=atof(argv[6]);
    frameskip=atoi(argv[7]);
  }
  printf("!Init: center %e radius %e amp %e\n",center,rad,str);

  int frame=0;
  float z[(N+1)*(N+1)];
  float v[(N+1)*(N+1)];
  float dx=2./(N+1);
  for (int i=0;
      i<=N;
      i++)    for (int j=0;j<=N;j++) 
  {
    float x=(float)((i-N/2)*2)/N;
    float y=(float)((j-N/2)*2)/N;
    v[i*(N+1)+j] = 0;
    if (mode == 'b') z[i*(N+1)+j] = init_bessel(x,y,m,n);
    if (mode == 't')
    {      
      z[i*(N+1)+j]=str*exp(-pow(hypot(x-center,y),2)/(rad*rad));
      z[i*(N+1)+j]+=str*exp(-pow(hypot(-x-center,y),2)/(rad*rad));
    }    
  }  
  double t;
  for (t;1;) 
  {
    for (int i=0;i<=N;i++) 
      for (int j=0;j<=N;j++)
	z[i*(N+1)+j] += v[i*(N+1)+j] * dt/2;
    for (int i=1;i<N;i++) 
      for (int j=1;j<N;j++)    
      {
	float x=(float)((i-N/2)*2)/N;
	float y=(float)((j-N/2)*2)/N;
	if (x*x+y*y>1) continue;
	float del;
	del = z[(i+1)*(N+1)+j] + z[(i-1)*(N+1)+j] + z[i*(N+1)+j+1] + z[i*(N+1)+j-1] - 4*z[i*(N+1)+j];
	v[i*(N+1)+j] += del / (dx*dx) * dt;
        v[i*(N+1)+j] -= v[i*(N+1)+j] * damping * dt;
      } 
    for (int i=0;i<=N;i++)
      for (int j=0;j<=N;j++)
	z[i*(N+1)+j] += v[i*(N+1)+j] * dt/2;
    frame++;
    if (frame%frameskip==0)
    {   
      for (int i=0;i<=N;i++)
	for (int j=0;j<=N;j++)
	{
	  float x=(float)((i-N/2)*2)/N;
	  float y=(float)((j-N/2)*2)/N;
	  if ((x*x+y*y)>1+dx*6) continue;
	  printf("C %.2e %.2e %.2e\n",0.5+z[(i)*(N+1)+j]*5,0.5,0.5-z[(i)*(N+1)+j]*5);
	  if ((x*x+y*y)<1+dx*1.5 && (x+dx)*(x+dx)+(y*y) < 1+dx*1.5 && i<N) 
	    if ((x*x+y*y)<1+dx*1.5 && (y+dx)*(y+dx)+(x*x) < 1+dx*1.5 && i<N) 
	 //     printf("c3 %.3e %.3e %.3e 0.01\n",x,y,z[(i)*(N+1)+j]);
	  	printf("q3 %.3e %.3e %.3e %.3e %.3e %.3e %.3e %.3e %.3e %.3e %.3e %.3e\n",x,y,      z[(i)*(N+1)+j]
	                                                  ,x,y+dx,   z[(i)*(N+1)+j+1]
	                                                   ,x+dx,y+dx,z[(i+1)*(N+1)+j+1]
	                                                   ,x+dx,y,   z[(i+1)*(N+1)+j]);
	}
      //      printf("c3 1 1 1 0.1\n");
      printf("C 1 1 1\n");
      printf("F\n");
    } 
  }
}
