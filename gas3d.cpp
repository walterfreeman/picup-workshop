#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include "vector.h"

int ilist_size=1000000; // how much memory to allocate for the list of possible interactions
int N=2000; // default number of particles (probably changed later)
double L=1; // box width
double r0=0.015; // interaction radius
int NG; // number of zones in each direction. Set later
double pad=2; // At what multiple of r0 do we assume the force is small and neglect the interactions?
double timer[10]; // Array shared by starttimer() and stoptimer() for timing things



// the following is a new thing y'all haven't learned yet. Ask me about in person; this is used to make
// "linked lists", which are a way of keeping track of which particles are where.
struct list_el {
  int val;
  struct list_el *next;
};
typedef struct list_el item;


double myclock(void)
{
  return clock()/1e6;
}

// Just a dopey timing function that tells you if n milliseconds have passed since the last time you called it
// "bool" = "boolean", true or false (1 or 0, basically)
bool istime2(int msec)
{
  static double lasttime=0;
  if (myclock() - lasttime > msec / 1000.0)
  {
    lasttime=myclock();
    return true;
  }
  else return false;
}

bool istime(int msec)
{
  static double lasttime=0;
  if (myclock() - lasttime > msec / 1000.0)
  {
    lasttime=myclock();
    return true;
  }
  else return false;
}

// given a particle's position, which zone is it in?
ivec getzone (vector v)
{
  ivec result;
  result.x=(int)((v.x+L/2)*NG/L);
  result.y=(int)((v.y+L/2)*NG/L);
  result.z=(int)((v.z+L/2)*NG/L);
  if (result.x>=NG) result.x=NG-1;
  if (result.x<0) result.x=0;
  if (result.y>=NG) result.y=NG-1;
  if (result.y<0) result.y=0;
  if (result.z>=NG) result.z=NG-1;
  if (result.z<0) result.z=0;
  return result;
}

// add/remove an item from linked lists -- will explain this in person
void add(item **list, int v)
{
  //   printf("!Adding item %d\n",v);
  item *curr;
  curr=(item *)malloc(sizeof(item));
  curr->val = v;
  curr->next = *list;
  *list=curr;
}

void del (item **list, int v)
{
  item *curr;
  curr=*list;
  if (curr == NULL) {return;}
  if (curr->val == v) {*list=curr->next; free(curr);}
  else
    while (curr -> next != NULL)
    {
      if (curr->next->val == v)
      {
	free(curr->next); curr->next=curr->next->next;
	break;
      }
      curr=curr->next;
    }
}



// cycle through all the particles and put them in the right zones
int check_zonelist(item **zonelist,ivec zone[], vector pos[])
{
  ivec cz; // the correct zone for the particle to be in

  for (int i=0;i<N;i++)
  {
    cz = getzone(pos[i]);
    if (cz.x != zone[i].x || cz.y != zone[i].y || cz.z != zone[i].z) // if we're not in the right zone, fix it
    {
      if (zone[i].x >= 0 && zone[i].x < NG && zone[i].y >= 0 && zone[i].y < NG && zone[i].z >= 0 && zone[i].x < NG) del(&zonelist[zone[i].x*NG*NG + zone[i].y * NG + zone[i].z],i);
      add(&zonelist[cz.x*NG*NG + cz.y*NG + cz.z],i); 
    }
    zone[i]=cz;
  }


}

// a function for V(r), used only to compute total energy
double V(double r)
{
  if (r > r0*pad) return V(r0*pad);
  r/=r0;
  static double r6;
  r6=r*r*r*r*r*r;
  return (4*r0*(1/(r6*r6)-1/r6));
}

double V2(double r)
{
  if (r*r > r0*r0*pad*pad) return V2(r0*r0*pad*pad);
  r/=(r0*r0);
  static double r6;
  r6=r*r*r;
  return (4*r0*(1/(r6*r6)-1/r6));
}


// derpy timing functions: this one starts a timer with index t...
void starttimer(int t)
{
  timer[t]=myclock();
}

// and this one stops it, returning how long it's been in microseconds since starttimer was called
double stoptimer(int t)
{
  return myclock()-timer[t];
}

// NOTE: Actually returns F/r to avoid having to do a square root in velocity_update
//       when dealing with the unit vectors. This is why it's r^14-r^8.
//       Square roots are 'spensive.
double force(double r2)
{
  static double r8,r14,r6;
  r6=r2*r2*r2; // doing this is substantially faster than using pow()
  r8=r6*r2;
  r14=r6*r8;
  return -24.*(2./r14 - 1./r8)  / r0; // the factor of 24 is in some canonical form of this. gods know why.
}

double drnd48(void) // just a way to get a centered random number
{
  return drand48()-0.5;
}

void draw_box(double L)
{
  // baaaaarf.
  printf("l3 -%e -%e -%e -%e -%e  %e\n",L,L,L,L,L,L);
  printf("l3 -%e -%e  %e -%e  %e  %e\n",L,L,L,L,L,L);
  printf("l3 -%e  %e  %e -%e  %e -%e\n",L,L,L,L,L,L);
  printf("l3 -%e  %e -%e -%e -%e -%e\n",L,L,L,L,L,L);
  printf("l3  %e -%e -%e  %e -%e  %e\n",L,L,L,L,L,L);
  printf("l3  %e -%e  %e  %e  %e  %e\n",L,L,L,L,L,L);
  printf("l3  %e  %e  %e  %e  %e -%e\n",L,L,L,L,L,L);
  printf("l3  %e  %e -%e  %e -%e -%e\n",L,L,L,L,L,L);
  printf("l3  %e -%e -%e -%e -%e -%e\n",L,L,L,L,L,L);
  printf("l3  %e  %e -%e -%e  %e -%e\n",L,L,L,L,L,L);
  printf("l3  %e -%e  %e -%e -%e  %e\n",L,L,L,L,L,L);
  printf("l3  %e  %e  %e -%e  %e  %e\n",L,L,L,L,L,L);
}



// This is the function I said I'd give you a copy of. What this does:
// * assumes that zonelist is an array of linked lists, one for each zone, that holds the particles in each zone
// * ... and that zone is an array of size N that holds which zone each particle is in
// * For each particle i in the simulation:
//   * For its own zone and all adjacent zones:
//     * For each particle j in those zones (iterated over with the while loop -- I'll explain this in person)
//       * If the distance from i to j is more than r0*pad, do nothing
//       * ... otherwise, add that pair of particles to ilist1 and ilist2, and their squared distance to r2list
// When we're done, ilist1 and ilist2 contain paired lists of all the interactions worth thinking about
// ... and r2list contains how far away those particles are, so you don't have to work it out again
// This returns the number of interactions that we found, so we know how far to take the for loops later

int build_interaction_list(vector pos[],item *zonelist[],ivec zone[],int ilist1[],int ilist2[], double r2list[])
{
  // need to make some internally-used interaction list arrays for the threads


//  starttimer(9);

  int num_found=0; // keep track of how many interactions we've found already 
  {  
    ivec myzone; //ivec = vector of integers (only used here)
    item *curr; 
    int i, j, k, l, m;

    for (i=0; i<N; i++)
    {
      double r2;

      myzone=getzone(pos[i]);
      for (k=myzone.x-1; k<=myzone.x+1; k++)
	for (l=myzone.y-1; l<=myzone.y+1; l++)
	  for (m=myzone.z-1; m<=myzone.z+1; m++)
	  {
	    // if we've driven off the edge of the simulation region, ignore and move to
	    // the next zone
	    if (k<0 || k>=NG || l<0 || l>=NG || m<0 || m>=NG) 
	      continue;
	    curr=zonelist[k*NG*NG + l*NG + m]; // this awkwardness avoids the issues with passing multidimensional arrays
	    // while there are particles left in this zone...
	    while (curr)
	    {
	      j=curr->val; 
	      // skip if pair misordered, same particle, or further apart than permitted
	      // this last condition may seem silly, but it allows us to consistently truncate
	      // the force so that the physics don't depend on the exact layout of the zone boundaries.
	      // (it also makes things faster.)
	      //r2=(pos[i]-pos[j])*(pos[i]-pos[j])/(r0*r0); // squared distance between, in units of r0
              r2 = ((pos[i].x - pos[j].x) * (pos[i].x - pos[j].x) +
                    (pos[i].y - pos[j].y) * (pos[i].y - pos[j].y) +
                    (pos[i].z - pos[j].z) * (pos[i].z - pos[j].z))/(r0*r0);
                 
	      if (i<=j || r2 > pad*pad) {curr=curr->next; continue;}
	      r2list[num_found] = r2; // it is very marginally faster to store a third array, containing the squared
	      // distances between particles, so we don't have to compute it again.
	      ilist1[num_found] = i;
	      ilist2[num_found] = j;
//	      	printf("!Thread %d found interaction between %d and %d, stored in %d (%d this thread so far)\n",tid,i,j,num_found[tid]+iliststride*tid,num_found[tid]);
		if (r2 < 0) printf("!ALARM: distance is %e\n",r2);
	      num_found++;
	      curr=curr->next;
	    }
	  }
    }
  }
  return num_found; 
}






void position_step(vector pos[], vector vel[], double dt)
{
  for (int i=0; i<N; i++)
    pos[i] += vel[i] * dt;
}









void velocity_step(vector *pos, vector *vel, double *m,
    int *ilist1, int *ilist2, double *r2list, int nint, double dt)
{
    double PE=0;
    vector sep;
    double r2;
for (int i=0; i<nint; i++) // traverse the interaction list in ilist1 and ilist2
  {
    vector F;
    int j,k;
    j=ilist1[i]; k=ilist2[i];
    sep=pos[j]-pos[k];
    //     r2 = sep*sep/r0*r0;
    //     if (r2 > pad*pad) continue;
    //     F=(force((sep*sep)/(r0*r0))*dt)*sep; // parentheses not absolutely necessary but may save a few flops
    F=(force((r2list[i]))*dt)*sep; // is it faster to look up the radius in the list, or compute it?
    // maybe for bigger problems that don't fit into cache the above is better
    vel[j] -= F/m[j];
    vel[k] += F/m[k];
  }
}


double kinetic(vector v[], double m[])
{
  double T=0;
  for (int i=0;i<N;i++)
    T+=0.5*m[i]*v[i]*v[i]; // note that v*v does a dot product
  return T;
}

double potential(vector pos[])
{
  double U=0;
  for (int i=0;i<N;i++)
    for (int j=i+1;j<N;j++)
    {
      U+=V(norm(pos[i]-pos[j]));
    }
  return U;
}









// dummy space




int main(int argc, char **argv)
{
  double time_intlist=0; 
  double time_velstep=0;
  double time_posstep=0;
  double time_check=0;
  double time_checklist=0;
  double time_anim=0;
  double lastchecksteps=0;
  double time_energy=0; // some variables for timing
  double dt=1e-3, vinit;
  double P=0, T=0, U, J=0; 
  int i=0,j=0,k=0,l=0;
  int ninttotal=0;
  int nint;
  int drawn;
  int *ilist1=NULL,*ilist2=NULL;
  double *r2list=NULL;
  double KE, PE,Tnow;
  ilist1=(int *)malloc(sizeof(int) * ilist_size);
  ilist2=(int *)malloc(sizeof(int) * ilist_size);
  r2list=(double *)malloc(sizeof(double) * ilist_size);
  double thermo_interval=10;
  double next_thermo=thermo_interval;
  double Taccum=0;
  double rate;
  if (argc < 5) // if they didn't give us command line parameters, die and print out the order
  {
    printf("!Usage: <this> <N> <dt> <vinit> <pad> <L>\n");
    printf("!Try: gas3d 1000 1e-3 0.25 1\n");
    exit(1);
  }
 
  // read command line parameters
  N=atoi(argv[1]);
  dt=atof(argv[2]);
  vinit=atof(argv[3]);
  L=atof(argv[4]);
  int drawms=N*80.0/1000.0;
  double smult=0;
  NG=1.0/(r0*pad);
  printf("!Read parameters. Drawing every %d ms\n",drawms);
  ivec zone[N];
  int lastframe=0;
  vector pos[N],v[N];
  double m[N];
  
  item *zonelist[NG*NG*NG];

  for (i=0;i<NG*NG*NG;i++) 
  {
    zonelist[i]=NULL;
  }

  // set up initial conditions
  
  double interval = r0 * pow(2,1./6.); // how far away to put the particles at the start
  for (i=0;i<N;i++)
  {
    j++;
    if (j>pow(N,1./3.)) {j=0;k++;}
    if (k>pow(N,1./3.)) {k=0;l++;}

    zone[i].x=-1; zone[i].y=-1; zone[i].z=-1; // set these to -1, we'll fix later
    pos[i].x=k*interval-interval * pow(N,1./3.)/2;
    pos[i].y=l*interval-interval * pow(N,1./3.)/2;
    pos[i].z=j*interval-interval * pow(N,1./3.)/2;
    v[i].x=drnd48()*vinit;
    v[i].y=drnd48()*vinit;
    v[i].z=drnd48()*vinit;
    m[i]=1;
  }

  check_zonelist(zonelist,zone,pos); // put particles in the right zones to start

  printf("!start main loop\n");

  int steps=0;
  for (double t=0; 1; t+=dt)
  {
       starttimer(9);
   if (istime2(4000))
     {
       KE=kinetic(v,m);
       PE=potential(pos);
     }
       time_energy+=stoptimer(9);

if (istime(drawms)) // anim time
   {
     starttimer(4);
     // it's actually super expensive to compute the total energy since we do it
     // pairwise to ensure absolute sanity, so do it only every 100 anim updates
          printf("C 0.5 1 0.5\n");
      // draw particles      
      for (i=0;i<N;i++)
      {
//        if (i==0) printf("C 1 1 1\n"); 
//        if (i==1) printf("C 0.5 1 0.5\n");
    //    if (i<1) printf("ct3 %d %.3e %.3e %.3e %.3e\n",i,pos[i].x,pos[i].y,pos[i].z,r0/4);
        printf("c3 %.4e %.4e %.4e %.4e\n",pos[i].x,pos[i].y,pos[i].z,r0/4);
      }
      printf("C 0.7 0.2 0.2\n");
      draw_box(L/2);
      printf("C 0.5 0.5 1\n");
      printf("T -0.9 0.85\nt=%.3f/%.2f, E=%.5e = %.5e + %.5e\n",t,next_thermo,KE+PE,KE,PE); stoptimer(0);
      printf("T -0.9 0.75\nPV = %.4e \t NkT = %.4e \t ratio = %f\n",P*L*L*L,N*T,P*L*L*L/(N*T));
      printf("T -0.9 0.65\nsteps/frame: %d\n",steps-lastframe);
      starttimer(0); lastframe=steps;
      
      printf("F\n");
//      fflush(stdout);
      time_anim += stoptimer(4);
   }
   
   // start force calculations

     // check zonelist
     // build list of interactions. pass in pointers
     // to interaction lists so it can realloc them
     // if need be
    position_step(pos,v,dt/2);

     starttimer(1);
     starttimer(2);
     check_zonelist(zonelist,zone,pos);
     time_checklist += stoptimer(2);
     nint=build_interaction_list(pos,zonelist,zone,ilist1,ilist2,r2list);
     time_intlist += stoptimer(1);
//     printf("!%d interactions\n",nint);
     

     starttimer(3);
     velocity_step(pos,v,m,ilist1,ilist2,r2list,nint,dt);
     position_step(pos,v,dt/2);

     double eps = 0.1786178958448;
     double lam = -0.2123418310626;
     double chi = -0.06626458266982;

    
     time_velstep += stoptimer(3);
     steps++;
     ninttotal+=nint;
   
     
 
     // all this code is just here to determine how much time different parts of the calculation are taking
     if (myclock() - 1.0 > time_check)
     {
       double curtime=myclock();
       printf("!%.2f force (%.2f ns/int), %.2f anim, %.2f ilist, %.2f checklist, %.2f energy (us/step): %.1f ints/atom, %.1f SPS\n",
            1e3*(float)time_velstep / (steps-lastchecksteps),
            1e6*(float)time_velstep/ (ninttotal),
            1e3*(float)time_anim / (steps-lastchecksteps),
            1e3*(float)time_intlist / (steps-lastchecksteps),
            1e3*(float)time_checklist / (steps-lastchecksteps),
            1e3*(float)time_energy / (steps-lastchecksteps),
            (float)ninttotal/N/(steps-lastchecksteps),
            (float)(steps-lastchecksteps)/(curtime-time_check));
       time_velstep=time_intlist=time_posstep=time_anim=time_energy=time_checklist=0;
       lastchecksteps=steps;
       time_check=myclock();
       ninttotal=0;
     }
     // do wall collisions and accumulate impulse (pressure calculation not in yet)
     for (i=0; i<N; i++)
     {   
       if ( (pos[i].y > L/2 && v[i].y > 0) || ((pos[i].y < -L/2 && v[i].y < 0) ) ) 
       {
         v[i].y=-v[i].y;
         J=J+2*fabs(v[i].y)*m[i];
       }
       if ( (pos[i].z > L/2 && v[i].z > 0) || ((pos[i].z < -L/2 && v[i].z < 0) ) ) 
       {
         v[i].z=-v[i].z;
         J=J+2*fabs(v[i].z)*m[i];
       }
       if ( (pos[i].x > L/2 && v[i].x > 0) || ((pos[i].x < -L/2 && v[i].x < 0) ) ) 
       {
         v[i].x=-v[i].x;
         J=J+2*fabs(v[i].x)*m[i];
       }
     }

     Tnow=0;
     for (i=0; i<N; i++)
     {
       Tnow += 0.5 * m[i] * v[i]*v[i] * dt;
     }
     Taccum+=Tnow;
     Tnow/=N;

     if (t > next_thermo)
     {
       P=J/(thermo_interval)/(6*L*L);
       T=Taccum * 2.0 / 3.0 / N / thermo_interval;
       J=Taccum=0;
       next_thermo = t + thermo_interval;
     }
 
  }
}
