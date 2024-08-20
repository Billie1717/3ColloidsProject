/******************************************************/
/******* FLUID MEMBRANE    DT                  ********/
/******* July 1 2003                           ********/
/******************************************************/
 
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define L 80  // it must be an even number
#define Lz L
// #define N 1922  // 1082, 3002 available 
// #define Ncoll 5
// #define D0 5

#define MCSTEPS 21844000
#define WRITE 2000 //5000
#define WRITE_CONF 2000 //5000
#define EPSILON 1e-8

int Ntri;
int N,Ncoll;
int counter;
double D0;
double phi1, phi2, phi3, theta1, theta2, theta3;
int seed;
int Tmax=10; // never allow more then Tmax-coordinated points. 
typedef struct{
  double pres;

  double sidex,sidey,Isidex,Isidey,side2x,side2y;
  double sidez,Isidez;
  
  double vol,Ivol,volC;
  int Lcelx,Lcely,Lcelz,Ncel;
  double celsidex,celsidey,Icelsidex,Icelsidey;
  double celsidez,Icelsidez;

  double etaC,etaR,etaR2;
  double sig1[9000+1],sig2,sig1_sq,sum_sig,um_sig_sq,q;
  double rv,rc,rb,rv2,rb2,Vgap;            
  double rb2_2D;
  double C0,C1,C2,CC;
  double eps_mc,eps_mv;  
  double kappa;
  double rad;
  double Svol;
  //double EnergyTotal;
  //double BindTotal;
  //double BendTotal;
  //double Work;
  //double EnCol;
  //double Work2;
  double EbendChange;
  double EbindChange;
  double EstretchChange;
  double springEnN;
  double springEnO;
  double springEnCh;
    
} SYSTEM;

typedef struct{
  double x,y,z;      /* Coordinates (x,y,z)         */
  double Vx,Vy,Vz;   /* Verlet Coordinates (x,y,z) for bead-bead interactions */
  double vx,vy,vz;   /* Verlet Coordinates (x,y,z) for colloid-bead interactions */
  int verl[2700],bond[2700],clus[2700];
  int list[9100+1];
  int nlist;
  int type;
  int Nverl,Nbond;
  int cellID;
  int before,after;
  int Ng,ng[15]; 
  int ngR[15];
  int ngL[15];
} PARTICLE;

typedef struct{
  int Nv,v[3+1];
  int Nt,t[3+1];
} TRIANGLE;

typedef struct 
{
  int ngb[30];
  int begin;  
} CELL;


typedef struct{
  double x,y;
} vec2D;

typedef struct{
  double x,y,z;
} vec3D;

typedef struct{
  double theta,phi;
} ANGOLO;

typedef struct{
 double sig;
 double dis_mc,rv,rv2,rc;
} COLL;


SYSTEM S;
PARTICLE *Elem;
TRIANGLE *Tr;
CELL *Cella;
COLL Coll;

char readname[100];

char Roname[100];
char Prob_name[100];

FILE *read,*wr1,*wr2,*prob;

void set_all(void);
void set_Elements(void);
void set_nb(void);
void set_triangles(void);
void set_Cells(void);
void set_neighbors(void);
void set_Verlet();
void switch_bond(int);
void b_coll_vlist(void);

int Interaction_1p(int);
int Interaction_ALL(void);
int Interaction_1coll(int);
double bending_1p(int,int);
double bending_edge(int,int);
double bind_1p(int);
double volume_change(int,int);
double volume_change_2tr(int,int,int,int,int,int);
double interaction_coll(int);
void MC_xyz(int);
void MC_xyz_coll(int);

void read_Elements(void);
void write(void);
void painter(void);
void painter_en(int,int,int,int,int);
void painter2(void);
void painter3(void);
void painter4(void);

vec3D normale(int,int,int);
void ppp(int,int,int,int);
double lattice,lattice_y;
double cutoff,cutoff2;
vec3D rotate(double , double ,double  ,double, double);
vec3D rotx(vec3D,double);
vec3D rotz(vec3D,double);
ANGOLO angle(double , double ,double);

double cm_phi,cm_theta;
vec3D Centre;

main(int argc, char* argv[]){
  int i,j,k,o,coll;
  int touch=0,contat;
  double en,PP,SS;
  ANGOLO cm;
  vec3D p;
 FILE *in0;
 counter = 0;
  lattice=1.35;
 S.EbendChange=0;
 in0=fopen("in.ves","r");
 fscanf(in0,"%lf %lf %d %d %lf %lf %lf %d\n",&S.kappa,&S.pres,&N,&Ncoll,&Coll.sig,&D0,&Coll.rc,&seed); 
 printf("Parameters %lf %lf %d %d %lf %lf %d\n",S.kappa,S.pres,N,Ncoll,Coll.sig,D0, seed); 
 Coll.dis_mc=0.1;
 
  Coll.rc=Coll.rc*(0.5+0.5*Coll.sig); 
   Coll.rv=1.6*Coll.rc;
   Coll.rv2=Coll.rv*Coll.rv;
 
  lattice_y=lattice*sqrt(3.)/2.;
  set_all();
   
   for (i=1;i<=MCSTEPS;i++){ //This is where the loop is
    //for (j=1;j<=Ncoll;j++) MC_xyz_coll(1+(int)(drand48()*(Ncoll+N)));   
    //for (j=1;j<=Ncoll;j++) MC_xyz_coll(N+(int)(drand48()*(Ncoll)));   
    for (j=1;j<=Ncoll;j++) MC_xyz_coll(N+1+(int)(drand48()*(Ncoll)));
    for (j=1;j<=N;j++) MC_xyz(1+(int)(drand48()*Ntri));   
    for (j=1;j<=N;j++) switch_bond(1+(int)(drand48()*Ntri));
       
        
    if (i%WRITE==0){ 
    S.etaC=4./3.*M_PI*(double)(N)/8./S.vol;
    printf("%d V=%f C=(%2.3f, %2.3f, %2.3f) \n",
	   i,S.Svol,Centre.x,Centre.y,Centre.z);  
    }
    if (i == 1 || i%WRITE_CONF==1 && i>10) {write();painter3();painter();painter_en(i,counter, N, Ncoll,Ntri);counter+=1;}
  }

  write();
  painter3();
}
/*---------------------------*/
void set_all(void){
  int i;

  srand48(seed);
  //srand48(time(NULL)*(2));
  S.rad=1.00;
  S.Svol=(4./3.)*M_PI*S.rad*S.rad*S.rad;

  cutoff=sqrt(3);
  cutoff2=cutoff*cutoff;

  S.sidex=(double)(L+2)*lattice;
  S.sidey=(double)(L+2)*lattice_y;
  
  S.side2x=S.sidex/2.;
  S.side2y=S.sidey/2.;
  
  S.sidez=S.sidex;

  S.vol=S.sidex*S.sidey*S.sidez;
  S.Isidex=1./S.sidex;
  S.Isidey=1./S.sidey;

  S.Isidez=1./S.sidez; 
  
  S.Ivol=1./S.vol;
  for (i=1;i<=N;i++){
     S.sig1[i]=1.0;
  }
 
  // ALWAYS SET rv=rb   
  S.rv=2.;  // verlet minimum 
  S.rb=2.;  // bond order minimum
  S.rc=1.78;  // cluster minimum 
 
  S.rb2_2D=1.5*1.5; // Bond Order Cut-OFF for 2D order parameter

  S.rv2=S.rv*S.rv;
  S.rb2=S.rb*S.rb;
  S.Vgap=.5*(S.rv-1);  // put the largest instead of 1.
  
  S.Lcelx=(int)(S.sidex/S.rb); 
  S.Lcely=(int)(S.sidey/S.rb);
  S.Lcelz=(int)(S.sidez/S.rb);
  S.Ncel=S.Lcelx*S.Lcely*S.Lcelz;

 
  S.celsidex=S.sidex/(int)(S.sidex/S.rb);
  S.celsidey=S.sidey/(int)(S.sidey/S.rb);
  S.celsidez=S.sidez/(int)(S.sidez/S.rb);
  
  S.Icelsidex=1./S.celsidex;
  S.Icelsidey=1./S.celsidey;
  S.Icelsidez=1./S.celsidez;
  
 
  if (S.Lcelx>L){
    printf("Lcel=%d> SIDE OF BOX=%d System too sparse increase rv and rb\n",
	   S.Lcelx,L);exit(-1);} 

  if (S.Lcely>L){
    printf("Lcel=%d> SIDE OF BOX=%d System too sparse increase rv and rb\n",
	   S.Lcely,L);exit(-1);} 

  if (S.Lcelz>Lz){
    printf("Lcelz=%d> SIDEZ OF BOX=%d System too sparse increase rv and rb\n",
	   S.Lcelz,Lz);exit(-1);} 
 
  if (S.celsidex<=S.rb) {printf("Lcel=%lf<Rb decrese Rb\n",S.celsidex);exit(-1);}
  if (S.celsidey<=S.rb) {printf("Lcel=%lf<Rb decrese Rb\n",S.celsidey);exit(-1);}
  if (S.celsidez<=S.rb) {printf("Lcelz=%lf<Rb decrese Rb\n",S.celsidez);
  exit(-1);}
 
 
  S.eps_mc=0.10;   
  S.eps_mv=0.0008;   
   
  Elem=(PARTICLE *) malloc((N+Ncoll+1)*sizeof(PARTICLE));
  Tr=(TRIANGLE *) malloc((8*N+1)*sizeof(TRIANGLE));  
  Cella=(CELL *) malloc((L*L*Lz+1)*sizeof(CELL));
  
  sprintf(readname,"conf_N%d_P%lf.dat",N,S.pres);
  read=fopen(readname,"r");
  if (!read)  set_Elements();
  else read_Elements();
    
  sprintf(Roname,"ro_N%d.dat",N);
  wr2=fopen(Roname,"w");
  fclose(wr2);

 
  set_Cells();
  set_Verlet();
  b_coll_vlist();
	
}
/*---------------------------*/
void set_Elements(void){
  int i,j,k,t=0,p;
  double ssside,disp,vvv;
  double dx,dy,dz,d2;
  double Ex,Ey,Ez;
  char to[100];
  FILE *in;
 
  //sprintf(to,"../icos_%d_%d_%2.1f.dat",N,Ncoll,Coll.sig); 
  sprintf(to,"icos_%d_%d_%d.dat",N,Ncoll,(int)(Coll.sig));
  printf("%s\n",to);
  in=fopen(to,"r");
  
  S.rad=1.0;
  S.Svol=(4./3.)*M_PI*S.rad*S.rad*S.rad;
  
  for (t=1;t<=N+Ncoll;t++){ 
    fscanf(in,"%d %lf %lf%lf\n",&Elem[t].type,&Elem[t].x,&Elem[t].y,&Elem[t].z);
    
    Elem[t].x*=S.rad; Elem[t].y*=S.rad; Elem[t].z*=S.rad;
 
    Elem[t].x+=.001*(.5-drand48());
    Elem[t].y+=.001*(.5-drand48());
    Elem[t].z+=.001*(.5-drand48());
  }  
    
  theta1 = atan((sqrt(Elem[N+Ncoll-2].x*Elem[N+Ncoll-2].x+Elem[N+Ncoll-2].y*Elem[N+Ncoll-2].y)/Elem[N+Ncoll-2].z));
  theta2 = atan((sqrt(Elem[N+Ncoll-1].x*Elem[N+Ncoll-1].x+Elem[N+Ncoll-1].y*Elem[N+Ncoll-1].y)/Elem[N+Ncoll-1].z));
  theta3 = atan((sqrt(Elem[N+Ncoll].x*Elem[N+Ncoll].x+Elem[N+Ncoll].y*Elem[N+Ncoll].y)/Elem[N+Ncoll].z));
  phi1 = atan(Elem[N+Ncoll-2].y/Elem[N+Ncoll-2].x);
  phi2 = atan(Elem[N+Ncoll-1].y/Elem[N+Ncoll-1].x);
  phi3 = atan(Elem[N+Ncoll].y/Elem[N+Ncoll].x);

  printf("theta1 2 3, phi1 2 3: %.2f %.2f %.2f %.2f %.2f %.2f\n",theta1,theta2,theta3, phi1, phi2, phi3);    
  Centre.x=0.;Centre.y=0.;Centre.z=0.;
  for (t=1;t<=N;t++){ 
    Centre.x+=Elem[t].x;
    Centre.y+=Elem[t].y;
    Centre.z+=Elem[t].z;
  }
  Centre.x/=(double)(N);
  Centre.y/=(double)(N);
  Centre.z/=(double)(N);
 
  printf("(%2.3f %2.3f %2.3f)\n",Centre.x,Centre.y,Centre.z);

  set_nb();
  set_triangles();
  //painter();- AS changed
  painter2();
}
/*---------------------------*/ 
void set_Cells(void){
  int i;
  int x,y,z;
  double facx,facy,facz;

  S.Lcelx=(int)(S.sidex/S.rb);
  S.Lcely=(int)(S.sidey/S.rb); 
  S.Lcelz=(int)(S.sidez/S.rb);
  if (S.Lcelz<S.rb) S.Lcelz=S.rb;
  S.Ncel=S.Lcelx*S.Lcely*S.Lcelz;
 
  S.celsidex=S.sidey/S.Lcelx;
  S.Icelsidex=1./S.celsidex;
  S.celsidey=S.sidey/S.Lcely;
  S.Icelsidey=1./S.celsidey;
  S.celsidez=S.sidez/S.Lcelz;
  S.Icelsidez=1./S.celsidez;
  

  facx=(double)(S.Lcelx)/S.sidex;
  facy=(double)(S.Lcely)/S.sidey;
  facz=(double)(S.Lcelz)/S.sidez;

  for (i=0;i<S.Ncel;i++){
    Cella[i].begin=0;
  }
  
  for (i=1;i<=N;i++){
    x = (int)((Elem[i].x+S.side2x)*facx);
    y = (int)((Elem[i].y+S.side2y)*facy);
    z = (int)((Elem[i].z+S.sidez/2.)*facz);
    
    Elem[i].cellID=x+y*S.Lcelx+z*S.Lcelx*S.Lcely;
    Elem[i].after=Cella[Elem[i].cellID].begin;
    Elem[ Elem[i].after ].before=i;
    Cella[ Elem[i].cellID].begin=i;
  }
  set_neighbors();
}
/*---------------------------*/
void set_neighbors(void){
  int k,s,i,j;
  int x,y,z;
  int x_P,x_M;
  int y_P,y_M;
  int z_P,z_M;
  
  for (x=0;x<S.Lcelx;x++){
    for (y=0;y<S.Lcely;y++){
      for (z=0;z<S.Lcelz;z++){
	
	k=x+y*S.Lcelx+z*S.Lcelx*S.Lcely;

	if (x==S.Lcelx-1) x_P=0; else x_P=x+1;
	if (x==0) x_M=S.Lcelx-1; else x_M=x-1;
	if (y==S.Lcely-1) y_P=0; else y_P=y+1;
	if (y==0) y_M=S.Lcely-1; else y_M=y-1;
	if (z==S.Lcelz-1) z_P=0; else z_P=z+1;
	if (z==0) z_M=S.Lcelz-1; else z_M=z-1; 
	
	Cella[k].ngb[0]=x+(y)*S.Lcelx+(z)*S.Lcelx*S.Lcely;
	Cella[k].ngb[1]=x_P+(y)*S.Lcelx+(z)*S.Lcelx*S.Lcely;
	Cella[k].ngb[2]=x_M+(y)*S.Lcelx+(z)*S.Lcelx*S.Lcely;
	Cella[k].ngb[3]=x+(y_P)*S.Lcelx+(z)*S.Lcelx*S.Lcely;
	Cella[k].ngb[4]=x+(y_M)*S.Lcelx+(z)*S.Lcelx*S.Lcely;
	Cella[k].ngb[5]=x_P+(y_P)*S.Lcelx+(z)*S.Lcelx*S.Lcely;
	Cella[k].ngb[6]=x_M+(y_P)*S.Lcelx+(z)*S.Lcelx*S.Lcely;
	Cella[k].ngb[7]=x_M+(y_M)*S.Lcelx+(z)*S.Lcelx*S.Lcely;
	Cella[k].ngb[8]=x_P+(y_M)*S.Lcelx+(z)*S.Lcelx*S.Lcely;
	
	Cella[k].ngb[9 ]=x_P+(y)*S.Lcelx+(z_P)*S.Lcelx*S.Lcely;
	Cella[k].ngb[10]=x_M+(y)*S.Lcelx+(z_P)*S.Lcelx*S.Lcely;
	Cella[k].ngb[11]=x+(y_P)*S.Lcelx+(z_P)*S.Lcelx*S.Lcely;
	Cella[k].ngb[12]=x+(y_M)*S.Lcelx+(z_P)*S.Lcelx*S.Lcely;
	Cella[k].ngb[13]=x_P+(y_P)*S.Lcelx+(z_P)*S.Lcelx*S.Lcely;
        Cella[k].ngb[14]=x_M+(y_P)*S.Lcelx+(z_P)*S.Lcelx*S.Lcely;
	Cella[k].ngb[15]=x_M+(y_M)*S.Lcelx+(z_P)*S.Lcelx*S.Lcely;
	Cella[k].ngb[16]=x_P+(y_M)*S.Lcelx+(z_P)*S.Lcelx*S.Lcely;
	
	Cella[k].ngb[17]=x_P+(y)*S.Lcelx+(z_M)*S.Lcelx*S.Lcely;
	Cella[k].ngb[18]=x_M+(y)*S.Lcelx+(z_M)*S.Lcelx*S.Lcely;
	Cella[k].ngb[19]=x+(y_P)*S.Lcelx+(z_M)*S.Lcelx*S.Lcely;
	Cella[k].ngb[20]=x+(y_M)*S.Lcelx+(z_M)*S.Lcelx*S.Lcely;
	Cella[k].ngb[21]=x_P+(y_P)*S.Lcelx+(z_M)*S.Lcelx*S.Lcely;
	Cella[k].ngb[22]=x_M+(y_P)*S.Lcelx+(z_M)*S.Lcelx*S.Lcely;
	Cella[k].ngb[23]=x_M+(y_M)*S.Lcelx+(z_M)*S.Lcelx*S.Lcely;
	Cella[k].ngb[24]=x_P+(y_M)*S.Lcelx+(z_M)*S.Lcelx*S.Lcely;
	
	Cella[k].ngb[25]=x+(y)*S.Lcelx+(z_P)*S.Lcelx*S.Lcely;
	Cella[k].ngb[26]=x+(y)*S.Lcelx+(z_M)*S.Lcelx*S.Lcely;
       
      }
    }
  }  
}
/*---------------------------*/
void set_Verlet(){  
  int i,j,neigh,k;
  double r2;
  double dx,dy,dz;
  
  for (i=1;i<=N;i++){
    Elem[i].Nverl=0;
    Elem[i].Nbond=0;
    
    Elem[i].Vx=Elem[i].x;
    Elem[i].Vy=Elem[i].y;
    Elem[i].Vz=Elem[i].z;
  }
  
  for (i=1;i<=N;i++){
    for (j=0;j<=26;j++){  // 26 is the number of neighboring cells to 1 cell
      neigh=Cella[ Elem[i].cellID ].ngb[j];
      k=Cella[neigh].begin;     
      while (k!=0){
        if (k!=i && k>=1){
	  dx=Elem[i].x-Elem[k].x;
	  dy=Elem[i].y-Elem[k].y;
	  dz=Elem[i].z-Elem[k].z;
	  
          r2=(dx*dx+ dy*dy+ dz*dz);
          if (r2<S.rv2){
	    Elem[i].Nverl++; 
	    Elem[i].verl[Elem[i].Nverl]=k;
          }
	  if (r2<S.rb2){
	    Elem[i].Nbond++; 
	    Elem[i].bond[Elem[i].Nbond]=k;
          }
	  
        }
	k=Elem[k].after;
      }
    }  
  }   
}
/*---------------------------*/  
int Interaction_1p(int k){
  int j;
  double dx,dy,dz,r2,rs,dd;
  
  for (j=1;j<=Elem[k].Nverl;j++){
    dx=Elem[k].x-Elem[ Elem[k].verl[j] ].x;
    dy=Elem[k].y-Elem[ Elem[k].verl[j] ].y;
    dz=Elem[k].z-Elem[ Elem[k].verl[j] ].z;
   
    r2=(dx*dx+ dy*dy+ dz*dz);
   
    if (r2<1) return 1;  // contact-reject
    
  }  

  return 0;
}
/*---------------------------*/  
int Interaction_ALL(void){
  int k,j;
  double dx,dy,dz,r2,rs,dd;
 
  
  for (k=1;k<=N;k++){
    for (j=1;j<=Elem[k].Nverl;j++){

      dx=Elem[k].x-Elem[ Elem[k].verl[j] ].x;
      dy=Elem[k].y-Elem[ Elem[k].verl[j] ].y;
      dz=Elem[k].z-Elem[ Elem[k].verl[j] ].z;
          
      r2=(dx*dx+ dy*dy+ dz*dz);
      

      if (r2<1) return 1;  // contact-reject
    }
  }  
  return 0;
}
/*---------------------------*/  
void MC_xyz(int oo){
  double dx,dy,dz,dx1,dy1,dz1,facx,facy,facz;
  double oldX,oldY,oldZ;
  double r2=0,r21=0,rs,rs1,DV1,DV2,DA1,DA2;
  int touch=0,ee,k;
  int i,s,x,y,z,oCel,nCel;
  double en,eo,bind_o,bind_n;
  double VOLD,AOLD;

  ee=(int)(drand48()*3.); // from 0 to 2//
  // pick a random vertex "ee" from the random triangle "oo"
 
  k=Tr[oo].v[ee];

  dx=Elem[k].x-Elem[k].Vx;
  dy=Elem[k].y-Elem[k].Vy;
  dz=Elem[k].z-Elem[k].Vz;
 
  r2=(dx*dx+ dy*dy+ dz*dz);
  
  rs=sqrt(r2);
  
  if (rs>S.Vgap){
    set_Verlet();
  }
 
 //coll-bead verlet list: 
  dx1=Elem[k].x-Elem[k].vx;
  dy1=Elem[k].y-Elem[k].vy;
  dz1=Elem[k].z-Elem[k].vz;
  r21=(dx1*dx1+ dy1*dy1+ dz1*dz1);
  rs1=sqrt(r21);
 
  if (rs1>((Coll.rv-Coll.rc)*0.5)) {
  b_coll_vlist();
  }

  oldX=Elem[k].x;
  oldY=Elem[k].y;
  oldZ=Elem[k].z;
  
  oCel=Elem[k].cellID;

  eo=bending_1p(oo,ee); // old energy (bending?)
  bind_o=bind_1p(k); 
  DV1=volume_change(oo,ee);
  VOLD=S.Svol;
   
  //S.EnergyTotal = S.kappa*(eo)+S.pres*(VOLD)+(bind_o);

  //S.BindTotal = bind_o;
  //S.BendTotal = S.kappa*(eo);
  //S.Work = S.pres*(VOLD);
    
  //painter_en(bind_o)
      
  Elem[k].x+=S.eps_mc*(0.5-drand48());
  Elem[k].y+=S.eps_mc*(0.5-drand48());
  Elem[k].z+=S.eps_mc*(0.5-drand48());
 
  // check no bonds are broken by this move
  for (i=1;i<=Elem[k].Ng;i++){
    s=Elem[k].ng[i];
    
    dx=Elem[s].x-Elem[k].x;
    dy=Elem[s].y-Elem[k].y;
    dz=Elem[s].z-Elem[k].z;
    
    r2=(dx*dx+ dy*dy+ dz*dz);
    if (r2 > cutoff2 ){ // reject here
      Elem[k].x=oldX;
      Elem[k].y=oldY;
      Elem[k].z=oldZ;
      return;
    }
  }
  // end check NO broken bonds


  dx=Elem[k].x-Elem[k].Vx;
  dy=Elem[k].y-Elem[k].Vy;
  dz=Elem[k].z-Elem[k].Vz;
  
  r2=(dx*dx+ dy*dy+ dz*dz);
  
  rs=sqrt(r2);
  if (rs>S.Vgap){
    set_Verlet();
  } 
 
  touch=Interaction_1p(k);
  
  if (touch==1){  // reject
    Elem[k].x=oldX;
    Elem[k].y=oldY;
    Elem[k].z=oldZ;
    return ;
  }

  dx1=Elem[k].x-Elem[k].vx; //coll-bead verlet list
  dy1=Elem[k].y-Elem[k].vy;
  dz1=Elem[k].z-Elem[k].vz;
  r21=(dx1*dx1+ dy1*dy1+ dz1*dz1);
  rs1=sqrt(r21);

  if (rs1>((Coll.rv-Coll.rc)*0.5)) {
  b_coll_vlist();
  }
    
  bind_n=bind_1p(k);
 
   if (bind_n==999){  // reject
    Elem[k].x=oldX;
    Elem[k].y=oldY;
    Elem[k].z=oldZ;
    return ;
  }
 
 else{ // accept
    en=bending_1p(oo,ee);
    
   DV2=volume_change(oo,ee);
    S.Svol=DV2;
    
    if (drand48()>exp(-S.kappa*(en-eo)-S.pres*(DV2-DV1)-(bind_n-bind_o)) ){
      Elem[k].x=oldX;
      Elem[k].y=oldY;
      Elem[k].z=oldZ;
      S.Svol=VOLD;
      return ;
    }
    
    S.EbendChange += S.kappa*(en-eo);
    S.EbindChange += (bind_n-bind_o);
    S.EstretchChange += S.pres*(DV2-DV1);
    
    Centre.x=Centre.x-oldX/(double)(N)+Elem[k].x/(double)(N);
    Centre.y=Centre.y-oldY/(double)(N)+Elem[k].y/(double)(N);
    Centre.z=Centre.z-oldZ/(double)(N)+Elem[k].z/(double)(N);

    
    facx=(double)(S.Lcelx)/S.sidex;
    facy=(double)(S.Lcely)/S.sidey;
    facz=(double)(S.Lcelz)/S.sidez;
    
    x=(int)((Elem[k].x+S.side2x)*facx);
    y=(int)((Elem[k].y+S.side2y)*facy);
    z=(int)((Elem[k].z+S.sidez/2.)*facz);

    
    nCel=x+y*S.Lcelx+z*S.Lcelx*S.Lcely;
    
    if (nCel!=oCel){
      if (Cella[oCel].begin==k){
	Cella[oCel].begin = Elem[k].after;
      }
      else{
	Elem[Elem[k].before].after=Elem[k].after;
	Elem[Elem[k].after].before=Elem[k].before;
      } 
      Elem[k].cellID=nCel;
      Elem[k].after=Cella[nCel].begin;
      Cella[nCel].begin=k;
      Elem[Elem[k].after].before=k;
    }
  } 
}
/*---------------------------*/  
void read_Elements(void){
  int i,j,npart,trash,t;
  
  printf("\n   READ ELEMENTS from file  ''%s''\n\n\n",readname);
 

  for (i=1;i<=N;i++){  
    fscanf(read,"%lf %lf %lf %d\n"
	   ,&Elem[i].x,&Elem[i].y,&Elem[i].z,&Elem[i].Ng);
    for (j=1;j<=Elem[i].Ng;j++){
      fscanf(read,"%d %d %d\n",&Elem[i].ng[j],
	     &Elem[i].ngR[j],&Elem[i].ngL[j]);
    }
  }  
  fclose(read);painter2();
}
/*----------------------------------------------------*/
void write(void){
  int i,j,npart,trash;
  double tmp1,tmp2,tmp3,tmp4,tmp5;
  FILE *wr1;
  
  wr1=fopen(readname,"w"); 
  for (i=1;i<=N;i++){  
    fprintf(wr1,"%4.15lf %4.15lf %4.15lf %d\n",
	    Elem[i].x,Elem[i].y,Elem[i].z,Elem[i].Ng);   
    for (j=1;j<=Elem[i].Ng;j++){
      fprintf(wr1,"%d %d %d\n",Elem[i].ng[j],
	      Elem[i].ngR[j],Elem[i].ngL[j]);
    }
    
  }
  fclose(wr1);  
} 
/*--------------------------------------------------*/ 
void set_nb(void){
  int i,j,k,w,s,keepL,keepR;
  double dx,dy,dz,dr;
  double COS,SIN,tm_L,tm_R;
  double low=100000;
  vec3D rr[N+1];

  vec3D n,n_jk;
  ANGOLO cm;

  for (i=1;i<=N;i++){
    Elem[i].Ng=0;rr[i].x=0;rr[i].y=0;rr[i].z=0;
    for (j=1;j<=N;j++){
    
      if (i!=j){
	dx=Elem[i].x-Elem[j].x;
	dy=Elem[i].y-Elem[j].y;
	dz=Elem[i].z-Elem[j].z;
	 
	
	dr=sqrt(dx*dx+dy*dy+dz*dz);
	
	if (dr<lattice+.2){
	  if (dr<low) low=dr;
	  Elem[i].Ng++;
	  Elem[i].ng[ Elem[i].Ng ]=j;
	}

      }
    }
    if (Elem[i].Ng<5||Elem[i].Ng>7 ) {printf("PROBLEMA\n"); exit(-1);}
    if (Elem[i].Ng==5) printf("5\n");if (Elem[i].Ng==7) printf("7\n");
  }

  if (low<0.99){ //should be 1.0
    printf("min dist <1 =%f RAD=%f set it at least to %f\n",
	   low,S.rad,S.rad/low);exit(-1);
  }
 


  // set bonds triagles

  for (k=1;k<=N;k++){    
    cm=angle(Elem[k].x,Elem[k].y,Elem[k].z);
    rr[k]=rotate(Elem[k].x,Elem[k].y,Elem[k].z,cm.theta,cm.phi);
    
    for (w=1;w<=Elem[k].Ng;w++){j=Elem[k].ng[w];
    rr[j]=rotate(Elem[j].x,Elem[j].y,Elem[j].z,cm.theta,cm.phi);
    }
    for (w=1;w<=Elem[k].Ng;w++){j=Elem[k].ng[w];
       
    dx=(rr[j].x-rr[k].x);
    dy=(rr[j].y-rr[k].y);
      
      
    dr=sqrt(dx*dx+dy*dy);

    n_jk.x=dx/dr; 
    n_jk.y=dy/dr;
    tm_L=-2.;
    tm_R=-2.;

    for (i=1;i<=Elem[k].Ng;i++){
      if (i!=w){
	  
	s=Elem[k].ng[i];
	dx=(rr[s].x-rr[k].x);
	dy=(rr[s].y-rr[k].y);
	  
	  
	dr=sqrt(dx*dx+dy*dy);
	n.x=dx/dr; 
	n.y=dy/dr;
	  
	COS=n.x*n_jk.x+n.y*n_jk.y;    
	SIN=n_jk.x*n.y-n_jk.y*n.x;
      
	if (SIN>0.){ // s is "up" --> counterclockwise w.r.t. jk
	  if (COS>tm_L){keepL=i;tm_L=COS;}
	}
	  
	if (SIN<0.){ // s is "down" --> clockwise w.r.t. jk
	  if (COS>tm_R){keepR=i;tm_R=COS;}
	}
      }
    }
    Elem[k].ngR[w]=keepR;
    Elem[k].ngL[w]=keepL;
    }
    
  }
}
/*--------------------------------------------------*/ 
void set_triangles(void){
  int i,j,k,flag,t;
  int tp1,tp2,tp3;
  
  t=0;
  for (i=1;i<=N;i++){
    for (j=1;j<=Elem[i].Ng;j++){
   
      tp1=i;
      tp2=Elem[i].ng[j];      
      tp3=Elem[i].ngL[j]; 
      // tp3 tells me where is the Left neighbor in the list
      // of site "i"  
      tp3=Elem[i].ng[tp3];
       //--------------- 
      for (k=1;k<=t;k++){ 
	// verify that the triangle is not already taken
	flag=0;
	if (tp1==Tr[k].v[0] ||tp1==Tr[k].v[1] ||tp1==Tr[k].v[2]) flag++;
	if (tp2==Tr[k].v[0] ||tp2==Tr[k].v[1] ||tp2==Tr[k].v[2]) flag++;
	if (tp3==Tr[k].v[0] ||tp3==Tr[k].v[1] ||tp3==Tr[k].v[2]) flag++;
	if (flag==3) break;
      }
      if (flag!=3){
	t++;
	Tr[t].v[0]=tp1;
	Tr[t].v[1]=tp2;
	Tr[t].v[2]=tp3;
	Tr[t].Nv=3;
      }
      //---------------   
       
    }
  }

  Ntri=t;
  printf("T=%d Ntri=%d\n",t,Ntri);
  painter3();
  painter4();
  
  // Now define the neighs of each triangle 
  // the neigh of v[1] must be the triangle t[1]
  // that is opposite to it; (they must share 2 points=1 line)


  for (i=1;i<=Ntri;i++){
 
   // Ngb 0 opposite to v[0]
    tp2=Tr[i].v[1];
    tp3=Tr[i].v[2];
    for (j=1;j<=Ntri;j++){
      if (i!=j){
	flag=0;
	if (tp2==Tr[j].v[0] ||tp2==Tr[j].v[1] ||tp2==Tr[j].v[2]) flag++;
	if (tp3==Tr[j].v[0] ||tp3==Tr[j].v[1] ||tp3==Tr[j].v[2]) flag++;
	if (flag==2) {// accept
	  Tr[i].t[0]=j;break;
	}	
      }
    }

    // Ngb 1 opposite to v[1]
    tp1=Tr[i].v[0];
    tp3=Tr[i].v[2];
    for (j=1;j<=Ntri;j++){
      if (i!=j){
	flag=0;
	if (tp1==Tr[j].v[0] ||tp1==Tr[j].v[1] ||tp1==Tr[j].v[2]) flag++;
	if (tp3==Tr[j].v[0] ||tp3==Tr[j].v[1] ||tp3==Tr[j].v[2]) flag++;
	if (flag==2) {// accept
	  Tr[i].t[1]=j;break;
	}	
      }
    }

   // Ngb 2 opposite to v[2]
    tp1=Tr[i].v[0];
    tp2=Tr[i].v[1];
    for (j=1;j<=Ntri;j++){
      if (i!=j){
	flag=0;
	if (tp1==Tr[j].v[0] ||tp1==Tr[j].v[1] ||tp1==Tr[j].v[2]) flag++;
	if (tp2==Tr[j].v[0] ||tp2==Tr[j].v[1] ||tp2==Tr[j].v[2]) flag++;
	if (flag==2) {// accept
	  Tr[i].t[2]=j;break;
	}	
      }
    }

    Tr[i].Nt=3;
    printf("i=%d) %d %d %d\n",i,Tr[i].v[0],Tr[i].v[1],Tr[i].v[2]);
    printf("t1=%d) %d %d %d\n",Tr[i].t[0],Tr[ Tr[i].t[0] ].v[0],
	   Tr[ Tr[i].t[0] ].v[1],Tr[ Tr[i].t[0] ].v[2]);
   
    printf("t2=%d) %d %d %d\n",Tr[i].t[1],Tr[ Tr[i].t[1] ].v[0],
	   Tr[ Tr[i].t[1] ].v[1],Tr[ Tr[i].t[1] ].v[2]);
    
    printf("t3=%d) %d %d %d\n",Tr[i].t[2],Tr[ Tr[i].t[2] ].v[0],
	   Tr[ Tr[i].t[2] ].v[1],Tr[ Tr[i].t[2] ].v[2]);
    
    printf("\n\n");
 
  }

  

}
/*--------------------------------------------------*/ 
void switch_bond(int k){
  int i,j,Ei,Es,p,w1,w2;
  int r1,r2,r3,r4;
  int tmpK[3];
  int tmpP[3];
  int FixA,FixB,Em,En;  
  int G1,G2;
  double d2,scalar,s1,s2,s3,s4,s5;
  double Energy_O,Energy_N;
  double VOLD,DV1,DV2;
  vec3D n2A,n4A,n2B,n4B;
  vec3D n1,n3,n5,n6;

  Ei=(int)(drand48()*3);  // pick a rand vertex "s" for triangle "K"
  p=Tr[k].t[Ei];          // this define the second triangle "P" given "K"

  if (Tr[p].t[0]==k) Es=0;
  if (Tr[p].t[1]==k) Es=1;
  if (Tr[p].t[2]==k) Es=2;
  //detect what is the opposite point in P to Ei

  // now I know Ei & Es  (see figure "triangles.eps")
   
  // compute distance
  r1=Tr[k].v[Ei]; 
  r2=Tr[p].v[Es];
  r3=Tr[k].v[(Ei+1)%3];
  r4=Tr[k].v[(Ei+2)%3];

  if ( (Elem[r1].Ng>=8)||(Elem[r2].Ng>=8)) return ;
  if ( (Elem[r3].Ng<=4)||(Elem[r4].Ng<=4)) return ;
 


  d2=( (Elem[r1].x-Elem[r2].x)*(Elem[r1].x-Elem[r2].x)+
       (Elem[r1].y-Elem[r2].y)*(Elem[r1].y-Elem[r2].y)+
       (Elem[r1].z-Elem[r2].z)*(Elem[r1].z-Elem[r2].z));
  
  if (d2>cutoff2) return; // reject as the candidate link is too long.

  // Check for the convexity condition of the triangle's pair to be changed

  n2A=normale(Tr[k].v[Ei],Tr[k].v[(Ei+1)%3],Tr[k].v[(Ei+2)%3]);
  n4A=normale(Tr[p].v[Es],Tr[p].v[(Es+1)%3],Tr[p].v[(Es+2)%3]);
  scalar=n2A.x*n4A.x + n2A.y*n4A.y + n2A.z*n4A.z;
  
  if (scalar<0) return; // reject
  
  n2B=normale(Tr[k].v[Ei],Tr[k].v[(Ei+1)%3],Tr[p].v[Es]);
  n4B=normale(Tr[p].v[Es],Tr[p].v[(Es+1)%3],Tr[k].v[Ei]);
  scalar=n2B.x*n4B.x + n2B.y*n4B.y + n2B.z*n4B.z;
  
  if (scalar<0) return; // reject
  
  // Check for metropolis -- energy (Check Fig. triangles2.eps)
  FixA=Tr[k].t[(Ei+1)%3];
  G1=Tr[k].t[(Ei+2)%3];

  FixB=Tr[p].t[(Es+1)%3];
  G2=Tr[p].t[(Es+2)%3];
  
  n1=normale(Tr[FixA].v[0],Tr[FixA].v[1],Tr[FixA].v[2]);
  n3=normale(Tr[G1].v[0],Tr[G1].v[1],Tr[G1].v[2]);
  n5=normale(Tr[G2].v[0],Tr[G2].v[1],Tr[G2].v[2]);
  n6=normale(Tr[FixB].v[0],Tr[FixB].v[1],Tr[FixB].v[2]);

  s1=n1.x*n2A.x+ n1.y*n2A.y+ n1.z*n2A.z;
  s2=n3.x*n2A.x+ n3.y*n2A.y+ n3.z*n2A.z;
  s3=n4A.x*n2A.x+ n4A.y*n2A.y+ n4A.z*n2A.z;
  s4=n5.x*n4A.x+ n5.y*n4A.y+ n5.z*n4A.z;
  s5=n6.x*n4A.x+ n6.y*n4A.y+ n6.z*n4A.z;
 
  Energy_O=(5.-(s1+s2+s3+s4+s5));

 s1=n1.x*n4B.x+ n1.y*n4B.y+ n1.z*n4B.z; 
  s2=n4B.x*n2B.x+ n4B.y*n2B.y+ n4B.z*n2B.z;
  s3=n3.x*n2B.x+ n3.y*n2B.y+ n3.z*n2B.z;
  s4=n6.x*n2B.x+ n6.y*n2B.y+ n6.z*n2B.z;
  s5=n5.x*n4B.x+ n5.y*n4B.y+ n5.z*n4B.z;
 

  Energy_N=(5.-(s1+s2+s3+s4+s5));
  VOLD=S.Svol;  
  DV1=volume_change_2tr(r1,r3,r4, r2,r4,r3);
  DV2=volume_change_2tr(r1,r3,r2, r2,r4,r1);
  
  S.Svol=S.Svol+(DV2-DV1);
    
  if (drand48()>exp(-S.kappa*(Energy_N-Energy_O)-S.pres*(S.Svol-VOLD))){ 
    S.Svol=VOLD;
    return; //reject
  }
  S.EbendChange += S.kappa*(Energy_N-Energy_O);
  S.EstretchChange += S.pres*(S.Svol-VOLD);
  // these definitions will be useful later 
  if (Tr[FixA].t[0]==k) Em=0;
  if (Tr[FixA].t[1]==k) Em=1;
  if (Tr[FixA].t[2]==k) Em=2;
 
  if (Tr[FixB].t[0]==p) En=0;
  if (Tr[FixB].t[1]==p) En=1;
  if (Tr[FixB].t[2]==p) En=2;
   

  // SWITCH TRIANGLES //
 
  // first switch vertices
  Tr[k].v[(Ei+2)%3]=r2;
  Tr[p].v[(Es+2)%3]=r1;
  
  // then triangle neighbors
  tmpK[Ei]=Tr[k].t[Ei];
  tmpK[(Ei+1)%3]=Tr[k].t[(Ei+1)%3];
  tmpK[(Ei+2)%3]=Tr[k].t[(Ei+2)%3];
  
  tmpP[Es]=Tr[p].t[Es];
  tmpP[(Es+1)%3]=Tr[p].t[(Es+1)%3];
  tmpP[(Es+2)%3]=Tr[p].t[(Es+2)%3];
   
  Tr[k].t[Ei]=tmpP[(Es+1)%3];
  Tr[k].t[(Ei+1)%3]=p;

  Tr[p].t[Es]=tmpK[(Ei+1)%3];
  Tr[p].t[(Es+1)%3]=k;
   
  Tr[FixA].t[Em]=p;
  Tr[FixB].t[En]=k;



  // finally add/subtruct bonds to the vertices
  
  // add
  Elem[r1].Ng++;
  Elem[r1].ng[Elem[r1].Ng]=r2;

  Elem[r2].Ng++;
  Elem[r2].ng[Elem[r2].Ng]=r1;
  
  // subtruct
  w1=0;
  for (i=1;i<=Elem[r3].Ng;i++){if (Elem[r3].ng[i]==r4) {w1=i;break;}}
  if (w1==0) {printf("Probbb 1\n"); exit(-1);}
  if (w1==Elem[r3].Ng){Elem[r3].Ng--;goto via1;}
  Elem[r3].ng[w1]=Elem[r3].ng[ Elem[r3].Ng ];
  Elem[r3].Ng--;

 via1:
  w2=0;
  for (i=1;i<=Elem[r4].Ng;i++){if (Elem[r4].ng[i]==r3) {w2=i;break;}}
  if (w2==0) {printf("Probbb 2\n"); exit(-1);}
  if (w2==Elem[r4].Ng){Elem[r4].Ng--;goto via2;}
  Elem[r4].ng[w2]=Elem[r4].ng[ Elem[r4].Ng ];
  Elem[r4].Ng--;
  
 via2:{}
    
  
}
/*--------------------------------------------------*/ 
void painter(void){
  double dx,dy,dz;
  int p,j,s,l,r;
  char iiiii[100];
  FILE *o;

  sprintf(iiiii,"out_%2.2lf.xyz",S.pres);

  o=fopen(iiiii,"a+");

  fprintf(o,"%d\n",(N+Ncoll));
  fprintf(o,"Atoms\n");
  
 for (p=1;p<=N+Ncoll;p++){ //printf("p=%d N=%d\n",p,Elem[p].Ng);
 
   if (p<=N) fprintf(o,"%d %f %f %f\n",
			     Elem[p].type,Elem[p].x,Elem[p].y,Elem[p].z);
    else fprintf(o,"%d %f %f %f\n", Elem[p].type,Elem[p].x,Elem[p].y,Elem[p].z);
       
    
  }

  fclose(o);
  /*
  for (p=1;p<=N;p++){ //printf("p=%d N=%d\n",p,Elem[p].Ng);
    for (j=1;j<=Elem[p].Ng;j++){
      s=Elem[p].ng[j];
      r=Elem[p].ngR[j];
      r=Elem[p].ng[r];
      l=Elem[p].ngL[j];
      l=Elem[p].ng[l];

      printf("(%d)--(%d) R=%d  L=%d\n",p,s,r,l);
      
    }
   printf("\n");
  }

  */

} 
void painter_en(int i, int counter, int N, int Ncoll, int Ntri){ //, bend, bind, vol
  char iiiii[100];
  FILE *o;
  int p, pe,pc;
  double EnCol, EnBend;
    
  EnCol=0;
  EnBend=0;
  for (p=1;p<=Ntri;p++){
      for (pe=0;pe<=2;pe++){
        if (bending_edge(p,pe) == bending_edge(p,pe)){ //getting rid of nan values
        EnBend +=(S.kappa)*bending_edge(p,pe);}
      
  }}
    
  for (pc=1;pc<=Ncoll;pc++){
      EnCol += bind_1p(N+pc);
      
  }
  sprintf(iiiii,"EnOut_%2.2lf.txt",S.pres);
  o=fopen(iiiii,"a+");

  //fprintf(o,"Timestep %d\n",(i));
  //fprintf(o,"Atoms\n");
  if (counter == 0) fprintf(o,"Timestep Binding Bending press*vol EbendChange EbindChange EstretchChange springEnN springEnO springEnCh; SYSTEM kappa = %f, Ncoll =  %d, Nmem = %d\n",S.kappa,Ncoll,N);
  //fprintf(o,"%f %f %f %f\n",(energytot, bend, bind, vol);
  fprintf(o,"%d %f %f %f %f %f %f %f %f %f\n",i,EnCol,EnBend,S.pres*S.Svol,S.EbendChange,S.EbindChange,S.EstretchChange,S.springEnN,S.springEnO,S.springEnCh);
  fclose(o);
  /*
  for (p=1;p<=N;p++){ //printf("p=%d N=%d\n",p,Elem[p].Ng);
    for (j=1;j<=Elem[p].Ng;j++){
      s=Elem[p].ng[j];
      r=Elem[p].ngR[j];
      r=Elem[p].ng[r];
      l=Elem[p].ngL[j];
      l=Elem[p].ng[l];

      printf("(%d)--(%d) R=%d  L=%d\n",p,s,r,l);
      
    }
   printf("\n");
  }

  */

} 
/*----------------------------------------------------*/
void ppp(int k,int j,int p,int q){
  int i;

  printf("%d--%d  (%d,%d)\n\n",k,j,p,q);
 
  for (i=1;i<=Elem[p].Ng;i++){
    printf(" E[%d].ng(%d)=%d (%d,%d)\n",p,i,Elem[p].ng[i],
	   Elem[p].ng[ Elem[p].ngR[i] ],
	   Elem[p].ng[ Elem[p].ngL[i] ]);
  }
  
  for (i=1;i<=Elem[q].Ng;i++){
    printf(" E[%d].ng(%d)=%d (%d,%d)\n",q,i,Elem[q].ng[i],
	   Elem[q].ng[ Elem[q].ngR[i] ],
	   Elem[q].ng[ Elem[q].ngL[i] ]);
  }
}
/*----------------------------------------------------*/
double bending_1p(int oo,int ee){
  int i,j,k,vo,op,nx;
  int here;
  double en=0.,tmp;
  vec3D nn,nA,nB;
 

  k=Tr[oo].v[ee];
  // triangle "oo" vertex "ee"    
 
  here=oo;
  vo=ee;
  for (i=1;i<=Tmax;i++){
    
    nn=normale(Tr[here].v[vo], Tr[here].v[(vo+1)%3],Tr[here].v[(vo+2)%3]); 
    op=Tr[here].t[vo]; // triangle opposite to point "k"
    nA=normale(Tr[op].v[0],Tr[op].v[1],Tr[op].v[2]);
    nx=Tr[here].t[(vo+1)%3]; // triangle next (C-Clwise) to point "k"
    nB=normale(Tr[nx].v[0],Tr[nx].v[1],Tr[nx].v[2]);
    
    tmp=2.-( (nn.x*nA.x+nn.y*nA.y+nn.z*nA.z)+
	  (nn.x*nB.x+nn.y*nB.y+nn.z*nB.z) );
    en+=tmp;
    /*
    printf("--Triangle    %d) p1=%d p2=%d p3=%d \n",
	   here,Tr[here].v[0],Tr[here].v[1],Tr[here].v[2] );
    printf("--Triangle OP %d) p1=%d p2=%d p3=%d \n",
	   op,Tr[op].v[0],Tr[op].v[1],Tr[op].v[2] );
    printf("--Triangle NX %d) p1=%d p2=%d p3=%d \n\n",
	   nx,Tr[nx].v[0],Tr[nx].v[1],Tr[nx].v[2] );
    */

    here=nx;  

    
    vo=4;
    if (Tr[here].v[0]==k) {vo=0;}
    if (Tr[here].v[1]==k) {vo=1;}
    if (Tr[here].v[2]==k) {vo=2;}
 
   if (vo==4){printf("BIG PROBLEM CHECK TRIANGLES\n Triangle %d) p1=%d p2=%d p3=%d looking for(%d)\n\n",
		 here,Tr[here].v[0],Tr[here].v[1],Tr[here].v[2],k );
   painter3();painter2();exit(-1);}
    
    if (here==oo) {break;}
  }
  return en;
}

double bending_edge(int oo,int ee){
  int i,j,k,vo,op,nx;
  int here;
  double en=0.,tmp;
  vec3D nn,nA;

  k=Tr[oo].v[ee];
  // triangle "oo" vertex "ee"    
  here=oo;
  vo=ee;
    
  nn=normale(Tr[here].v[vo], Tr[here].v[(vo+1)%3],Tr[here].v[(vo+2)%3]); 
  op=Tr[here].t[vo]; // triangle opposite to point "k"
  nA=normale(Tr[op].v[0],Tr[op].v[1],Tr[op].v[2]);
    
  tmp=1.-(nn.x*nA.x+nn.y*nA.y+nn.z*nA.z);
  en+=tmp/2;
  return en;
}
/*----------------------------------------------------*/
vec3D normale(int o, int i,int j){
  double lg;
  vec3D a,b,n;
  
  a.x=Elem[i].x-Elem[o].x;
  a.y=Elem[i].y-Elem[o].y;
  a.z=Elem[i].z-Elem[o].z;
  
  if (a.x>S.side2x)    a.x-=S.sidex;
  if (a.x<-S.side2x)   a.x+=S.sidex;
  if (a.y>S.side2y)    a.y-=S.sidey;
  if (a.y<-S.side2y)   a.y+=S.sidey;
  if (a.z>S.sidez/2.)  a.z-=S.sidez;
  if (a.z<-S.sidez/2.) a.z+=S.sidez;
  
 
  b.x=Elem[j].x-Elem[o].x;
  b.y=Elem[j].y-Elem[o].y;
  b.z=Elem[j].z-Elem[o].z;
  
  if (b.x>S.side2x)    b.x-=S.sidex;
  if (b.x<-S.side2x)   b.x+=S.sidex;
  if (b.y>S.side2y)    b.y-=S.sidey;
  if (b.y<-S.side2y)   b.y+=S.sidey;
  if (b.z>S.sidez/2.)  b.z-=S.sidez;
  if (b.z<-S.sidez/2.) b.z+=S.sidez;
  
  n.x=a.y*b.z-a.z*b.y;
  n.y=a.z*b.x-a.x*b.z;
  n.z=a.x*b.y-a.y*b.x;

  lg=sqrt(n.x*n.x+
	  n.y*n.y+
	  n.z*n.z);
 
  n.x/=lg; n.y/=lg; n.z/=lg;
  //  printf("nx=%f,ny=%f nz=%f\n",n.x,n.y,n.z);
  return n;
}
/******************************************************/
vec3D rotate(double x,double y,double z,double theta, double phi)
{
  vec3D first,th;
  vec3D second;
  th.x=x;th.y=y;th.z=z;

  first=rotz(th,phi);
  second=rotx(first,theta);
  return second;
}
/******************************************************/
vec3D rotx(vec3D a,double Aangle)
{
  vec3D rx;
  
  rx.x=a.x;
  rx.y=a.y*cos(Aangle)-a.z*sin(Aangle);
  rx.z=a.y*sin(Aangle)+a.z*cos(Aangle);

  return rx;
}
/******************************************************/
vec3D rotz(vec3D a,double Aangle)
{
  vec3D rz;
  
  rz.x=a.x*cos(Aangle)-a.y*sin(Aangle);
  rz.y=a.x*sin(Aangle)+a.y*cos(Aangle);
  rz.z=a.z; 
  
  return rz;
}
/******************************************************/
ANGOLO angle(double x,double y, double z)
{
  int i;
  double length,d0,sign;
  ANGOLO pp;
  vec3D nor;
  nor.x=x;nor.y=y;nor.z=z;

  length=sqrt(nor.x*nor.x+nor.y*nor.y+nor.z*nor.z);
  nor.x/=length;
  nor.y/=length;
  nor.z/=length;
  
  
  if ( nor.z>1.0-EPSILON ) nor.z=1.-EPSILON;
  if ( nor.z<-1.+EPSILON ) nor.z=-1.+EPSILON;
  pp.theta=acos(nor.z);

  if (fabs(nor.y)<EPSILON && fabs(nor.x)<EPSILON ) {pp.phi=0.;goto here;} 
  if (fabs(nor.y)<EPSILON && fabs(nor.x)>EPSILON ) {pp.phi=M_PI/2.;}
  
  pp.phi=atan(nor.x/nor.y);
  
 here:
  {  
    if (pp.phi!=0.)
      {
	if (nor.x<0 && nor.y<0) pp.phi+=M_PI;
	if (nor.x>0 && nor.y<0) pp.phi+=M_PI;
      }
  }
  
  return pp;
}
/*--------------------------------------------------*/ 
void painter2(void){
  double dx,dy,dz;
  double dax,day,daz;
  int p,j,s,r;
  FILE *o;
  o=fopen("out.list","w");

  fprintf(o,"LIST \n");
  for (p=1;p<=N;p++){ //printf("p=%d N=%d\n",p,Elem[p].Ng);
    for (j=1;j<=Elem[p].Ng;j++){
      s=Elem[p].ng[j];
      r=Elem[p].ngL[j];
      r=Elem[p].ng[r];

      fprintf(o," appearance{ material {kd 1 diffuse 1 1 .5  alpha 1.}}\n");
      fprintf(o," OFF\n3 1 3\n");
      fprintf(o," %f %f %f\n %f %f %f\n %f %f %f\n3 0 1 2 \n",
	      Elem[p].x,Elem[p].y,Elem[p].z,
	      Elem[s].x,Elem[s].y,Elem[s].z,
	      Elem[r].x,Elem[r].y,Elem[r].z);
      
    }
   
  }

  fclose(o);
 

} 

/*--------------------------------------------------*/ 
void painter3(void){
  double dx,dy,dz;
  double dax,day,daz;
  int p,j,s,r,a,b,c;
  FILE *o;
  o=fopen("Xout.list","w");

  fprintf(o,"LIST \n");
  for (p=1;p<=Ntri;p++){ //printf("p=%d N=%d\n",p,Elem[p].Ng);
    a=Tr[p].v[0];
    b=Tr[p].v[1];
    c=Tr[p].v[2];
    
    fprintf(o," appearance{ material {kd 1 diffuse .9 .9 .6  alpha 1.}}\n");
    fprintf(o," OFF\n3 1 3\n");
    fprintf(o," %f %f %f\n %f %f %f\n %f %f %f\n3 0 1 2 \n",
	    Elem[a].x,Elem[a].y,Elem[a].z,
	    Elem[b].x,Elem[b].y,Elem[b].z,
	    Elem[c].x,Elem[c].y,Elem[c].z);
    
  }

  fclose(o);
} 
/*--------------------------------------------------*/ 
void painter4(void){
  double dx,dy,dz;
  double dax,day,daz;
  int a,b,c,p,j,s,r;
  FILE *o;
  o=fopen("Xout.xyz","w");

  fprintf(o,"%d \n\n",Ntri);
  for (p=1;p<=Ntri;p++){ //printf("p=%d N=%d\n",p,Elem[p].Ng);
    dx=(Elem[Tr[p].v[0]].x+Elem[Tr[p].v[1]].x+Elem[Tr[p].v[2]].x)/3.;
    dy=(Elem[Tr[p].v[0]].y+Elem[Tr[p].v[1]].y+Elem[Tr[p].v[2]].y)/3.;
    dz=(Elem[Tr[p].v[0]].z+Elem[Tr[p].v[1]].z+Elem[Tr[p].v[2]].z)/3.;
   
    fprintf(o,"H0 %f %f %f\n",
	    dx,dy,dz);
	            
      
  }

  fclose(o);
 

} 
/*--------------------------------------------------*/ 
double volume_change(int oo,int ee){ 
  int p,q,j,k,i;
  int here,vo;
  double h,area,ss,volume=0.;
  double l1,l2,l3,ln;
  vec3D r,u1,u2,u3,n;
  vec3D dir;

  // "oo" is the triangle and "ee" is the 
  // vertex of that triangle that has been moved

  k=Tr[oo].v[ee];
    
  here=oo;
  vo=ee;
  for (i=1;i<=Tmax;i++){    
    p=Tr[here].v[(vo+1)%3];
    q=Tr[here].v[(vo+2)%3];
    
    u1.x=Elem[p].x-Elem[k].x;
    u1.y=Elem[p].y-Elem[k].y;
    u1.z=Elem[p].z-Elem[k].z;
 
    u2.x=Elem[q].x-Elem[k].x;
    u2.y=Elem[q].y-Elem[k].y;
    u2.z=Elem[q].z-Elem[k].z;

    u3.x=Elem[q].x-Elem[p].x;
    u3.y=Elem[q].y-Elem[p].y;
    u3.z=Elem[q].z-Elem[p].z;
    
    l1=sqrt(u1.x*u1.x +u1.y*u1.y +u1.z*u1.z);
    l2=sqrt(u2.x*u2.x +u2.y*u2.y +u2.z*u2.z);
    l3=sqrt(u3.x*u3.x +u3.y*u3.y +u3.z*u3.z);
      
    ss=.5*(l1+l2+l3);
    area=(sqrt(ss*(ss-l1)*(ss-l2)*(ss-l3)));
 /* 
    n.x=u1.y*u2.z-u1.z*u2.y;
    n.y=-(u1.x*u2.z-u1.z*u2.x);
    n.z=u1.x*u2.y-u1.y*u2.x;
      
    ln=sqrt(n.x*n.x+n.y*n.y+n.z*n.z);
    n.x/=ln;n.y/=ln;n.z/=ln;
    
    dir.x=Elem[k].x-Centre.x;
    dir.y=Elem[k].y-Centre.y;
    dir.z=Elem[k].z-Centre.z;
    

    h=n.x*dir.x+n.y*dir.y+n.z*dir.z;
 */ 
    volume+=area;

    
    here=Tr[here].t[(vo+1)%3];;  
    vo=4;
    if (Tr[here].v[0]==k) {vo=0;}
    if (Tr[here].v[1]==k) {vo=1;}
    if (Tr[here].v[2]==k) {vo=2;}
 
    if (vo==4){printf("BIG PROBLEMAAAAAA!!!!\n");exit(-1);}
    if (here==oo) {break;}
  }

  return volume;
}
/*--------------------------------------------------*/ 
double volume_change_2tr(int a1,int b1,int c1,int a2,int b2,int c2){
  double h,area,ss,volume=0.;
  double l1,l2,l3,ln;
  vec3D r,u1,u2,u3,n;
  vec3D dir;

    u1.x=Elem[b1].x-Elem[a1].x;
    u1.y=Elem[b1].y-Elem[a1].y;
    u1.z=Elem[b1].z-Elem[a1].z;
 
    u2.x=Elem[c1].x-Elem[a1].x;
    u2.y=Elem[c1].y-Elem[a1].y;
    u2.z=Elem[c1].z-Elem[a1].z;

    u3.x=Elem[c1].x-Elem[b1].x;
    u3.y=Elem[c1].y-Elem[b1].y;
    u3.z=Elem[c1].z-Elem[b1].z;
    
    l1=sqrt(u1.x*u1.x +u1.y*u1.y +u1.z*u1.z);
    l2=sqrt(u2.x*u2.x +u2.y*u2.y +u2.z*u2.z);
    l3=sqrt(u3.x*u3.x +u3.y*u3.y +u3.z*u3.z);
    
    ss=.5*(l1+l2+l3);
    area=(sqrt(ss*(ss-l1)*(ss-l2)*(ss-l3)));
/*    
    n.x=u1.y*u2.z-u1.z*u2.y;
    n.y=-(u1.x*u2.z-u1.z*u2.x);
    n.z=u1.x*u2.y-u1.y*u2.x;
      
    ln=sqrt(n.x*n.x+n.y*n.y+n.z*n.z);
    n.x/=ln;n.y/=ln;n.z/=ln;
    
    dir.x=Elem[a1].x-Centre.x;
    dir.y=Elem[a1].y-Centre.y;
    dir.z=Elem[a1].z-Centre.z;
    
    h=n.x*dir.x+n.y*dir.y+n.z*dir.z;
 */   
    volume+=area;
    
    /// second part
    
    u1.x=Elem[b2].x-Elem[a2].x;
    u1.y=Elem[b2].y-Elem[a2].y;
    u1.z=Elem[b2].z-Elem[a2].z;
 
    u2.x=Elem[c2].x-Elem[a2].x;
    u2.y=Elem[c2].y-Elem[a2].y;
    u2.z=Elem[c2].z-Elem[a2].z;

    u3.x=Elem[c2].x-Elem[b2].x;
    u3.y=Elem[c2].y-Elem[b2].y;
    u3.z=Elem[c2].z-Elem[b2].z;
    
    l1=sqrt(u1.x*u1.x +u1.y*u1.y +u1.z*u1.z);
    l2=sqrt(u2.x*u2.x +u2.y*u2.y +u2.z*u2.z);
    l3=sqrt(u3.x*u3.x +u3.y*u3.y +u3.z*u3.z);
    
    ss=.5*(l1+l2+l3);
    area=(sqrt(ss*(ss-l1)*(ss-l2)*(ss-l3)));
/*    
    n.x=u1.y*u2.z-u1.z*u2.y;
    n.y=-(u1.x*u2.z-u1.z*u2.x);
    n.z=u1.x*u2.y-u1.y*u2.x;
      
    ln=sqrt(n.x*n.x+n.y*n.y+n.z*n.z);
    n.x/=ln;n.y/=ln;n.z/=ln;
    
    dir.x=Elem[a2].x-Centre.x;
    dir.y=Elem[a2].y-Centre.y;
    dir.z=Elem[a2].z-Centre.z;
    
    h=n.x*dir.x+n.y*dir.y+n.z*dir.z;
*/    
    volume+=area;

    return volume;
}
/*---------------------------------------------------------*/
void MC_xyz_coll(int coll){
int touch=0;
double dx,dy,dz,dx1,dy1,dz1,r2=0,r21=0,rs1;
double oldx,oldy,oldz,en_old,en_new;
double dtheta,dphi, springEn_new, springEn_old;
double thetacoll, phicoll, phi, theta, kspring;
double dphiN, phinew, dthetaN, thetanew;
//FILE *wrens;
    
oldx=Elem[coll].x;
oldy=Elem[coll].y;
oldz=Elem[coll].z;
if (coll == N+1){;
    theta = theta1;
    phi = phi1;}
if (coll == N+2){;
    theta = theta2;
    phi = phi2;}
if (coll == N+3){;
    theta = theta3;
    phi = phi3;}

thetacoll = atan((sqrt(oldx*oldx+oldy*oldy)/oldz));
phicoll = atan(oldy/oldx);
  
dtheta =  thetacoll-theta;
dphi = phicoll - phi;
kspring = 10000;
springEn_old = kspring*(dtheta*dtheta+dphi*dphi); //spring constant is 0.6 (equivalent of 3 in units of length)
  
  dx1=Elem[coll].x-Elem[coll].vx; //coll-bead verlet list
  dy1=Elem[coll].y-Elem[coll].vy;
  dz1=Elem[coll].z-Elem[coll].vz;
  r21=(dx1*dx1+ dy1*dy1+ dz1*dz1);
  rs1=sqrt(r21);

  if (rs1>((Coll.rv-Coll.rc)*0.5)) {
  b_coll_vlist();
  }

en_old=bind_1p(coll);

Elem[coll].x+=Coll.dis_mc*(0.5-drand48());
Elem[coll].y+=Coll.dis_mc*(0.5-drand48());
Elem[coll].z+=Coll.dis_mc*(0.5-drand48());

touch=Interaction_1coll(coll);

 if (touch==1){  // reject
    Elem[coll].x=oldx;
    Elem[coll].y=oldy;
    Elem[coll].z=oldz;
    return ;
  }
  
  dx1=Elem[coll].x-Elem[coll].vx; //coll-bead verlet list
  dy1=Elem[coll].y-Elem[coll].vy;
  dz1=Elem[coll].z-Elem[coll].vz;
  r21=(dx1*dx1+ dy1*dy1+ dz1*dz1);
  rs1=sqrt(r21);

  if (rs1>((Coll.rv-Coll.rc)*0.5)) {
  b_coll_vlist();
  }
    
thetanew = atan((sqrt(Elem[coll].x*Elem[coll].x+Elem[coll].y*Elem[coll].y)/Elem[coll].z));
phinew = atan(Elem[coll].y/Elem[coll].x);
  
dthetaN =  thetanew-theta;
dphiN = phinew - phi;
    
springEn_new = kspring*(dthetaN*dthetaN+dphiN*dphiN); //spring constant is 0.6 (equivalent of 3 in units of length)
en_new=bind_1p(coll);

if (en_new==999) {
Elem[coll].x=oldx;
Elem[coll].y=oldy;
Elem[coll].z=oldz;
return;
}
/*if (springEn>0.014) { //maximum that each of phi and theta are off by 3 degrees (0.052 rads)
Elem[coll].x=oldx;
Elem[coll].y=oldy;
Elem[coll].z=oldz;
return;
}*/
//wrens=fopen("springen.dat","a"); 
//fprintf(wrens,"(en_new-en_old): %.2f (springEn_new-springEn_old): %.2f\n",en_new-en_old,springEn_new-springEn_old);
//if (en_new>en_old) {
if (drand48()>exp(-(en_new-en_old)-(springEn_new-springEn_old))) {
Elem[coll].x=oldx;
Elem[coll].y=oldy;
Elem[coll].z=oldz;
return;
    
}
//}
S.EbindChange += (en_new-en_old);
S.springEnN +=(springEn_new);
S.springEnO +=(springEn_old);
S.springEnCh +=(springEn_new-springEn_old);//-springEn_old);
    
}
/*---------------------------------------------------------*/
double bind_1p(int k){
int j,m;
double dx,dy,dz,r2=0,invr2,invr6;
double ee=0.0;

for (j=1;j<=Elem[k].nlist;j++) {
m=Elem[k].list[j];
dx=Elem[k].x-Elem[m].x;
dy=Elem[k].y-Elem[m].y;
dz=Elem[k].z-Elem[m].z;

r2=(dx*dx+dy*dy+dz*dz);
if (r2<(0.5*(1+Coll.sig)*0.5*(1+Coll.sig))) return 999;
else{
if (r2<(Coll.rc*Coll.rc)) {
invr2=0.5*(1+Coll.sig)*0.5*(1+Coll.sig)/(r2);
invr6=invr2*invr2*invr2;
ee+=-invr6;
}
}
}
ee=D0*ee;
return ee;
}
/*---------------------------------------------------------*/
int Interaction_1coll(int k){
  int j;
  double dx,dy,dz,r2=0,rs,dd;

  for (j=N+1;j<=N+Ncoll;j++){
   if (k!=j) { 
   dx=Elem[k].x-Elem[j].x;
    dy=Elem[k].y-Elem[j].y;
    dz=Elem[k].z-Elem[j].z;

    r2=(dx*dx+ dy*dy+ dz*dz);

    if (r2<(Coll.sig*Coll.sig)) return 1;  // contact-reject

  }
  }
  return 0;
}
//-----------------b_coll_vlist--------------
void b_coll_vlist(void) {
int j,i;
double dx,dy,dz,r2;

for (i=1;i<=N+Ncoll;i++) {
Elem[i].nlist=0;
Elem[i].vx=Elem[i].x;
Elem[i].vy=Elem[i].y;
Elem[i].vz=Elem[i].z;
}

for (i=1;i<=N;i++) {
for (j=N+1;j<=N+Ncoll;j++) {

dx=Elem[i].x-Elem[j].x;
dy=Elem[i].y-Elem[j].y;
dz=Elem[i].z-Elem[j].z;

r2=(dx*dx+dy*dy+dz*dz);

if (r2<Coll.rv2) {
Elem[i].nlist+=1;
Elem[j].nlist+=1;
Elem[i].list[Elem[i].nlist]=j;
Elem[j].list[Elem[j].nlist]=i;
}
}
}
}
