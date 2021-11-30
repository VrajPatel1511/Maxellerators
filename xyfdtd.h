#include <cuda.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <malloc.h>
#include <sys/time.h>
#include <omp.h>
#define Min(X, Y) (((X) < (Y)) ? (X) : (Y))
#define Max(X, Y) (((X) > (Y)) ? (X) : (Y))
int nproc,ierr;
int tag,llim[101],ulim[101],rank;

/*--------------------------------------------------------------------*/
int i,j,n,factor,nold,mold,chc,option;
int xend,yend,xstar,ystar,xstarfix,ystarfix,yendfix,xstarol,xstarnew,ystarol,ystarnew,yendol,yendnew,inibxs,inibys,xinc,yincs,yincn; //! Mesh refinement region
double xs,ys,xe,ye,istop,bxsize,bysize; //!variables to take input of mesh refinement region
double incrx,incry; //! increament factor for mesh refinement region
double c,pi,eps0,xmu0,qe,cmasse,akb,inter_value;
double xpos,ypos;
double xd1,xd2,yd1,yd2,xxpr,timp,xxpr1;
double t,dt,ds,dte,dtm,dteds,dtmds,dtsi,DTAC,newt;
double c1,c2,c3,c4;
double E0,FREQ,OMEG;
double om2snu,QSM,crec,dnma,xfront;
double FNUM,RECOMB,EMOB,EDIF,ETEM,nu2omeg2;
double PARC,ACCEL,CORNUI,ndifmax;
double TIMD,PRESSURE,TEMP0;
double DENG0,fnucor,amu,nu_omeg,e2_epsm;
double powt,fnuimax,cene;
double difa,dife;
int IORDER,IABSOR,naccel;
int nmaxwell,nx,ny,nperdt,nx1,ny1,nec;
int isig,nlamb;
int ani,ansetup2,amsetup2,amsetup3;
int nstep,nini,inul;
double inv_nperdt,inv_c,inv_cmasse, inv_ds2;
double AA1,BB1,AA2,BB2,AA3,BB3,EE3;
double **ERMS;
double **erms2;
double **exrms;
double **exrms2;
double **eyrms;
double **eyrms2;
double **erms;
double **den;
double **dent;
double **exi;
double **eyi;
double **exi1;
double **eyi1;
double **exs;
double **eys;
double **hzi;
double **hzs;
double **vx;
double **vy;
double **ext;
double **eyt;
double **exs1;
double **exs2;
double **eys1;
double **eys2;
double *xmid;
double *ymid;
double *sgdx0;
double *sgdy0;
double *DINI;
double **DIFFUSION;
double **SIO;
double **SIOT;
double **pow1;
double **frqio;
double **puis;
double **denp;
double **frq;
int *done;
int *imid;
int *jmid;

int c_nx,c_ny,c_ani,c_n;
double c_dt,c_ds,c_dte,c_dtm,c_dteds,c_dtmds,c_dtsi,c_DTAC,c_nperdt,c_inv_nperdt;
double c_c1,c_c2,c_c3,c_c4;

double **c_ERMS;
double **c_erms2;
double **c_exrms;
double **c_exrms2;
double **c_eyrms;
double **c_eyrms2;
double **c_erms;
double **c_den;
double **c_dent;
double **c_exi;
double **c_eyi;
double **c_exi1;
double **c_eyi1;
double **c_exs;
double **c_eys;
double **c_hzi;
double **c_hzs;
double **c_vx;
double **c_vy;
double **c_ext;
double **c_eyt;
double **c_exs1;
double **c_exs2;
double **c_eys1;
double **c_eys_old;
double **c_exs_old;
double *c_xmid;
double *c_ymid;
double *c_sgdx0;
double *c_sgdy0;
double *c_DINI;
double **c_DIFFUSION;
double **c_SIO;
double **c_SIOT;
double **c_pow1;
double **c_frqio;
double **c_puis;
double **c_denp;
double **c_denpold;
double **c_exold;
double **c_eyold;
double **c_vxold;
double **c_vyold;
double **c_hzold;
//double **c_hzi;
double **c_frq;
int *c_done;
int *c_imid;
int *c_jmid;

//New Added

double **ERMSp;
//double **c_eys_old;
int c_startx;
int c_starty;
int c_endx;
int c_endy;
double **exs_old;
double **eys_old;
double **vx_old;
double **vy_old;


//Added by Vishrut

int KELEC,KRMS,nmod,NECRIR,icpling,K,c_KELEC;
double output[421][2012];   //ulimp -> llimp:ulimp ??
double slambx,slamby,tstop, radius;
int length,istart,iend;
double sine, sine1, x, x0,c_x0,x01,c_x01,sine2,sine3;
int i0,i01;
double aa,betax,betay,alpha,gamma1,extk,eytk,const7,const8;
double c_qmdt, c_aa, c_alpha, c_gamma1,c_const7,c_const8,c_const4;
double omp2x,omp2y,qmdt,const3,const4,const5x,const6x,const5y,const6y;
      
FILE *fptr, *file_xelec,*fptr2;

struct node *root_elec, *child_elec;
struct node *root_mag, *child_mag,*dev_mag;
struct node *root_den, *child_den;

double **root_mesh_elec, **child_mesh_elec;
double **root_mesh_mag, **child_mesh_mag;
double **root_mesh_den, **child_mesh_den;

//Functions
struct node
{
	int m,n,locx,locy,level;
	double **mesh;
	struct node* children[10];
	struct node *parent;
};
void ZERO();
void anim();
void animE();
void free_all();
struct node *newnode(int m, int n ,int locx, int locy, int level);
void mem_allocate();
void read_and_assign();
int dimGrid, dimBlock;
double **dev_exs;
double **dev_eys;
double dev_dtmds;
struct node * dev_root_elec,*dev_child_elec,*dev_den;
double dev_x0,dev_OMEG,dev_newt,dev_inv_c,dev_c_dt,dev_sine,dev_sine1,dev_x,dev_c,dev_c_ds;
int dev_ny,dev_nx;
double dev_xxi,dev_ds,dev_ardix,dev_yyj,dev_ardiy,dev_xd0,dev_yd0,dev_dinig;  
double * dev_sgdx0,*dev_sgdy0,*dev_DINI; 
double dev_z1,dev_z2,dev_inv_nperdt,dev_E0;
double **dev_ext,**dev_eyt,**dev_ERMSp,**dev_erms2,**dev_temp_rms,**dev_c_exs,**dev_c_eys;
int *dev_KRMS,*dev_K;
