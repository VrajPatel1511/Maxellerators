%%cu
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
int tag0,llim[101],ulim[101],rank0;

/*--------------------------------------------------------------------*/
int i,j,n,factor,nold,mold,chc,option;
int xend,yend,xstar,ystar,xstarfix,ystarfix,yendfix,xstarol,xstarnew,ystarol,ystarnew,yendol,yendnew,inibxs,inibys,xinc,yincs,yincn; //! Mesh refinement region
double xs,ys,xe,ye,istop,bxsize,bysize; //!variables to take input of mesh refinement region
double incrx,incry; //! increament factor for mesh refinement region
double inter_value;
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
double slambx,slamby,tstop;
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
double *dev_dtmds;
struct node * dev_root_elec,*dev_child_elec,*dev_den;
double *dev_x0,dev_OMEG,dev_newt,dev_inv_c,dev_c_dt,dev_sine,dev_sine1,dev_x,dev_c,dev_c_ds;
int *dev_ny,*dev_nx;
double *dev_xxi,*dev_ds,*dev_ardix,*dev_yyj,*dev_ardiy,*dev_xd0,*dev_yd0,*dev_dinig;  
double * dev_sgdx0,*dev_sgdy0,*dev_DINI; 
double dev_z1,dev_z2,*dev_inv_nperdt,*dev_E0;
double **dev_ext,**dev_eyt,**dev_ERMSp,**dev_erms2,**dev_temp_rms,**dev_c_exs,**dev_c_eys;
int *dev_KRMS,*dev_K;


//xyfdtd.cu



int tag=1;
double c=2.99795e8,pi=3.14159265,eps0=8.854e-12,xmu0=1.55063706e-6,qe=1.602176487e-19,cmasse=9.10938215e-31,akb=1.3806503e-23;
double radius=0;
int rank =0;
 int sizar=3000;
struct node *newnode(int m, int n ,int locx, int locy, int level)
{
	struct node *temp = (struct node*)malloc(sizeof(struct node));
	temp->m = m;
	temp->n = n;
	temp->locx = locx;
	temp->locy = locy;
	temp->level = level;
	return temp;
}

void mem_allocate()
{
	// ERMS = malloc(sizeof(double *) * sizar);
	// if (ERMS){
	// 	for (i = 0; i < sizar; i++){
	// 		ERMS[i] = malloc(sizeof(double) * sizar);
	// 	}
	// }

	ERMSp = (double **)malloc(sizeof(double *) * sizar);
	if (ERMSp){
		for (i = 0; i < sizar; i++){
			ERMSp[i] = (double*)malloc(sizeof(double) * sizar);
		}
	}
	erms2 = (double **)malloc(sizeof(double *) * sizar);
	if (erms2){
		for (i = 0; i < sizar; i++){
			erms2[i] = (double *)malloc(sizeof(double) * sizar);
		}
	}
	// den = malloc(sizeof(double *) * sizar);
	// if (den){
	// 	for (i = 0; i < sizar; i++){
	// 		den[i] = malloc(sizeof(double) * sizar);
	// 	}
	// }
	
	exi = (double **)malloc(sizeof(double *) * sizar);
	if (exi){
		for (i = 0; i < sizar; i++){
			exi[i] = (double *)malloc(sizeof(double) * sizar);
		}
	}
	eyi = (double **)malloc(sizeof(double *) * sizar);
	if (eyi){
		for (i = 0; i < sizar; i++){
			eyi[i] = (double *)malloc(sizeof(double) * sizar);
		}
	}
	exi1 = (double **)malloc(sizeof(double *) * sizar);
	if (exi1){
		for (i = 0; i < sizar; i++){
			exi1[i] = (double *)malloc(sizeof(double) * sizar);
		}
	}
	eyi1 = (double **)malloc(sizeof(double *) * sizar);
	if (eyi1){
		for (i = 0; i < sizar; i++){
			eyi1[i] = (double *)malloc(sizeof(double) * sizar);
		}
	}
	exs = (double **)malloc(sizeof(double *) * sizar);
	if (exs){
		for (i = 0; i < sizar; i++){
			exs[i] = (double *)malloc(sizeof(double) * sizar);
		}
	}
	eys = (double **)malloc(sizeof(double *) * sizar);
	if (eys){
		for (i = 0; i < sizar; i++){
			eys[i] = (double *)malloc(sizeof(double) * sizar);
		}
	}
	hzi = (double **)malloc(sizeof(double *) * sizar);
	if (hzi){
		for (i = 0; i < sizar; i++){
			hzi[i] = (double *)malloc(sizeof(double) * sizar);
		}
	}
	vx = (double **)malloc(sizeof(double *) * sizar);
	if (vx){
		for (i = 0; i < sizar; i++){
			vx[i] = (double *)malloc(sizeof(double) * sizar);
		}
	}
	vy = (double **)malloc(sizeof(double *) * sizar);
	if (vy){
		for (i = 0; i < sizar; i++){
			vy[i] = (double *)malloc(sizeof(double) * sizar);
		}
	}
	ext = (double **)malloc(sizeof(double *) * sizar);
	if (ext){
		for (i = 0; i < sizar; i++){
			ext[i] = (double *)malloc(sizeof(double) * sizar);
		}
	}
	eyt = (double **)malloc(sizeof(double *) * sizar);
	if (eyt){
		for (i = 0; i < sizar; i++){
			eyt[i] = (double *)malloc(sizeof(double) * sizar);
		}
	}
	exs1 = (double **)malloc(sizeof(double *) * sizar);
	if (exs1){
		for (i = 0; i < sizar; i++){
			exs1[i] = (double *)malloc(sizeof(double) * sizar);
		}
	}
	eys1 = (double **)malloc(sizeof(double *) * sizar);
	if (eys1){
		for (i = 0; i < sizar; i++){
			eys1[i] = (double *)malloc(sizeof(double) * sizar);
		}
	}
	xmid = (double *)malloc(sizeof(double) * sizar);
	ymid = (double *)malloc(sizeof(double) * sizar);
	sgdx0 = (double *)malloc(sizeof(double) * sizar);
	sgdy0 = (double *)malloc(sizeof(double) * sizar);
	DINI = (double *)malloc(sizeof(double) * sizar);
	DIFFUSION = (double **)malloc(sizeof(double *) * sizar);
	if (DIFFUSION){
		for (i = 0; i < sizar; i++){
			DIFFUSION[i] = (double *)malloc(sizeof(double) * sizar);
		}
	}
	frqio = (double **)malloc(sizeof(double *) * sizar);
	if (frqio){
		for (i = 0; i < sizar; i++){
			frqio[i] = (double *)malloc(sizeof(double) * sizar);
		}
	}
	denp = (double **)malloc(sizeof(double *) * sizar);
	if (denp){
		for (i = 0; i < sizar; i++){
			denp[i] = (double *)malloc(sizeof(double) * sizar);
		}
	}
	
	imid = (int*)malloc(sizeof(int) * sizar);
	jmid = (int*)malloc(sizeof(int) * sizar);

	exs_old = (double **)malloc(sizeof(double *) * sizar);
	if (exs_old){
		for (i = 0; i < sizar; i++){
			exs_old[i] = (double *)malloc(sizeof(double) * sizar);
		}
	}
	eys_old = (double **)malloc(sizeof(double *) * sizar);
	if (eys_old){
		for (i = 0; i < sizar; i++){
			eys_old[i] = (double *)malloc(sizeof(double) * sizar);
		}
	}
	vx_old = (double **)malloc(sizeof(double *) * sizar);
	if (vx_old){
		for (i = 0; i < sizar; i++){
			vx_old[i] = (double *)malloc(sizeof(double) * sizar);
		}
	}
	vy_old = (double **)malloc(sizeof(double *) * sizar);
	if (vy_old){
		for (i = 0; i < sizar; i++){
			vy_old[i] = (double *)malloc(sizeof(double) * sizar);
		}
	}

	root_mesh_elec = (double **)malloc(sizeof(double)*sizar);
	child_mesh_elec = (double **)malloc(sizeof(double)*sizar);
	if(root_mesh_elec && child_mesh_elec)
	{
		for (i = 0; i < sizar; ++i)
		{
			root_mesh_elec[i] = (double *)calloc(sizar, sizeof(double));
			child_mesh_elec[i] = (double *)calloc(sizar, sizeof(double));	
		}
	}

	root_mesh_mag = (double **)malloc(sizeof(double)*sizar);
	child_mesh_mag = (double **)malloc(sizeof(double)*sizar);
	if(root_mesh_mag && child_mesh_mag)
	{
		for (i = 0; i < sizar; ++i)
		{
			root_mesh_mag[i] = (double *)calloc(sizar, sizeof(double));
			child_mesh_mag[i] = (double *)calloc(sizar, sizeof(double));	
		}
	}

	root_mesh_den = (double **)malloc(sizeof(double)*sizar);
	child_mesh_den = (double **)malloc(sizeof(double)*sizar);
	if(root_mesh_den && child_mesh_den)
	{
		for (i = 0; i < sizar; ++i)
		{
			root_mesh_den[i] = (double *)calloc(sizar, sizeof(double));
			child_mesh_den[i] = (double *)calloc(sizar, sizeof(double));	
		}
	}

	c_exs = (double **)malloc(sizeof(double *) * sizar);
	if (c_exs){
		for (i = 0; i < sizar; i++){
			c_exs[i] = (double *)calloc(sizar, sizeof(double));
		}
	}
	c_eys = (double **)malloc(sizeof(double *) * sizar);
	if (c_eys){
		for (i = 0; i < sizar; i++){
			c_eys[i] = (double *)calloc(sizar, sizeof(double));
		}
	}
	c_eyi = (double **)malloc(sizeof(double *) * sizar);
	if (c_eyi){
		for (i = 0; i < sizar; i++){
			c_eyi[i] = (double *)calloc(sizar, sizeof(double));
		}
	}
	c_eyi1 = (double **)malloc(sizeof(double *) * sizar);
	if (c_eyi1){
		for (i = 0; i < sizar; i++){
			c_eyi1[i] = (double *)calloc(sizar, sizeof(double));
		}
	}

	c_erms2 = (double **)malloc(sizeof(double *) * sizar);
	if (c_erms2){
		for (i = 0; i < sizar; i++){
			c_erms2[i] = (double *)malloc(sizeof(double) * sizar);
		}
	}
	// den = malloc(sizeof(double *) * sizar);
	// if (den){
	// 	for (i = 0; i < sizar; i++){
	// 		den[i] = malloc(sizeof(double) * sizar);
	// 	}
	// }
	c_vx = (double **)malloc(sizeof(double *) * sizar);
	if (c_vx){
		for (i = 0; i < sizar; i++){
			c_vx[i] = (double *)malloc(sizeof(double) * sizar);
		}
	}
	c_vy = (double **)malloc(sizeof(double *) * sizar);
	if (c_vy){
		for (i = 0; i < sizar; i++){
			c_vy[i] = (double *)malloc(sizeof(double) * sizar);
		}
	}
	c_exi = (double **)malloc(sizeof(double *) * sizar);
	if (c_exi){
		for (i = 0; i < sizar; i++){
			c_exi[i] = (double *)malloc(sizeof(double) * sizar);
		}
	}
	c_exi1 = (double **)malloc(sizeof(double *) * sizar);
	if (c_exi1){
		for (i = 0; i < sizar; i++){
			c_exi1[i] = (double *)malloc(sizeof(double) * sizar);
		}
	}
	c_ext = (double **)malloc(sizeof(double *) * sizar);
	if (c_ext){
		for (i = 0; i < sizar; i++){
			c_ext[i] = (double *)malloc(sizeof(double) * sizar);
		}
	}
	c_eyt = (double **)malloc(sizeof(double *) * sizar);
	if (c_eyt){
		for (i = 0; i < sizar; i++){
			c_eyt[i] = (double *)malloc(sizeof(double) * sizar);
		}
	}
	c_exs1 = (double **)malloc(sizeof(double *) * sizar);
	if (c_exs1){
		for (i = 0; i < sizar; i++){
			c_exs1[i] = (double *)malloc(sizeof(double) * sizar);
		}
	}
	c_eys1 = (double **)malloc(sizeof(double *) * sizar);
	if (c_eys1){
		for (i = 0; i < sizar; i++){
			c_eys1[i] = (double *)malloc(sizeof(double) * sizar);
		}
	}
	c_xmid = (double *)malloc(sizeof(double) * sizar);
	c_ymid = (double *)malloc(sizeof(double) * sizar);
	c_sgdx0 = (double *)malloc(sizeof(double) * sizar);
	c_sgdy0 = (double *)malloc(sizeof(double) * sizar);
	c_DINI = (double *)malloc(sizeof(double) * sizar);
	c_DIFFUSION = (double **)malloc(sizeof(double *) * sizar);
	if (c_DIFFUSION){
		for (i = 0; i < sizar; i++){
			c_DIFFUSION[i] = (double *)malloc(sizeof(double) * sizar);
		}
	}
	c_frqio = (double **)malloc(sizeof(double *) * sizar);
	if (c_frqio){
		for (i = 0; i < sizar; i++){
			c_frqio[i] = (double *)malloc(sizeof(double) * sizar);
		}
	}
	c_denp = (double **)malloc(sizeof(double *) * sizar);
	if (c_denp){
		for (i = 0; i < sizar; i++){
			c_denp[i] = (double *)malloc(sizeof(double) * sizar);
		}
	}

  c_denpold = (double **)malloc(sizeof(double *) * sizar);
	if (c_denpold){
		for (i = 0; i < sizar; i++){
			c_denpold[i] = (double *)malloc(sizeof(double) * sizar);
		}
	}
  
  c_exold = (double **)malloc(sizeof(double *) * sizar);
	if (c_exold){
		for (i = 0; i < sizar; i++){
			c_exold[i] = (double *)malloc(sizeof(double) * sizar);
		}
	}
  
   c_eyold = (double **)malloc(sizeof(double *) * sizar);
	if (c_eyold){
		for (i = 0; i < sizar; i++){
			c_eyold[i] = (double *)malloc(sizeof(double) * sizar);
		}
	}
 
  c_vxold = (double **)malloc(sizeof(double *) * sizar);
	if (c_vxold){
		for (i = 0; i < sizar; i++){
			c_vxold[i] = (double *)malloc(sizeof(double) * sizar);
		}
	} 
 
  c_vyold = (double **)malloc(sizeof(double *) * sizar);
	if (c_vyold){
		for (i = 0; i < sizar; i++){
			c_vyold[i] = (double *)malloc(sizeof(double) * sizar);
		}
	}
  
  c_hzold = (double **)malloc(sizeof(double *) * sizar);
	if (c_hzold){
		for (i = 0; i < sizar; i++){
			c_hzold[i] = (double *)malloc(sizeof(double) * sizar);
		}
	}
   
	c_hzi = (double **)malloc(sizeof(double *) * sizar);
	if (c_hzi){
		for (i = 0; i < sizar; i++){
			c_hzi[i] = (double *)malloc(sizeof(double) * sizar);
		}
	}
	
	c_imid = (int*)malloc(sizeof(int) * sizar);
	c_jmid = (int*)malloc(sizeof(int) * sizar);
}

void read_and_assign()
{	
	fptr=fopen("xstart.dat","r");
	//x0b=(nx-1)*ds;
    AA1=15.0*1.0e2;
    BB1=365.0*1.0e2;
    AA2=8.8050*1.0e2;
    BB2=258.450*1.0e2;
    AA3=0.0050*1.0e2;
    BB3=200.0*1.0e2;
    EE3=3.125e3;

    char *line = NULL;
    size_t len = 0;
    getline(&line, &len, fptr);
    getline(&line, &len, fptr);
    len=10;
    getline(&line, &len, fptr);
    tstop=atof(line);
    len=0;

    getline(&line, &len, fptr);
    getline(&line, &len, fptr);
    len=10;
    getline(&line, &len, fptr);
    nlamb=atof(line);
    len=0;

    getline(&line, &len, fptr);
    getline(&line, &len, fptr);
    len=10;
    getline(&line, &len, fptr);
    slambx=atof(line);
    len=0;

    getline(&line, &len, fptr);
    getline(&line, &len, fptr);
    len=10;
    getline(&line, &len, fptr);
    slamby=atof(line);
    len=0;

    getline(&line, &len, fptr);
    getline(&line, &len, fptr);
    len=10;
    getline(&line, &len, fptr);
    E0=atof(line);
    len=0;

    getline(&line, &len, fptr);
    getline(&line, &len, fptr);
    len=10;
    getline(&line, &len, fptr);
    PRESSURE=atof(line);
    len=0;

    getline(&line, &len, fptr);
    getline(&line, &len, fptr);
    len=10;
    getline(&line, &len, fptr);
    FREQ=atof(line);
    len=0;


    getline(&line, &len, fptr);
    getline(&line, &len, fptr);
    len=10;
    getline(&line, &len, fptr);
    IABSOR=atof(line);
    len=0;

    if(nini<=0)
    	nini=1;  
    
    K=1;
    
    for(K;K<=nini;K++)
    {
	    getline(&line, &len, fptr);
	    getline(&line, &len, fptr);
	    len=10;
	    getline(&line, &len, fptr);
	    xmid[K]=atof(line);
	    // printf("xmid[%d]=%f\n",K,xmid[K] );
	    len=0;

	    getline(&line, &len, fptr);
	    getline(&line, &len, fptr);
	    len=10;
	    getline(&line, &len, fptr);
	    ymid[K]=atof(line);
	    // printf("ymid[%d]=%f\n",K,ymid[K] );
	    len=0;

	    getline(&line, &len, fptr);
	    getline(&line, &len, fptr);
	    len=10;
	    getline(&line, &len, fptr);
	    sgdx0[K]=atof(line);
	    // printf("sgdx0[%d]=%f\n",K,sgdx0[K] );
	    len=0;

	    getline(&line, &len, fptr);
	    getline(&line, &len, fptr);
	    len=10;
	    getline(&line, &len, fptr);
	    sgdy0[K]=atof(line);
	    // printf("sgdy0[%d]=%f\n",K,sgdy0[K] );
	    len=0;

	    getline(&line, &len, fptr);
	    getline(&line, &len, fptr);
	    len=10;
	    getline(&line, &len, fptr);
	    DINI[K]=atof(line);
	    len=0; 
    }
      
    getline(&line, &len, fptr);
    getline(&line, &len, fptr);
    len=10;
    getline(&line, &len, fptr);
    icpling=atof(line);
    len=0;

    getline(&line, &len, fptr);
    getline(&line, &len, fptr);
    len=10;
    getline(&line, &len, fptr);
    nmaxwell=atof(line);
    len=0;

    getline(&line, &len, fptr);
    getline(&line, &len, fptr);
    len=10;
    getline(&line, &len, fptr);
    ndifmax=atof(line);
    len=0;

    getline(&line, &len, fptr);
    getline(&line, &len, fptr);
    len=10;
    getline(&line, &len, fptr);
    naccel=atof(line);
    len=0;

    getline(&line, &len, fptr);
    getline(&line, &len, fptr);
    len=10;
    getline(&line, &len, fptr);
    NECRIR=atof(line);
    len=0; 
      
    getline(&line, &len, fptr);
    getline(&line, &len, fptr);
    len=10;
    getline(&line, &len, fptr);
    crec=atof(line);
    len=0;

    getline(&line, &len, fptr);
    getline(&line, &len, fptr);
    len=10;
    getline(&line, &len, fptr);
    cene=atof(line);
    len=0; 

    getline(&line, &len, fptr);
    getline(&line, &len, fptr);
    len=10;
    getline(&line, &len, fptr);
    factor=atof(line);
    len=0; 
    
    getline(&line, &len, fptr);
    getline(&line, &len, fptr);
    len=10;
    getline(&line, &len, fptr);
    xpos=atof(line);
    len=0; 
    
    getline(&line, &len, fptr);
    getline(&line, &len, fptr);
    len=10;
    getline(&line, &len, fptr);
    ypos=atof(line);
    len=0; 
    
    getline(&line, &len, fptr);
    getline(&line, &len, fptr);
    len=10;
    getline(&line, &len, fptr);
    xs=atof(line);
    len=0;
    
    getline(&line, &len, fptr);
    getline(&line, &len, fptr);
    len=10;
    getline(&line, &len, fptr);
    ys=atof(line);
    len=0;
    
    getline(&line, &len, fptr);
    getline(&line, &len, fptr);
    len=10;
    getline(&line, &len, fptr);
    xe=atof(line);
    len=0;
    
    getline(&line, &len, fptr);
    getline(&line, &len, fptr);
    len=10;
    getline(&line, &len, fptr);
    ye=atof(line);
    len=0;
    
    getline(&line, &len, fptr);
    getline(&line, &len, fptr);
    len=10;
    getline(&line, &len, fptr);
    incrx=atof(line);
    len=0;
    
    getline(&line, &len, fptr);
    getline(&line, &len, fptr);
    len=10;
    getline(&line, &len, fptr);
    incry=atof(line);
    len=0;
    
    getline(&line, &len, fptr);
    getline(&line, &len, fptr);
    len=10;
    getline(&line, &len, fptr);
    istop=atof(line);
    len=0;
    
    getline(&line, &len, fptr);
    getline(&line, &len, fptr);
    len=10;
    getline(&line, &len, fptr);
    bxsize=atof(line);
    len=0;
    
    getline(&line, &len, fptr);
    getline(&line, &len, fptr);
    len=10;
    getline(&line, &len, fptr);
    bysize=atof(line);
    len=0;
   
    getline(&line, &len, fptr);
    getline(&line, &len, fptr);
    len=10;
    getline(&line, &len, fptr);
    option=atof(line);
    len=0;
    
    printf("xpos=%f ypos=%f xs=%f ys=%f xe=%f ye=%f incrx=%f incry=%f istop=%f bxsize=%f bysize=%f option=%d\n",xpos,ypos,xs,ys,xe,ye,incrx,incry,istop,bxsize,bysize,option);
    fclose(fptr);
}


void ZERO()
{
    int i;
    for(i=0;i<=nx;i++)
    {   
        for(j=0;j<=ny;j++)
        {
            // ERMS[i][j]=0.0;
            ERMSp[i][j]=0.0;
            erms2[i][j]=0.0;
            // den[i][j]=0.0;
            denp[i][j]=0.0;
            exi[i][j]=0.0;
            eyi[i][j]=0.0;
            exi1[i][j]=0.0;
            eyi1[i][j]=0.0;
            exs[i][j]=0.0;
            eys[i][j]=0.0;
            hzi[i][j]=0.0;
            vx[i][j]=0.0;
            vy[i][j]=0.0;
            ext[i][j]=0.0;
            eyt[i][j]=0.0;
            exs1[i][j]=0.0;
            eys1[i][j]=0.0;
            frqio[i][j]=0.0;
            DIFFUSION[i][j]=0.0;
        }
    }
    PARC=0.0;

    for(i=0;i<=nx*factor;i++)
    {   
        for(j=0;j<=ny*factor;j++)
        {	
        	// printf("%d %d\n", i,j);
		    c_erms2[i][j]=0.0;
            // den[i][j]=0.0;
            c_denp[i][j]=0.0;
            c_denpold[i][j]=0.0;
            c_exold[i][j]=0.0;
            c_eyold[i][j]=0.0;
            c_vxold[i][j]=0.0;
            c_vyold[i][j]=0.0;
            c_hzold[i][j]=0.0;
            c_eyi[i][j]=0.0;            
            c_eyi1[i][j]=0.0;
            c_exi[i][j]=0.0;            
            c_exi1[i][j]=0.0;
            // c_ext[i][j]=0.0
            c_exs[i][j]=0.0;
            c_eys[i][j]=0.0;
            // hzs[i][j]=0.0;
            c_vx[i][j]=0.0;
            c_vy[i][j]=0.0;
            c_ext[i][j]=0.0;
            c_eyt[i][j]=0.0;
            c_exs1[i][j]=0.0;
            c_eys1[i][j]=0.0;
            c_frqio[i][j]=0.0;
            c_DIFFUSION[i][j]=0.0;
            c_hzi[i][j]=0.0;
        }
    }
}

void free_all()
{
    // free(ERMS);
    free(ERMSp);
    free(erms2);
    // free(den);
    free(exi);
    free(eyi);
    free(exi1);
    free(eyi1);
    free(exs);
    free(eys);
    free(hzi);
    free(vx);
    free(vy);
    free(ext);
    free(eyt);
    free(exs1);
    free(eys1);
    free(xmid);
    free(ymid);
    free(sgdx0);
    free(sgdy0);
    free(DINI);
    free(DIFFUSION);
    free(frqio);
    free(denp);
    free(imid);
    free(jmid);

    free(c_erms2);
    // free(den);
    free(c_exi);
    free(c_eyi);
    free(c_exi1);
    free(c_eyi1);
    free(c_exs);
    free(c_eys);
    // free(hzs);
    free(c_vx);
    free(c_vy);
    free(c_ext);
    free(c_eyt);
    free(c_exs1);
    free(c_eys1);
    free(c_xmid);
    free(c_ymid);
    free(c_sgdx0);
    free(c_sgdy0);
    free(c_DINI);
    free(c_DIFFUSION);
    free(c_frqio);
    free(c_denp);
    free(c_denpold);
    free(c_exold);
    free(c_eyold);
    free(c_vxold);
    free(c_vyold);
    free(c_hzold);
    free(c_imid);
    free(c_jmid);
    free(c_hzi);
}

void anim()
{
    //printf("In anim\n");
	double e_total;
	char fil[50];
	char* buffer = (char *)malloc(sizeof(int));
	snprintf(buffer, sizeof(buffer) - 1, "%d", ani);
	strcpy(fil,"anim/");
	strcat(fil, buffer);
	strcat(fil,".dat");

	file_xelec=fopen(fil,"w");
	
	for(j=0; j<root_den->n; j++)
	{
		for(i=0; i<root_den->m; i++)
		{
			fprintf(file_xelec,"%e ",(root_den->mesh[i][j]));
			if (root_den->mesh[i][j]!=0)
			{
				// printf("%.18f\n",child_den->mesh[i][j]);
			}
		}
		fprintf(file_xelec,"\n");
	}
	// printf("Okay, printed\n");
	fclose(file_xelec);
}

void panim()
{
    //printf("In anim\n");
	double e_total;
	char fil[50];
	char* buffer = (char *)malloc(sizeof(int));
	snprintf(buffer, sizeof(buffer) - 1, "%d", ansetup2);
	strcpy(fil,"panim/");
	strcat(fil, buffer);
	strcat(fil,".dat");

	file_xelec=fopen(fil,"w");
	
	for(j=0; j<root_den->n; j++)
	{
		for(i=0; i<root_den->m; i++)
		{
			fprintf(file_xelec,"%e ",(root_den->mesh[i][j]));
			if (root_den->mesh[i][j]!=0)
			{
				// printf("%.18f\n",child_den->mesh[i][j]);
			}
		}
		fprintf(file_xelec,"\n");
	}
	
	fclose(file_xelec);
}

void canim()
{
    //printf("In anim\n");
	double e_total;
	char fil[50];
	char* buffer = (char *)malloc(sizeof(int));
	snprintf(buffer, sizeof(buffer) - 1, "%d", ansetup2);
	strcpy(fil,"canim/");
	strcat(fil, buffer);
	strcat(fil,".dat");

	file_xelec=fopen(fil,"w");
	
	for(j=0; j<child_den->n; j++)
	{
		for(i=0; i<child_den->m; i++)
		{
			fprintf(file_xelec,"%e ",(child_den->mesh[i][j]));
			if (root_den->mesh[i][j]!=0)
			{
				// printf("%.18f\n",child_den->mesh[i][j]);
			}
		}
		fprintf(file_xelec,"\n");
	}
	
	fclose(file_xelec);
}

void canim2()
{
    //printf("In anim\n");
	double e_total;
	char fil[50];
	char* buffer = (char *)malloc(sizeof(int));
	snprintf(buffer, sizeof(buffer) - 1, "%d", amsetup2);
	strcpy(fil,"canim2/");
	strcat(fil, buffer);
	strcat(fil,".dat");

	file_xelec=fopen(fil,"w");
	
	for(j=0; j<child_den->n; j++)
	{
		for(i=0; i<child_den->m; i++)
		{
			fprintf(file_xelec,"%e ",(child_den->mesh[i][j]));
			if (root_den->mesh[i][j]!=0)
			{
				// printf("%.18f\n",child_den->mesh[i][j]);
			}
		}
		fprintf(file_xelec,"\n");
	}
	
	fclose(file_xelec);
}

void canim3()
{
    //printf("In anim\n");
	double e_total;
	char fil[50];
	char* buffer = (char *)malloc(sizeof(int));
	snprintf(buffer, sizeof(buffer) - 1, "%d", amsetup3);
	strcpy(fil,"canim3/");
	strcat(fil, buffer);
	strcat(fil,".dat");

	file_xelec=fopen(fil,"w");
	
	for(j=0; j<child_den->n; j++)
	{
		for(i=0; i<child_den->m; i++)
		{
			fprintf(file_xelec,"%e ",(child_den->mesh[i][j]));
			if (root_den->mesh[i][j]!=0)
			{
				// printf("%.18f\n",child_den->mesh[i][j]);
			}
		}
		fprintf(file_xelec,"\n");
	}
	
	fclose(file_xelec);
} 

void panimE()
{
	char fil[50];
    char* buffer = (char *)malloc(sizeof(int));
    snprintf(buffer, sizeof(buffer) - 1, "%d", ansetup2);
    strcpy(fil,"panimE/");
    strcat(fil, buffer);
    strcat(fil,".dat");

    file_xelec=fopen(fil,"w");
    
    for(j=0; j<root_elec->n; j++)
	{
		for(i=0; i<root_elec->m; i++)
		{
			fprintf(file_xelec,"%e ",(root_elec->mesh[i][j]));
		}
		fprintf(file_xelec,"\n");
	}
    fclose(file_xelec);
}

void canimE()
{
	char fil[50];
    char* buffer = (char *)malloc(sizeof(int));
    snprintf(buffer, sizeof(buffer) - 1, "%d", ansetup2);
    strcpy(fil,"canimE/");
    strcat(fil, buffer);
    strcat(fil,".dat");

    file_xelec=fopen(fil,"w");
    
    for(j=0; j<child_elec->n; j++)
	{
		for(i=0; i<child_elec->m; i++)
		{
			fprintf(file_xelec,"%e ",(child_elec->mesh[i][j]));
		}
		fprintf(file_xelec,"\n");
	}
    fclose(file_xelec);
}


void denytempo()
{
  	
	char fil[50];
  int num=1;
  char *nam1="denytempo1";
  char *nam2="denytempo2";
	char* buffer = (char *)malloc(sizeof(int));
	snprintf(buffer, sizeof(buffer) - 1, "%d", num);
  sprintf(fil,"%s.dat",nam1);

	file_xelec=fopen(fil,"a");
 //!j=0.74*ny,0.26*ny, i=(0.35+0.03)*nx
	i=((xpos-0.4)+0.03)*nx;
  j=(ypos+0.24)*ny;
  
  fprintf(file_xelec,"%lf\t%e",(t*1.0e9),(root_den->mesh[i][j]));
  fprintf(file_xelec,"\n");
  fclose(file_xelec);
  
  sprintf(fil,"%s.dat",nam2);
  file_xelec=fopen(fil,"a");
  
  i=((xpos-0.4)+0.03)*nx;
  j=(ypos-0.24)*ny;
  
  fprintf(file_xelec,"%lf\t%e",(t*1.0e9),(root_den->mesh[i][j]));
  fprintf(file_xelec,"\n");
  fclose(file_xelec);
  
}

void Eytempo()
{
	  char fil[50];
    int num=1;
    char *nam1="Eytempo1" ;
    char *nam2="Eytempo2" ;
    char *nam3="Eytempo3" ;
    char *nam4="Eytempo4" ;
    char* buffer = (char *)malloc(sizeof(int));
    snprintf(buffer, sizeof(buffer) - 1, "%d", num);
    sprintf(fil,"%s.dat",nam1);

    file_xelec=fopen(fil,"a");
    //!i=0.75*nx,0.55*nx,0.35*nx,0.15*nx, j=0.5*ny
    j=ypos*ny;
    i= (xpos)*nx;
    
    fprintf(file_xelec,"%lf\t%e",(t*1.0e9),eyt[i][j]);
    fprintf(file_xelec,"\n");
    fclose(file_xelec);
    
    sprintf(fil,"%s.dat",nam2);
    file_xelec=fopen(fil,"a");
    j=ypos*ny;
    i= (xpos-0.2)*nx;
    
    fprintf(file_xelec,"%lf\t%e",(t*1.0e9),eyt[i][j]);
    fprintf(file_xelec,"\n");
    fclose(file_xelec);
    
    sprintf(fil,"%s.dat",nam3);
    file_xelec=fopen(fil,"a");
    j=ypos*ny;
    i= (xpos-0.4)*nx;
    
    fprintf(file_xelec,"%lf\t%e",(t*1.0e9),eyt[i][j]);
    fprintf(file_xelec,"\n");
    fclose(file_xelec);
    
    sprintf(fil,"%s.dat",nam4);
    file_xelec=fopen(fil,"a");
    j=ypos*ny;
    i= (xpos-0.6)*nx;
    
    fprintf(file_xelec,"%lf\t%e",(t*1.0e9),eyt[i][j]);
    fprintf(file_xelec,"\n");
    fclose(file_xelec);
    
}

void denytempocent()
{
  	
	char fil[50];
  int num=1;
  char *nam1="denyce1";
  char *nam2="denyce2";
  char *nam3="denyce3";
  char *nam4="denyce4";
	char* buffer = (char *)malloc(sizeof(int));
	snprintf(buffer, sizeof(buffer) - 1, "%d", num);
  sprintf(fil,"%s.dat",nam1);

	file_xelec=fopen(fil,"a");
 //j=0.74*ny,0.26*ny, i=(0.35+0.03)*nx
	//i=((xpos-0.4)+0.03)*nx;
  i=((xpos))*nx;
  j=(ypos)*ny;
  
  fprintf(file_xelec,"%lf\t%e",(t*1.0e9),(root_den->mesh[i][j]));
  fprintf(file_xelec,"\n");
  fclose(file_xelec);
  
  sprintf(fil,"%s.dat",nam2);
  file_xelec=fopen(fil,"a");
  
  i=((xpos))*nx;
  j=(ypos+0.15)*ny;
  
  fprintf(file_xelec,"%lf\t%e",(t*1.0e9),(root_den->mesh[i][j]));
  fprintf(file_xelec,"\n");
  fclose(file_xelec);
  
  sprintf(fil,"%s.dat",nam3);
  file_xelec=fopen(fil,"a");
  
  i=((xpos)+0.02)*nx;
  j=(ypos+0.2)*ny;
  
  fprintf(file_xelec,"%lf\t%e",(t*1.0e9),(root_den->mesh[i][j]));
  fprintf(file_xelec,"\n");
  fclose(file_xelec);
  
  sprintf(fil,"%s.dat",nam4);
  file_xelec=fopen(fil,"a");
  
  i=((xpos-0.05))*nx;
  j=(ypos)*ny;
  
  fprintf(file_xelec,"%lf\t%e",(t*1.0e9),(root_den->mesh[i][j]));
  fprintf(file_xelec,"\n");
  fclose(file_xelec);

}

void Eytempocent()
{
    char fil[50];
    int num=1;
    char *nam1="Eyce1" ;
    char *nam2="Eyce2" ;
    char *nam3="Eyce3" ;
    char *nam4="Eyce4" ;
    char* buffer = (char *)malloc(sizeof(int));
    snprintf(buffer, sizeof(buffer) - 1, "%d", num);
    sprintf(fil,"%s.dat",nam1);

    file_xelec=fopen(fil,"a");
    
    //j=ypos*ny;
    //i= (xpos)*nx;
    i=((xpos))*nx;
    j=(ypos)*ny;
    
    fprintf(file_xelec,"%lf\t%e",(t*1.0e9),fabs(eyt[i][j]));
    fprintf(file_xelec,"\n");
    fclose(file_xelec);
    
    sprintf(fil,"%s.dat",nam2);
    file_xelec=fopen(fil,"a");
    //j=ypos*ny;
    //i= (xpos-0.2)*nx;
    i=((xpos))*nx;
    j=(ypos+0.15)*ny;
    
    fprintf(file_xelec,"%lf\t%e",(t*1.0e9),fabs(eyt[i][j]));
    fprintf(file_xelec,"\n");
    fclose(file_xelec);
    
    sprintf(fil,"%s.dat",nam3);
    file_xelec=fopen(fil,"a");
    //j=ypos*ny;
    //i= (xpos-0.4)*nx;
    i=((xpos)+0.02)*nx;
    j=(ypos+0.2)*ny;
    
    fprintf(file_xelec,"%lf\t%e",(t*1.0e9),fabs(eyt[i][j]));
    fprintf(file_xelec,"\n");
    fclose(file_xelec);
    
    sprintf(fil,"%s.dat",nam4);
    file_xelec=fopen(fil,"a");
    //j=ypos*ny;
    //i= (xpos-0.6)*nx;
    i=((xpos-0.05))*nx;
    j=(ypos)*ny;
    
    fprintf(file_xelec,"%lf\t%e",(t*1.0e9),fabs(eyt[i][j]));
    fprintf(file_xelec,"\n");
    fclose(file_xelec);
    
}

void denx()
{
    //printf("In anim\n");
	
	char fil[50];
  int num=1;
  char *nam="denx";
	char* buffer = (char *)malloc(sizeof(int));
	snprintf(buffer, sizeof(buffer) - 1, "%d", num);
	sprintf(fil,"%s.dat",nam);
  
	file_xelec=fopen(fil,"a");
 
	j=ypos*ny;
	
		for(i=0; i<root_den->m; i++)
		{
			fprintf(file_xelec,"%e ",(root_den->mesh[i][j]));
			
		}
			fprintf(file_xelec,"\n");
	
	fclose(file_xelec);
}
void deny()
{
    
	
	char fil[50];
  int num=1;
  char *nam="deny";
	char* buffer = (char *)malloc(sizeof(int));
	snprintf(buffer, sizeof(buffer) - 1, "%d", num);
  sprintf(fil,"%s.dat",nam);

	file_xelec=fopen(fil,"a");
 
	i=(xpos*nx)+4;
	for(j=0; j<root_den->n; j++)
	{
		
			fprintf(file_xelec,"%e ",(root_den->mesh[i][j]));
			
	}
		fprintf(file_xelec,"\n");
	
	
	fclose(file_xelec);
}
void diffux()
{
    
	char fil[50];
  int num=1;
  char *nam="diffux";
	char* buffer = (char *)malloc(sizeof(int));
	snprintf(buffer, sizeof(buffer) - 1, "%d", num);
	sprintf(fil,"%s.dat",nam);

	file_xelec=fopen(fil,"a");
	j=ypos*ny;
 
		for(i=0; i<root_den->m; i++)
		{
			fprintf(file_xelec,"%f ",(DIFFUSION[i][j]));
			
		}
		fprintf(file_xelec,"\n");
	
	fclose(file_xelec);
}
void diffuy()
{
    
	char fil[50];
  int num=1;
  char *nam="diffuy";
	char* buffer = (char *)malloc(sizeof(int));
	snprintf(buffer, sizeof(buffer) - 1, "%d", num);
	sprintf(fil,"%s.dat",nam);

	file_xelec=fopen(fil,"a");
   i=(xpos*nx)+4;
	
	for(j=0; j<root_den->n; j++)
	{
		
			fprintf(file_xelec,"%f ",(DIFFUSION[i][j]));
		
	}
		fprintf(file_xelec,"\n");
	
	fclose(file_xelec);
}
void ionfreqx()
{
    
	char fil[50];
  int num=1;
  char *nam="ionfrx";
	char* buffer = (char *)malloc(sizeof(int));
	snprintf(buffer, sizeof(buffer) - 1, "%d", num);
	sprintf(fil,"%s.dat",nam);

	file_xelec=fopen(fil,"a");
	j=ypos*ny;
		for(i=0; i<root_den->m; i++)
		{
			fprintf(file_xelec,"%e ",((root_den->mesh[i][j])*frqio[i][j]));
			
		}
		fprintf(file_xelec,"\n");
	
	fclose(file_xelec);
}
void ionfreqy()
{
    
	char fil[50];
  int num=1;
  char *nam="ionfry";
	char* buffer = (char *)malloc(sizeof(int));
	snprintf(buffer, sizeof(buffer) - 1, "%d", num);
	sprintf(fil,"%s.dat",nam);

	file_xelec=fopen(fil,"a");
	i=(xpos*nx)+4;
 
	for(j=0; j<root_den->n; j++)
	{
	
			fprintf(file_xelec,"%e ",((root_den->mesh[i][j])*frqio[i][j]));
			
	}
		fprintf(file_xelec,"\n");
	
	fclose(file_xelec);
}
void animE()
{
	  char fil[50];
    char* buffer = (char *)malloc(sizeof(int));
    snprintf(buffer, sizeof(buffer) - 1, "%d", ani);
    strcpy(fil,"animE/");
    strcat(fil, buffer);
    strcat(fil,".dat");

    file_xelec=fopen(fil,"w");
    
    for(j=0; j<root_elec->n; j++)
	{
		for(i=0; i<root_elec->m; i++)
		{
			fprintf(file_xelec,"%e ",(root_elec->mesh[i][j]));
		}
		fprintf(file_xelec,"\n");
	}
    fclose(file_xelec);
}
void Ermsx()
{
	  char fil[50];
    int num=1;
    char *nam="Ermsx" ;
    char* buffer = (char *)malloc(sizeof(int));
    snprintf(buffer, sizeof(buffer) - 1, "%d", num);
    sprintf(fil,"%s.dat",nam);

    file_xelec=fopen(fil,"a");
    j=ypos*ny;
    
		for(i=0; i<root_elec->m; i++)
		{
			fprintf(file_xelec,"%e ",(root_elec->mesh[i][j]));
		}
		fprintf(file_xelec,"\n");
	
    fclose(file_xelec);
}
void Ermsy()
{
	  char fil[50];
    int num=1;
    char *nam="Ermsy";
    char* buffer = (char *)malloc(sizeof(int));
    snprintf(buffer, sizeof(buffer) - 1, "%d", num);
    sprintf(fil,"%s.dat",nam);

    file_xelec=fopen(fil,"a");
    i=(xpos*nx)+4;
    
    for(j=0; j<root_elec->n; j++)
	{
		
			fprintf(file_xelec,"%e ",(root_elec->mesh[i][j]));
	}
		fprintf(file_xelec,"\n");
	
    fclose(file_xelec);
}

void bothdenEx()
{
	char fil[50];
    char* buffer = (char *)malloc(sizeof(int));
    snprintf(buffer, sizeof(buffer) - 1, "%d", ani);
    sprintf(fil,"bothdenEx/denE%d.txt",ani);

    file_xelec=fopen(fil,"w");
  	j=ypos*ny; 
		
	      for(i=0; i<root_elec->m; i++)
		{
			fprintf(file_xelec,"%d %d %e %e \n",i,j,(root_den->mesh[i][j]),(root_elec->mesh[i][j]));
		}
    fclose(file_xelec);
}

void bothdenEy()
{
	char fil[50];
    char* buffer = (char *)malloc(sizeof(int));
    snprintf(buffer, sizeof(buffer) - 1, "%d", ani);
    sprintf(fil,"bothdenEy/denE%d.txt",ani);

    file_xelec=fopen(fil,"w");
  	i=(xpos*nx)+4;
		
	       for(j=1; j<root_elec->n; j++)
		{
		 

			fprintf(file_xelec,"%d %d %e %e \n",j,i,(root_den->mesh[i][j]),(root_elec->mesh[i][j]));
		 } 
    fclose(file_xelec);
}

void parentgrid()
{
	char fil[50];
    char* buffer = (char *)malloc(sizeof(int));
    snprintf(buffer, sizeof(buffer) - 1, "%d", ani);
    sprintf(fil,"parentgrid/par%d.txt",ani);
   // strcpy(fil,"bothdenE/");
  //  strcat(fil, buffer);
  //  strcat(fil,".dat");

    file_xelec=fopen(fil,"w");
  	//i=(1*nx/4)+4; 
   fprintf(file_xelec,"i j locx locy m n\n");
  for(j=0; j<=root_elec->n; j++)
	{
		for(i=0; i<=root_elec->m; i++)
		{
			fprintf(file_xelec,"%d %d %d %d %d %d\n",i,j,(root_elec->locx),(root_elec->locy),(root_elec->m),(root_elec->n));
		}
		fprintf(file_xelec,"\n");
	}
    fclose(file_xelec);
}

void childgrid()
{
  //int factor=1;
	char fil[50];
    char* buffer = (char *)malloc(sizeof(int));
    snprintf(buffer, sizeof(buffer) - 1, "%d", ani);
    sprintf(fil,"childgrid/chi%d.txt",ani);
   // strcpy(fil,"bothdenE/");
  //  strcat(fil, buffer);
  //  strcat(fil,".dat");

    file_xelec=fopen(fil,"w");
  //	i=(1*nx/4)+4; 
  fprintf(file_xelec,"i j xlen ylen locx locy m n\n");
   //for(j=0; j<=(root_elec->children[0])->n; j++)
    for(j=0; j<=child_elec->n; j+=1)
	{
		for(i=0; i<=child_elec->m; i+=1)
		{
		        
     // fprintf(file_xelec,"%d %d %f %f %d %d %d %d\n",i,j,(root_elec->children[0])->locx+i/4.0,(root_elec->children[0])->locy+j/4.0,(root_elec->children[0])->locx,(root_elec->children[0])->locy,(root_elec->children[0])->m, (root_elec->children[0])->n);
     fprintf(file_xelec,"%d %d %d %d %d %d %d %d\n",i,j, (child_elec->locx)*factor+i, (child_elec->locy)*factor+j, child_elec->locx, child_elec->locy, child_elec->m, child_elec->n);
		}
		fprintf(file_xelec,"\n");
	}
    fclose(file_xelec);
}

//! the following C-code numerically solves the Maxwell-Plasma continuity equation
//! to model the HPM air/gas breakdown and induced plasma and its interaction with the EM wave (high power HF microwave)
//! The code provides analysis of the phenomenon in terms of the plasmoid-->streamer growth --> discrete filaments 
//! propagates towards the Microwave source.
//! The study reduces the huge simulation time by implementing a staic Mesh Refinement (MR) and further a dynamic MR technique 
//! as it's extension.
//! various subroutines are used for calcuting the Maxwells and Plasma continuity equation both for the coarse and Fine mesh 

//! sequence of files that are used in the building the executable file static_amr 
//! main file static_amr.c, the inputs (the dynamic arrays, the inputs are read from file xstart.dat as well as constants) 
//! in xyfdtd.c
//! The postproceesing subroutines are written in the xyfdtd.c file.
//! All the global variables and arrays are declared in the xyfdtd.h . 
//! the parent mesh or coarse mesh data structures represented without c_ or child_ (fine mesh) 




void SETUP(); //! has been commented in the code 
//Not used
void SETUP2();  //! currently Activated: defines the initial E-fields, the Gaussian plasma density. 
                //! Defines the initial small Mesh Refined region (Coarse and Fine). Next inside the time iteration, the Mesh expansion 
                //! uses this subroutine when certain threshold condition is satisfied (dynamic MR).
//naive implementation done
void HFIELD();  //! updates the Parent (Coarse Mesh) Magnetic field
//naive implementation done
void RMS(int k);    //! The RMS (root mean squared E-field is calculated) (Coarse mesh)
//naive implementation done
void MR_MUR(int row);   //! to implement the boundary conditions
//loop carried dependency so not possible
double FIONIZ(double EE,double PR); //! the ionization subroutine
//no loop is present
void ELEC_DENS();       //! the Electron density update subroutine (Coarse mesh)
//too many if conditions ay lead to thread divergence
void EFIELD();          //! the Electric field update subroutine   (Coarse mesh)
//boundary conition and too many variables
void child_HFIELD();    //! the child magnetic field subroutine
//naive implementation done
void child_RMS(int k);  //! child RMS E-field is calcultated here
//naive implementation done
//! avoid the following child MUR subroutines only parent MUR boundary is required for the scattered wave absorption
void child_MR_MUR_0(int row);  
void child_MR_MUR_1(int row);
void child_MR_MUR_2(int row);
void child_MR_MUR_3(int row);
//------------------------------
void child_ELEC_DENS();         //! subroutine calculates child electron density
//same as ELEC_DENS
void child_EFIELD();            //! subroutine calculates child electric field
//same as EFIELD()
//! subroutines to interpolate the Parent(Coarse mesh) data for E-,H-,Electron density to child (fine mesh)
void interpolatex(double a);   // ! along x-direction row wise ( only on the Parent-child boundary)
void interpolatey(double a);    //! along y-direction columnwise    ( only on the Parent-child boundary)
void interpolatecorners(double a); //! along the corners 
void interpolatecornersnew(double a); //! avoid
void interpolatexall(double a);     //! avoid
void interpolatexinitial(double a); //! for obtaining the initial parent gaussian density profile to child mesh
void interpolatexnew(double a);     //! interpolate all the initial parent data to child along x expansion of box (mesh)
void interpolateynew(double a);     //! interpolate all the initial parent data to child along y expansion of box (mesh)     
void interpolateyall(double a);     //! interpolate all the initial parent data to child along y (parent gaussian density)
void c2p();                         //! to transfer the child (fine mesh ) data to parent mesh

//! declare all the time variables (structure type) as well as double variables to store times 
struct timeval begin, end, total_start, total_end,program_start,program_end,t_panimstar,t_panimend,telcend;

// varaibled for noting time of Calculations for each subroutine 
double t_cal_hfield,t_cal_efield,t_cal_rms,t_cal_elec_dens,t_cal_anim,t_cal_c2p,t_cal_interpolatex,t_cal_interpolatexinitial,t_cal_interpolatexnew,t_cal_interpolateynew;
double t_cal_child_hfield,t_cal_child_efield,t_cal_child_rms,t_cal_child_elec_dens;

double t_panimcal;
double t1, t_efield_hfield,t_elec_dens,t_anim,t_rms,t_zero,t_vel_x,t_vel_y;
//! declare all the file pointers used/to be used in the code  
FILE *fptr3,*fx1,*fx2,*fx3,*fy1,*fy2,*fy3,*fxt,*frefine,*canitimxy,*velp1,*velp2,*velp3;

// __global__ void HFIELD(struct node * grid,double ** exs,double ** eys, double *dtmds)
// {
//   int i = blockIdx.x*blockDim.x+threadIdx.x;
//  int j = blockIdx.y*blockDim.y+threadIdx.y;
 
   
//     if ( i < grid->m && j<grid->n )
//     {
//         grid->mesh[i][j]+= (-(eys[i+1][j]-eys[i][j])+(exs[i][j+1]-exs[i][j]))*(*dtmds);
//     }
// }

// __global__ void RMS(struct node * root_elec,double z1,double z2,double *inv_nperdt,double **ext,double **eyt,double **ERMSp,double **erms2,int *k)
// {
//     int i = blockIdx.x*blockDim.x+threadIdx.x;
//     int j = blockIdx.y*blockDim.y+threadIdx.y;
//     if ( i < root_elec->m && j<root_elec->n )
//     {
//         z1=(ext[i][j]*ext[i][j]+ext[i-1][j]*ext[i-1][j])*.5f;   //! avg of the two scattered field is required (E_eff) for the density update
//         z2=(eyt[i][j]*eyt[i][j]+eyt[i][j-1]*eyt[i][j-1])*.5f;   //! avg of the two scattered field is required (E_eff) for the density update
//         ERMSp[i][j] = erms2[i][j];
//         erms2[i][j]=erms2[i][j]+(z1+z2)*(*inv_nperdt);     //! time updates and averages (parent)
//         if(*k==2)
//         {
//             if (erms2[i][j]<0)
//             	{
//             		printf("Alert!!\n");
//             	}
//                 root_elec->mesh[i][j] = sqrt(erms2[i][j]);  //! completes a period and then squre root the time avg data (parent)
//                 erms2[i][j]=0.0f;
//         }
//     }
// }

// __global__ void child_RMS(struct node * root_elec,double z1,double z2,double *inv_nperdt,double **ext,double **eyt,double **erms2,int *k)
// {
//     int i = blockIdx.x*blockDim.x+threadIdx.x;
//     int j = blockIdx.y*blockDim.y+threadIdx.y;
//     if ( i < root_elec->m && j<root_elec->n )
//     {
//         z1=(ext[i][j]*ext[i][j]+ext[i-1][j]*ext[i-1][j])*.5f;   //! avg of the two scattered field is required (E_eff) for the density update
//         z2=(eyt[i][j]*eyt[i][j]+eyt[i][j-1]*eyt[i][j-1])*.5f;   //! avg of the two scattered field is required (E_eff) for the density update
//         erms2[i][j]=erms2[i][j]+(z1+z2)*(*inv_nperdt);     //! time updates and averages (parent)
//         if(*k==2)
//         {
//             if (erms2[i][j]<0)
//             	{
//             		printf("Alert!!\n");
//             	}
//                 root_elec->mesh[i][j] = sqrt(erms2[i][j]);  //! completes a period and then squre root the time avg data (parent)
//                 erms2[i][j]=0.0f;
//         }
//     }
// }

__global__ void setup_init(struct node *dev_root_elec,struct node * dev_den,double *E0)
{
    int i = threadIdx.x + blockIdx.x * blockDim.x; 
    int j = threadIdx.y + blockIdx.y * blockDim.y; 
    if(i<dev_root_elec->m && j>dev_root_elec->n)
    {
        dev_den->mesh[i][j] = 0.0;
  		dev_root_elec->mesh[i][j] = (*E0)/sqrt(2.0);
    }
}

// __global__ void setup_init1(struct node * root_den,int *ny,int *nx,double *xxi,double *ds,double *ardix,double *yyj,double *ardiy,double *xd0,double *yd0,double *dinig,double * sgdx0,double *sgdy0,double *DINI, int *K)
// {
//     int i = threadIdx.x + blockIdx.x * blockDim.x; 
//     int j = threadIdx.y + blockIdx.y * blockDim.y; 
//     if(i<*nx && j<*ny)
//     {
//         *xxi=(*ds)*i;
// 	    *ardix=0.0;
// 	    if(sgdx0[1]>0)
// 	        *ardix=(-pow(((*xxi)-(*xd0)),2))/2.0/sgdx0[1]/sgdx0[1];
//         *yyj=(*ds)*j;
//         *ardiy=0.0;
//         if(sgdy0[1]>0) 
//             *ardiy=-pow((*(yyj)-(*yd0)),2)/2.0/sgdy0[*K]/sgdy0[*K];
//             *dinig=DINI[*K]*exp((*ardix)+(*ardiy));
//             if(*dinig<=1.0e13)
//                 *dinig=0;
                 
//         root_den->mesh[i][j] = root_den->mesh[i][j]+ (*dinig);
//     }
// }

int main()
{

  gettimeofday(&program_start, NULL);
	mem_allocate();        //! inside the xyfdtd.c file : all the dynamic arrays are declared 
	read_and_assign();     //! inside the xyfdtd.c file: all the inputs are read from the xstart.dat file and assigned to corresponding global varibles 
   
    
    //!------------------new  ( the portion of code to note the start time )
  time_t tkl; 
  struct tm *info;  
  char buffer10[64];
  
    gettimeofday(&begin, NULL);
    t1 = begin.tv_usec;
    gettimeofday(&total_start, NULL);
    tkl=total_start.tv_sec; //!new
    info = localtime(&tkl);//!new
    //!new
  printf("%s",asctime (info));
  strftime (buffer10, sizeof buffer10, "Today is %A, %B %d.\n", info);
  printf("%s",buffer10);
  strftime (buffer10, sizeof buffer10, "The time is %I:%M %p.\n", info);
  printf("%s",buffer10);
  //!----------------------new ( also print the starting time, the start time: for runtime of complete simulation and intermediates)
    nini=1;
    

    gettimeofday(&end, NULL);
    printf("Reading Input time: %f s\n", ((end.tv_sec - begin.tv_sec) + ((end.tv_usec - begin.tv_usec)/1000000.0)));

    //======================Caculating Constants=====================

	gettimeofday(&begin, NULL);
	t1 = begin.tv_usec;

  //! Grid sizes along x and y
	nx=nlamb*slambx;     
	ny=nlamb*slamby;
  //! peiod / number of steps to complete a period
	nperdt=2.0*nlamb;
  //! for requirement in calculations as constants
	inv_nperdt=1.0f/(double)nperdt;
	inv_c=1.0f/c;
	inv_cmasse=1.0f/cmasse;
	//!--------------
  t=0;
	ani=0;         //! not required now
  ansetup2=0;   //! only required in the code for printing simulation time for parent density plot
  amsetup2=0;   //! avoid
  amsetup3=0;   //! avoid
	//nstep=0;
  n=0;    //! number of iterations (temporal) 
	printf("nx %d, ny %d, nlamb %d, slambx %f, slamby %f\n",nx,ny,nlamb,slambx, slamby );

    //!=============================================================

    printf("Here\n");
    gettimeofday(&begin, NULL);
    ZERO();       //! initialize all the arrays to zeros 
   // SETUP();
    SETUP2();     //! to create the initial Mesh refined region, the initial gaussian density and the fileds are defined
                  //! inside do loop it is used to expand the mesh refined region dynamically (expand box)
    gettimeofday(&end, NULL);
    t_zero += ((end.tv_sec - begin.tv_sec) + ((end.tv_usec - begin.tv_usec)/1000000.0));

    gettimeofday(&end, NULL);

    printf("Initialization time: %f s\n", ((end.tv_sec - begin.tv_sec) + ((end.tv_usec - begin.tv_usec)/1000000.0)));


      ELEC_DENS();       //! the electron density update on parent 

    child_ELEC_DENS();  //! the electron density update on child
    
    printf("computing up to time  %e \n",tstop);

    //!=============================================================
    fptr=fopen("output.txt","w");      

    //--------------- Initial parameters ouput----------------------

    fprintf(fptr,"nx=%d   ny=%d\n",nx,ny);
    fprintf(fptr,"DS=%f\tC_DS=%f\t DT=%f\tC_DT=%f\n",ds,c_ds,dt*1.0e15,c_dt*1.0e15);
    fprintf(fptr,"Freq=%f   Omega=%f\n",FREQ,OMEG);
    fprintf(fptr,"Time period=%f\n",1.0/FREQ);
    fprintf(fptr,"Lambda=%f\n",3.0e8/FREQ);
       printf("fprintfs\n");
    //!put 0    
    fprintf(fptr,"Collision Freq=%f\n",FNUM);
    fprintf(fptr,"Recombination Coef=%f\n",RECOMB);
    fprintf(fptr,"Mobility=%f\n",EMOB);
    fprintf(fptr,"Electron Temp= (eV)%f\n",ETEM);
    fprintf(fptr,"DIFFUSION Coef=%f\n",EDIF);
    fprintf(fptr,"Initial Gas/Neutral density=%f\n",DENG0);
    fprintf(fptr,"Cutoff-density=%f\n",(eps0*cmasse/pow(qe,2))*(pow(OMEG,2)+pow(FNUM,2)));// check formula
    fprintf(fptr, "Tstop %.10lf\n",tstop); 
    fclose(fptr);
    printf("second bunch of fprintfs\n");
    
    //!=============================================================
    KELEC=0;    //! to print the electron-density/gas-density/gas-temperature/electron-temperature data on parents at discrete periods
    c_KELEC=0;  //! to print the electron-density/gas-density/gas-temperature/electron-temperature data on child at discrete periods
    n=0;        //! parent iteration (temporal)
    c_n=0;      //! child iterations (temporal)
    int f=0;

    //!=============================================================

    printf("some fopens\n");
    //! constants/coefficients used in E-field updates in parent mesh
    qmdt=qe*inv_cmasse*dt;  
  	aa=FNUM*dt/2.0; 
  	alpha=(1.0-aa)/(1.0+aa);
  	gamma1=1+aa;
  	const3=.50*qe*qe*inv_cmasse/eps0;
  	const4=dt*dt/4.0/gamma1;
  	const7=.25*dte*(1.0+alpha);
  	const8=1.0/2.0/gamma1;
  	//!-------------------
    //! for forward wave on parent (root)
    i0=2;
  	x0=(i0-1)*ds;
    //! for backward wave on parent (root)
    i01=nx;
    x01=(i01-1)*ds;
	  //!-------------------

	  //! constants/coefficients used in E-field updates in child mesh
    c_qmdt=qe*inv_cmasse*c_dt;
    c_aa=FNUM*c_dt/2.0; 
    c_alpha=(1.0-c_aa)/(1.0+c_aa);
    c_gamma1=1+c_aa;
    c_const4=c_dt*c_dt/4.0/c_gamma1;	
    c_const7=.25*c_dte*(1.0+c_alpha);
    c_const8=1.0/2.0/c_gamma1;
    //!-----------------------
    c_x0 = c_ds;                //! forward wave on child
    c_x01= (i01-1)*factor*c_ds; //! backward wave on child
    c_nperdt = nperdt*dt/c_dt;
    c_inv_nperdt = 1.0f/(double)c_nperdt;
    
    
    int icent,jcent,iloc1,iloc2,iloc3,iloc4,jloc1,jloc2,jloc3,k,iend,p1,ix,py1,py2,py3,py4,py5,py6,py7,py8;
    int count_cen,countx1,countx2,countx3,county1,county2,county3,county4,n_ini,lup,countrank,count_run,count_cen1;
    double velx1,velx2,velx3,vely1,vely2,vely3,tx1,tx2,tx3,ty1,ty2,ty3,lambda,tstart,t0,tnew;
    inter_value = 1.0/factor;   //! decide fraction for interpolated data
    
    icent=(int)(xpos*nx); //the x-center (interms of cell number)
    //! avoid for now 
    iloc1=(icent-(0.15*nx));iloc2=(icent-(0.25*nx));iloc3=(icent-(0.35*nx)); iloc4=(icent-(0.55*nx));
    
        
    //!-- ycenter (interms of cell number)
    jcent=(int)(ypos*ny);
    //! avoid for now
    jloc1=(jcent+(0.15*ny));jloc2=(jcent+(0.25*ny));jloc3=(jcent+(0.35*ny));
    //!-------------------
   
    p1=(int)(0.2*nlamb);
    //!--- location on y to detect streamer growth and expand box (mesh) along y----
    //! 10% top/bottom, 30% top/bottom, 40% top/bottom, 50% top/bottom,
    py1=jcent+(int)(0.1*jcent);
    py2=jcent-(int)(0.1*jcent);
    py3=jcent+(int)(0.3*jcent);
    py4=jcent-(int)(0.3*jcent); 
    py5=jcent+(int)(0.4*jcent);
    py6=jcent-(int)(0.4*jcent);
    py7=jcent+(int)(0.5*jcent);
    py8=jcent-(int)(0.5*jcent); 
   
   
    iend=(int)(istop*nx);   //!same as xs*nx (max mesh region extent towards source on left)
    lambda=c/FREQ;tstart=0; tnew=0;
    k=1;
    

    //!avoid now
    count_cen=-1;countx1=-1;countx2=-1;countx3=-1;
    //! to detect the streamer crossing a point on the line once to dynamically expand the mesh  
    county1=-1;county2=-1;county3=-1;county4=-1;
    //! avoid
    count_run=-1;count_cen1=-1;
 
    printf("icent=%d\tjcent=%d\tiloc1=%d\tiloc2=%d\tiloc3=%d\tjloc1=%djloc2=%d\tjloc3=%d\n",icent,jcent,iloc1,iloc2,iloc3,jloc1,jloc2,jloc3);
    printf("iend=%d \n",iend);
    
    //! print the initial mesh informations
    frefine=fopen("refregionxy.txt","a");
      fprintf(frefine,"time ");
      fprintf(frefine,"KELEC ");
      fprintf(frefine,"Tcellsrefx ");
      fprintf(frefine,"Tcellsrefy ");
      fprintf(frefine,"xstar ");
      fprintf(frefine,"xend ");
      fprintf(frefine,"ystar ");
      fprintf(frefine,"yend \n");
      fclose(frefine);
      
      frefine=fopen("refregionxy.txt","a");
      fprintf(frefine,"%lf ",t*1.0e9);
      fprintf(frefine,"%d ",KELEC);
      fprintf(frefine,"%d ",(xend-xstar)*factor);
      fprintf(frefine,"%d ",(yend-ystar)*factor);
      fprintf(frefine,"%d ",xstar);
      fprintf(frefine,"%d ",xend);
      fprintf(frefine,"%d ",ystar);
      fprintf(frefine,"%d \n",yend);
      fclose(frefine);
    //!--------------------------
      //! print the simulation and corresponding runtime for expanding mesh
      canitimxy=fopen("canitimxy.txt","a");
      fprintf(canitimxy,"KELEC ");
      fprintf(canitimxy,"sim.time(ns) ");
      fprintf(canitimxy,"run time(s) \n");
      fclose(canitimxy);

    /*    //! activate to print the front location for calculating front velocity
          velp1=fopen("velplot1.txt","a");
        	fprintf(velp1,"i ");  
          fprintf(velp1,"distance(lamb) ");
          fprintf(velp1,"time(ns) "); 
          fprintf(velp1,"density\n"); 
          fclose(velp1);
    */
    do
    {
        n=n+1;
                
        
      
     if(option==1)  //! activates only when option is 1 in xstart.dat ( for dynamic MR /AMR) else for only MR it is deactivated
     {  
       //!================= Increament Mesh size (dynamic MR) along x and y ======================
       if((root_den->mesh[xstar][jcent])>1.0e16) //! checks the density at the xstart of initial box exceeding threshold to expand along x
        {   
            chc=1;
            if(xstar<=iend)         //! stops growing as the xstart reaches equal or less than the iend=istop*nx (or xs*nx)
            {
              printf(" fine mesh has reached the boundary\n");
              break;
            }
            else{
                 
            printf("Entered SETUP2\n");
            printf("KELEC=%d\n",(int)(n/nperdt));
            
            SETUP2();           //! for creating expanding mesh from initial mesh along x-direction
            
           printf("leave SETUP2\n");
           printf("KELEC=%d\n",(int)(n/nperdt));
           
            frefine=fopen("refregionxy.txt","a");
            fprintf(frefine,"%lf ",t*1.0e9);
            fprintf(frefine,"%d ",KELEC);
            fprintf(frefine,"%d ",(xend-xstar)*factor);
            fprintf(frefine,"%d ",(yend-ystar)*factor);
            fprintf(frefine,"%d ",xstar);
            fprintf(frefine,"%d ",xend);
            fprintf(frefine,"%d ",ystar);
            fprintf(frefine,"%d \n",yend);
            fclose(frefine);
           }
        }
        
        for (ix=0;ix<root_den->m;ix++)
        {
            //! first threshold line (10% up/ down jcent (or ycent))
         if((root_den->mesh[ix][py1])>1.0e16 || (root_den->mesh[ix][py2])>1.0e16) //! density check at both ystart/yend of initial box exceeding threshold to expand along y
         {  
           county1+=1;
           if(county1==0){ 
            chc=2;
            
            if(ystar<=(ys*ny)|| yend>=(ye*ny) )   //! stops growing as the ystart/yend reaches equal or less than the (ys*ny)/ greater than (ye*ny)
            {
                printf(" fine mesh reached the boundary\n");
             // break;
            }
            else{
                 
            printf("Entered SETUP2\n");
            printf("KELEC=%d\n",(int)(n/nperdt));
            
            SETUP2();          //! for creating expanding mesh from initial mesh along y-direction
                        
            printf("leave SETUP2\n");
            printf("KELEC=%d\n",(int)(n/nperdt));
            
           frefine=fopen("refregionxy.txt","a");
            fprintf(frefine,"%lf ",t*1.0e9);
            fprintf(frefine,"%d ",KELEC);
            fprintf(frefine,"%d ",(xend-xstar)*factor);
            fprintf(frefine,"%d ",(yend-ystar)*factor);
            fprintf(frefine,"%d ",xstar);
            fprintf(frefine,"%d ",xend);
            fprintf(frefine,"%d ",ystar);
            fprintf(frefine,"%d \n",yend);
            fclose(frefine);
            }
           }
         }
        }
        
        for (ix=0;ix<root_den->m;ix++)
        {
          //! first threshold line (30% up/ down jcent (or ycent))
         if((root_den->mesh[ix][py3])>1.0e16 || (root_den->mesh[ix][py4])>1.0e16) //! density check at both ystart/yend of initial box exceeding threshold to expand along y
         {  
           county2+=1;
           if(county2==0){ 
             chc=2;
            
            if(ystar<=(ys*ny)|| yend>=(ye*ny) ) //! stops growing as the ystart/yend reaches equal or less than the (ys*ny)/ greater than (ye*ny)
            {
                printf(" fine mesh reached the boundary\n");
             // break;
            }
            else{
                 
            printf("Entered SETUP2\n");
            printf("KELEC=%d\n",(int)(n/nperdt));
            
            SETUP2();                    //! for creating expanding mesh from initial mesh along y-direction
                        
            printf("leave SETUP2\n");
            printf("KELEC=%d\n",(int)(n/nperdt));
            
           frefine=fopen("refregionxy.txt","a");
            fprintf(frefine,"%lf ",t*1.0e9);
            fprintf(frefine,"%d ",KELEC);
            fprintf(frefine,"%d ",(xend-xstar)*factor);
            fprintf(frefine,"%d ",(yend-ystar)*factor);
            fprintf(frefine,"%d ",xstar);
            fprintf(frefine,"%d ",xend);
            fprintf(frefine,"%d ",ystar);
            fprintf(frefine,"%d \n",yend);
            fclose(frefine);
            }
           }
         }
        }
        
         for (ix=0;ix<root_den->m;ix++)
        {
          //! first threshold line (40% up/ down jcent (or ycent))
         if((root_den->mesh[ix][py5])>1.0e16 || (root_den->mesh[ix][py6])>1.0e16)  //! density check at both ystart/yend of initial box exceeding threshold to expand along y   
         {   
            county3+=1;
           if(county3==0){ 
            chc=2;
            
            if(ystar<=(ys*ny)|| yend>=(ye*ny) )   //! stops growing as the ystart/yend reaches equal or less than the (ys*ny)/ greater than (ye*ny)
            {
                printf(" fine mesh reached the boundary\n");
             // break;
            }
            else{
                 
            printf("Entered SETUP2\n");
            printf("KELEC=%d\n",(int)(n/nperdt));
            
            SETUP2();           //! for creating expanding mesh from initial mesh along y-direction
                        
            printf("leave SETUP2\n");
            printf("KELEC=%d\n",(int)(n/nperdt));
            
           frefine=fopen("refregionxy.txt","a");
            fprintf(frefine,"%lf ",t*1.0e9);
            fprintf(frefine,"%d ",KELEC);
            fprintf(frefine,"%d ",(xend-xstar)*factor);
            fprintf(frefine,"%d ",(yend-ystar)*factor);
            fprintf(frefine,"%d ",xstar);
            fprintf(frefine,"%d ",xend);
            fprintf(frefine,"%d ",ystar);
            fprintf(frefine,"%d \n",yend);
            fclose(frefine);
            }
           }
         }
        }
        
         for (ix=0;ix<root_den->m;ix++)
        {
          //! first threshold line (50% up/ down jcent (or ycent))
         if((root_den->mesh[ix][py7])>1.0e16 || (root_den->mesh[ix][py8])>1.0e16) //! density check at both ystart/yend of initial box exceeding threshold to expand along y
         {   
            county4+=1;
           if(county4==0){ 
            chc=2;
            
            if(ystar<=(ys*ny)|| yend>=(ye*ny) )     //! stops growing as the ystart/yend reaches equal or less than the (ys*ny)/ greater than (ye*ny)
            {
                printf(" fine mesh reached the boundary\n");
             // break;
            }
            else{
                 
            printf("Entered SETUP2\n");
            printf("KELEC=%d\n",(int)(n/nperdt));
            
            SETUP2();                           //! for creating expanding mesh from initial mesh along y-direction
                      
            printf("leave SETUP2\n");
            printf("KELEC=%d\n",(int)(n/nperdt));
            
           frefine=fopen("refregionxy.txt","a");
            fprintf(frefine,"%lf ",t*1.0e9);
            fprintf(frefine,"%d ",KELEC);
            fprintf(frefine,"%d ",(xend-xstar)*factor);
            fprintf(frefine,"%d ",(yend-ystar)*factor);
            fprintf(frefine,"%d ",xstar);
            fprintf(frefine,"%d ",xend);
            fprintf(frefine,"%d ",ystar);
            fprintf(frefine,"%d \n",yend);
            fclose(frefine);
            }
           }
         }
        }
      }
      //!===================================================== 
        //!----------- velocity ----------
        /*
          velp1=fopen("velplot1.txt","a");    
      		
      			if((root_den->mesh[icent-k][jcent])>1.0e19)
            {
                         
             fprintf(velp1,"%d %f %lf %e \n",icent-k,((1.0/nlamb)*(float)k),t*1.0e9,(root_den->mesh[icent-k][jcent]));
             k+=1;
      		  }
             
          fclose(velp1);
        */
        
        //!------------------------------------
        
                
      	
        if((t>tstop)||(root_den->mesh[iend][jcent])>1.0e16) //! use to terminate the full code while simulating
        //if((t>tstop)||KELEC==10520)      //! use particular value of KELEC to stop code while experimenting               
        //if((t>tstop)||KELEC==4020)
        {  
	         
    	    printf("code has successfully ended\n");  
          gettimeofday(&total_end,NULL);   
           tkl=total_end.tv_sec; //!new
           info = localtime(&tkl);//!new
           //! gives info on the clock time when code ended
           printf("%s",asctime (info)); 
          t1 = ((total_end.tv_sec - total_start.tv_sec) + ((total_end.tv_usec - total_start.tv_usec)/1000000.0));
          break;
	     
        }
        gettimeofday(&begin, NULL);

        
        HFIELD();   //! call H-field at n   t=0
        // cudaMalloc((void**)&dev_mag, sizeof(root_mag));
        // cudaMalloc((void**)&dev_exs, sizeof(exs));
        // cudaMalloc((void**)&dev_eys, sizeof(eys));
        // cudaMalloc((void**)&dev_dtmds, sizeof(dtmds));
        // cudaMemcpy(dev_mag, root_mag, sizeof(root_mag), cudaMemcpyHostToDevice);
        // cudaMemcpy(dev_exs, exs, sizeof(exs), cudaMemcpyHostToDevice);
        // cudaMemcpy(dev_eys, eys, sizeof(eys), cudaMemcpyHostToDevice);
        // cudaMemcpy(dev_dtmds, &dtmds, sizeof(dtmds), cudaMemcpyHostToDevice);
        // HFIELD<<<(ceil(root_mag->m/32),ceil(root_mag->n/32)),(32,32)>>>(dev_mag,dev_exs,dev_eys,dev_dtmds);
        // cudaMemcpy(root_mag, dev_mag, sizeof(root_mag), cudaMemcpyDeviceToHost);
        // cudaMemcpy(hzi, dev_mag->mesh, sizeof(root_mag), cudaMemcpyDeviceToHost);

        // cudaFree(dev_mag);
        // cudaFree(dev_exs);
        // cudaFree(dev_eys);
        // cudaFree(dev_dtmds);

        EFIELD();   //! call E-field at n+1/2 time steps  t=1/2dt
        
        // cudaMalloc((void**)&dev_root_elec, sizeof(root_elec));
        // cudaMalloc((void**)&dev_x0, sizeof(double));
        // cudaMalloc((void**)&dev_OMEG, sizeof(double));
        // cudaMalloc((void**)&dev_newt, sizeof(double));
        // cudaMalloc((void**)&dev_inv_c, sizeof(double));
        // cudaMalloc((void**)&dev_c_dt, sizeof(double));
        // cudaMalloc((void**)&dev_sine, sizeof(double));
        // cudaMalloc((void**)&dev_sine1, sizeof(double));
        // cudaMalloc((void**)&dev_x, sizeof(double));
        // cudaMalloc((void**)&dev_c, sizeof(double));

        // cudaMemcpy(dev_root_elec, root_elec, sizeof(root_elec), cudaMemcpyHostToDevice);
        // cudaMemcpy(dev_x0, x0, sizeof(double), cudaMemcpyHostToDevice);
        // cudaMemcpy(dev_OMEG, OMEG, sizeof(double), cudaMemcpyHostToDevice);
        // cudaMemcpy(dev_newt, newt, sizeof(double), cudaMemcpyHostToDevice);
        // cudaMemcpy(dev_inv_c, inv_c, sizeof(double), cudaMemcpyHostToDevice);
        // cudaMemcpy(dev_inv_c, inv_c, sizeof(double), cudaMemcpyHostToDevice);
        // cudaMemcpy(dev_c, c, sizeof(double), cudaMemcpyHostToDevice);
        
        // EFIELD<<<dimGrid, dimBlock>>>(dev_elec,dev_x0,dev_OMEG,dev_newt,dev_inv_c,dev_c_dt,dev_sine,dev_sine1,dev_x,dev_c);
        
        // cudaMemcpy(sine, dev_sine, sizeof(dev_sine), cudaMemcpyDeviceToHost);
        // cudaMemcpy(sine1, dev_sine1, sizeof(dev_sine1), cudaMemcpyDeviceToHost);

        // cudaFree(dev_elec);
        // cudaFree(dev_x0);
        // cudaFree(dev_OMEG);
        // cudaFree(dev_newt);
        // cudaFree(dev_inv_c);
        // cudaFree(dev_c_dt);
        // cudaFree(dev_sine);
        // cudaFree(dev_sine1);
        // cudaFree(dev_x);
        // cudaFree(dev_c);

        gettimeofday(&end, NULL);
        t_efield_hfield += ((end.tv_sec - begin.tv_sec) + ((end.tv_usec - begin.tv_usec)/1000000.0));
        // printf("Efield hfield done\n");
        t = t + dt;

		KRMS=1;
        //!=============================================================
        if(fmod(n,nperdt)==0)         //! when n is integral multiples of period
        	KRMS=2;

        gettimeofday(&begin, NULL);
        RMS(KRMS);                     //! finally the RMS calculated as KRMS=2 on each period
        
        // cudaMalloc((void**)&dev_root_elec, sizeof(root_elec));
        // cudaMalloc((void**)&dev_z1, sizeof(double));
        // cudaMalloc((void**)&dev_z2, sizeof(double));
        // cudaMalloc((void**)&dev_inv_nperdt, sizeof(inv_nperdt));
        // cudaMalloc((void**)&dev_ext, sizeof(ext));
        // cudaMalloc((void**)&dev_eyt, sizeof(eyt));
        // cudaMalloc((void**)&dev_ERMSp, sizeof(ERMSp));
        // cudaMalloc((void**)&dev_erms2, sizeof(erms2));
        // cudaMalloc((void**)&dev_KRMS, sizeof(KRMS));

        // cudaMemcpy(dev_root_elec, root_elec, sizeof(root_elec), cudaMemcpyHostToDevice);
        // cudaMemcpy(dev_inv_nperdt, &inv_nperdt, sizeof(inv_nperdt), cudaMemcpyHostToDevice);
        // cudaMemcpy(dev_ext, ext, sizeof(ext), cudaMemcpyHostToDevice);
        // cudaMemcpy(dev_eyt, eyt, sizeof(eyt), cudaMemcpyHostToDevice);
        // cudaMemcpy(dev_ERMSp, ERMSp, sizeof(ERMSp), cudaMemcpyHostToDevice);
        // cudaMemcpy(dev_erms2, erms2, sizeof(erms2), cudaMemcpyHostToDevice);
        // cudaMemcpy(dev_KRMS, &KRMS, sizeof(int), cudaMemcpyHostToDevice);

        // RMS<<<(ceil(root_elec->m/32),ceil(root_elec->n/32)),(32,32)>>>(dev_root_elec, dev_z1,dev_z2,dev_inv_nperdt,dev_ext,dev_eyt,dev_ERMSp,dev_erms2,dev_KRMS);
        
        // cudaMemcpy(ERMSp, dev_ERMSp, sizeof(dev_ERMSp), cudaMemcpyDeviceToHost);
        // cudaMemcpy(erms2, dev_erms2, sizeof(dev_erms2), cudaMemcpyDeviceToHost);
        // if(KRMS==2)
        // {
        //     cudaMemcpy(root_elec->mesh, dev_root_elec->mesh, sizeof(dev_root_elec->mesh), cudaMemcpyDeviceToHost);
        // }
       
        // cudaFree(dev_root_elec);
        // //cudaFree(dev_z1);
        // //cudaFree(dev_z2);
        // cudaFree(dev_inv_nperdt);
        // cudaFree(dev_ext);
        // cudaFree(dev_eyt);
        // cudaFree(dev_ERMSp);
        // cudaFree(dev_erms2);
        // cudaFree(dev_KRMS);

        gettimeofday(&end, NULL);
        t_rms += ((end.tv_sec - begin.tv_sec) + ((end.tv_usec - begin.tv_usec)/1000000.0));
        //!=============================================================
	     /* 
       if(n>=1)
       {
            ani++;
            denx();
            deny();
            diffux();
            diffuy();
            ionfreqx();
            ionfreqy();
            Ermsx();
            Ermsy();
       
       }*/
       
        //!%%%%%%% print the temporal profile of Ey-field (Total)  at i=0.75*nx,0.55*nx,0.35*nx,0.15*nx, j=0.5*ny for all time instants %%%%%%%%%%
       //Eytempo();
       //!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       //!%%%%%%% print the temporal profile of density at j=0.74*ny,0.26*ny, i=(0.38)*nx for all time instants %%%%%%%%%%
      // denytempo();
       //!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       
	     
        if(fmod(n,nperdt*nmaxwell)==0)
        {
            if(icpling!=0)
            { 
                gettimeofday(&begin, NULL);
                ELEC_DENS();                  //! Electron density is called on start of each period under RMS field
                gettimeofday(&telcend, NULL);
                t_elec_dens += ((end.tv_sec - begin.tv_sec) + ((end.tv_usec - begin.tv_sec)/1000000.0));
            }
            else
            {
                DTAC=1.0/FREQ;
                TIMD=(double)(n)*inv_nperdt*nmaxwell*DTAC;
            }
            //!=============================================================XXX
            KELEC=KELEC+1;
           if(n==1){
           //parentgrid();
             //   childgrid();   
           }
          
            
                 gettimeofday(&begin, NULL);

   		 //if(KELEC==25||KELEC==50||KELEC==100||KELEC==200||KELEC==250||KELEC==500||KELEC==1000||KELEC==1200||KELEC==1500||KELEC==1600||KELEC==1825||KELEC==1850||KELEC==1925||KELEC==1950||KELEC==2000){
       // if(KELEC==25||KELEC==50||KELEC==100||KELEC==200||KELEC==500||KELEC==1000||KELEC==2000||KELEC==3000||KELEC==4000||KELEC==4200||KELEC==4400||KELEC==4800||KELEC==5000||KELEC==5200||KELEC==5400){
       //if(KELEC==25||KELEC==50||KELEC==100||KELEC==200||KELEC==500||KELEC==1000||KELEC==2000||KELEC==3000){
 //      if(KELEC==1||KELEC==50||KELEC==500||KELEC==1000||KELEC==2000||KELEC==3000||KELEC==4000||KELEC==4200||KELEC==4400||KELEC==4800||KELEC==5000||KELEC==5200||KELEC==5400||KELEC==6000||KELEC==6500||KELEC==7000||KELEC==7500||KELEC==8000||KELEC==8500||KELEC==9000||KELEC==9500||KELEC==10500||KELEC==11100){
      //if(KELEC==1||KELEC==50||KELEC==100||KELEC==200||KELEC==352||KELEC==400||KELEC==500||KELEC==650||KELEC==700||KELEC==800||KELEC==850||KELEC==900||KELEC==950||KELEC==1000||KELEC==1200||KELEC==1600||KELEC==1800||KELEC==2000||KELEC==2200||KELEC==2500||KELEC==2800||KELEC==3000){
      //if(KELEC==1||KELEC==50||KELEC==100||KELEC==200||KELEC==352||KELEC==400||KELEC==500||KELEC==650||KELEC==700||KELEC==800||KELEC==850||KELEC==900||KELEC==950||KELEC==1000||KELEC==1200||KELEC==1600||KELEC==1800||KELEC==2000||KELEC==2200||KELEC==2500||KELEC==2800||KELEC==3000||KELEC==3400||KELEC==3600||KELEC==3800||KELEC==4000||KELEC==4200||KELEC==4400||KELEC==4600||KELEC==4800||KELEC==5000||KELEC==5200||KELEC==5400){
      // if(KELEC==1||KELEC==50||KELEC==900||KELEC==905||KELEC==910||KELEC==915||KELEC==920||KELEC==925||KELEC==930||KELEC==935||KELEC==940||KELEC==945||KELEC==950){
      if(KELEC==1||KELEC==50||KELEC==100||KELEC==200||KELEC==352||KELEC==400||KELEC==500||KELEC==650||KELEC==700||KELEC==800||KELEC==850||KELEC==900||KELEC==950||KELEC==1000||KELEC==1200||KELEC==1600||KELEC==1800||KELEC==2000||KELEC==2200||KELEC==2500||KELEC==2800||KELEC==3000||KELEC==3400||KELEC==3600||KELEC==3800||KELEC==4000||KELEC==4200||KELEC==4400||KELEC==4600||KELEC==4800||KELEC==5000||KELEC==5200||KELEC==5400||KELEC==6000||KELEC==6500||KELEC==7000||KELEC==7500||KELEC==8000||KELEC==8500||KELEC==9000||KELEC==9500||KELEC==10500||KELEC==11100){
      		  // ani++;    
           ansetup2++;   
      		 printf("%d %lf\n",ansetup2,t*1.0e9);
           panim();
           gettimeofday(&t_panimend, NULL); 
           tkl=t_panimend.tv_sec; //!new
           info = localtime(&tkl);//!new
           
            printf("%s",asctime (info));//!new
            
           canitimxy=fopen("canitimxy.txt","a"); //! writes the runtime and simulation time correspond to each KELEC, to detect dynamic mesh: dynamic runtime 
           fprintf(canitimxy,"%d ",KELEC);
           fprintf(canitimxy,"%lf ",t*1.0e9);     //! simulation time
           t_panimcal = ((t_panimend.tv_sec - total_start.tv_sec) + ((t_panimend.tv_usec - total_start.tv_usec)/1000000.0));  
           fprintf(canitimxy,"%lf \n",t_panimcal);  //!runtime
           fclose(canitimxy); 
            
           // anim();
           // animE();
	          //bothdenEx();
      	   // bothdenEy();
            
             //  canim();
		      //canimE();  
              }
            //ani++;
           // printf("%d KELEC=%d t=%lf time=%lf\n",ani,KELEC,t,t*1.0e12);
          /*
            denx();
	          diffux();
	          ionfreqx();	
             Ermsx();   
           */     
           /*
            deny();
            diffux();
            diffuy();
            ionfreqx();
            ionfreqy();
            Ermsx();
            Ermsy();
          */
            
          //!%%%%%%% print the temporal profile of Ey-field (Total)  %%%%%%%%%%
          //Eytempo();
         //!%%%%%%% print the temporal profile of density %%%%%%%%%%
         // denytempo();
         //!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         // denytempocent();
          //!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         // Eytempocent();
          //!%%%%%%%%%%%%%%%%%%%%%
                gettimeofday(&end, NULL);
                t_anim += ((end.tv_sec - begin.tv_sec) + ((end.tv_usec - begin.tv_usec)/1000000.0));
	  
	     //!=============================================================XXXXX
           
            if(KELEC==1)
                fptr=fopen("RES.out","w");
            if(KELEC>1)
            { 
                fptr=fopen("RES.out","a");
                if(fptr!=NULL)
                    fprintf(fptr,"%d %f %f %f %f",n,TIMD,DTAC,ACCEL,dnma);
                else
                    goto l100;
                fclose(fptr);
            }
        }

        l100:

        //!=============================================================XXXXX

        nmod=2*nperdt;

        t=t-dt; //! going back in t=0, for completing the factor time child updates
         for(lup=1;lup<factor;lup++)
         {
            c_n++;
		
            child_HFIELD();
            
            // cudaMalloc((void**)&dev_mag, sizeof(child_mag));
            // cudaMalloc((void**)&dev_c_exs, sizeof(c_exs));
            // cudaMalloc((void**)&dev_c_eys, sizeof(c_eys));
            // cudaMalloc((void**)&dev_dtmds, sizeof(dtmds));
            // cudaMemcpy(dev_mag, child_mag, sizeof(child_mag), cudaMemcpyHostToDevice);
            // cudaMemcpy(dev_exs, c_exs, sizeof(c_exs), cudaMemcpyHostToDevice);
            // cudaMemcpy(dev_eys, c_eys, sizeof(c_eys), cudaMemcpyHostToDevice);
            // cudaMemcpy(dev_dtmds, &dtmds, sizeof(dtmds), cudaMemcpyHostToDevice);
            // HFIELD<<<(ceil(child_mag->m/32),ceil(child_mag->n/32)),(32,32)>>>(dev_mag,dev_exs,dev_eys,dev_dtmds);
            // cudaMemcpy(child_mag, dev_mag, sizeof(child_mag), cudaMemcpyDeviceToHost);
            // cudaMemcpy(c_hzi, dev_mag->mesh, sizeof(child_mag), cudaMemcpyDeviceToHost);

            // cudaFree(dev_mag);
            // cudaFree(dev_c_exs);
            // cudaFree(dev_c_eys);
            // cudaFree(dev_dtmds);
        child_EFIELD(); 
        // cudaMalloc((void**)&dev_elec, sizeof(child_elec));
        // cudaMalloc((void**)&dev_x0, sizeof(double));
        // cudaMalloc((void**)&dev_OMEG, sizeof(double));
        // cudaMalloc((void**)&dev_newt, sizeof(double));
        // cudaMalloc((void**)&dev_inv_c, sizeof(double));
        // cudaMalloc((void**)&dev_c_dt, sizeof(double));
        // cudaMalloc((void**)&dev_sine, sizeof(double));
        // cudaMalloc((void**)&dev_sine1, sizeof(double));
        // cudaMalloc((void**)&dev_x, sizeof(double));
        // cudaMalloc((void**)&dev_c, sizeof(double));
        // cudaMalloc((void**)&dev_c_ds, sizeof(double));

        // cudaMemcpy(dev_elec, child_elec, sizeof(child_elec), cudaMemcpyHostToDevice);
        // cudaMemcpy(dev_x0, x0, sizeof(double), cudaMemcpyHostToDevice);
        // cudaMemcpy(dev_OMEG, OMEG, sizeof(double), cudaMemcpyHostToDevice);
        // cudaMemcpy(dev_newt, newt, sizeof(double), cudaMemcpyHostToDevice);
        // cudaMemcpy(dev_inv_c, inv_c, sizeof(double), cudaMemcpyHostToDevice);
        // cudaMemcpy(dev_c_dt, c_dt, sizeof(double), cudaMemcpyHostToDevice);
        // cudaMemcpy(dev_c, c, sizeof(double), cudaMemcpyHostToDevice);
        // cudaMemcpy(dev_c_ds, c_ds, sizeof(double), cudaMemcpyHostToDevice);
        
        
        // child_EFIELD<<<dimGrid, dimBlock>>>(dev_elec,dev_x0,dev_OMEG,dev_newt,dev_inv_c,dev_c_dt,dev_sine,dev_sine1,dev_x,dev_c,dev_c_ds);
        
        // cudaMemcpy(sine, dev_sine, sizeof(dev_sine), cudaMemcpyDeviceToHost);
        // cudaMemcpy(sine1, dev_sine1, sizeof(dev_sine1), cudaMemcpyDeviceToHost);

        // cudaFree(dev_elec);
        // cudaFree(dev_x0);
        // cudaFree(dev_OMEG);
        // cudaFree(dev_newt);
        // cudaFree(dev_inv_c);
        // cudaFree(dev_c_dt);
        // cudaFree(dev_sine);
        // cudaFree(dev_sine1);
        // cudaFree(dev_x);
        // cudaFree(dev_c);
        // cudaFree(dev_c_ds);

            //! interpolate corners, along x along y for E-,H-,electron density for child mesh from parent
            //! inter_value decides the fraction of data contributed from before and updated parent data (provided factor-1 times child update)
        interpolatecorners(inter_value);      
        interpolatex(inter_value);
        interpolatey(inter_value); 
                        
        newt = t + c_dt;            
        KRMS=1;

        child_RMS(KRMS);

        // cudaMalloc((void**)&dev_child_elec, sizeof(child_elec));
        // cudaMalloc((void**)&dev_z1, sizeof(double));
        // cudaMalloc((void**)&dev_z2, sizeof(double));
        // cudaMalloc((void**)&dev_inv_nperdt, sizeof(c_inv_nperdt));
        // cudaMalloc((void**)&dev_ext, sizeof(c_ext));
        // cudaMalloc((void**)&dev_eyt, sizeof(c_eyt));
        // cudaMalloc((void**)&dev_erms2, sizeof(c_erms2));
        // cudaMalloc((void**)&dev_KRMS, sizeof(KRMS));

        // cudaMemcpy(dev_child_elec, child_elec, sizeof(child_elec), cudaMemcpyHostToDevice);
        // cudaMemcpy(dev_inv_nperdt, &c_inv_nperdt, sizeof(c_inv_nperdt), cudaMemcpyHostToDevice);
        // cudaMemcpy(dev_ext, c_ext, sizeof(c_ext), cudaMemcpyHostToDevice);
        // cudaMemcpy(dev_eyt, c_eyt, sizeof(c_eyt), cudaMemcpyHostToDevice);
        // //cudaMemcpy(dev_ERMSp, c_ERMSp, sizeof(c_ERMSp), cudaMemcpyHostToDevice);
        // cudaMemcpy(dev_erms2, c_erms2, sizeof(c_erms2), cudaMemcpyHostToDevice);
        // cudaMemcpy(dev_KRMS, &KRMS, sizeof(KRMS), cudaMemcpyHostToDevice);

        // child_RMS<<<(ceil(child_elec->m/32),ceil(child_elec->n/32)),(32,32)>>>(dev_child_elec, dev_z1,dev_z2,dev_inv_nperdt,dev_ext,dev_eyt,dev_erms2,dev_KRMS);
        
        // cudaMemcpy(c_erms2, dev_erms2, sizeof(dev_erms2), cudaMemcpyDeviceToHost);
        // if(KRMS==2)
        // {
        //     cudaMemcpy(child_elec->mesh, dev_child_elec->mesh, sizeof(dev_child_elec->mesh), cudaMemcpyDeviceToHost);
        // }
       
        // cudaFree(dev_child_elec);
        // //cudaFree(dev_z1);
        // //cudaFree(dev_z2);
        // cudaFree(dev_inv_nperdt);
        // cudaFree(dev_ext);
        // cudaFree(dev_eyt);
        // cudaFree(dev_erms2);
        // cudaFree(dev_KRMS);

        }

	//!-----------------------(final)
	
  c_n++;

        child_HFIELD();
        // cudaMalloc((void**)&dev_mag, sizeof(child_mag));
        // cudaMalloc((void**)&dev_c_exs, sizeof(c_exs));
        // cudaMalloc((void**)&dev_c_eys, sizeof(c_eys));
        // cudaMalloc((void**)&dev_dtmds, sizeof(dtmds));
        // cudaMemcpy(dev_mag, child_mag, sizeof(child_mag), cudaMemcpyHostToDevice);
        // cudaMemcpy(dev_c_exs, c_exs, sizeof(c_exs), cudaMemcpyHostToDevice);
        // cudaMemcpy(dev_c_eys, c_eys, sizeof(c_eys), cudaMemcpyHostToDevice);
        // cudaMemcpy(dev_dtmds, &dtmds, sizeof(dtmds), cudaMemcpyHostToDevice);
        // HFIELD<<<(ceil(child_mag->m/32),ceil(child_mag->n/32)),(32,32)>>>(dev_mag,dev_c_exs,dev_c_eys,dev_dtmds);
        // cudaMemcpy(child_mag, dev_mag, sizeof(child_mag), cudaMemcpyDeviceToHost);
        // cudaMemcpy(c_hzi, dev_mag->mesh, sizeof(child_mag), cudaMemcpyDeviceToHost);

        // cudaFree(dev_mag);
        // cudaFree(dev_c_exs);
        // cudaFree(dev_c_eys);
        // cudaFree(dev_dtmds);
        child_EFIELD();
        // cudaMalloc((void**)&dev_elec, sizeof(child_elec));
        // cudaMalloc((void**)&dev_x0, sizeof(double));
        // cudaMalloc((void**)&dev_OMEG, sizeof(double));
        // cudaMalloc((void**)&dev_newt, sizeof(double));
        // cudaMalloc((void**)&dev_inv_c, sizeof(double));
        // cudaMalloc((void**)&dev_c_dt, sizeof(double));
        // cudaMalloc((void**)&dev_sine, sizeof(double));
        // cudaMalloc((void**)&dev_sine1, sizeof(double));
        // cudaMalloc((void**)&dev_x, sizeof(double));
        // cudaMalloc((void**)&dev_c, sizeof(double));
        // cudaMalloc((void**)&dev_c_ds, sizeof(double));

        // cudaMemcpy(dev_elec, child_elec, sizeof(child_elec), cudaMemcpyHostToDevice);
        // cudaMemcpy(dev_x0, x0, sizeof(double), cudaMemcpyHostToDevice);
        // cudaMemcpy(dev_OMEG, OMEG, sizeof(double), cudaMemcpyHostToDevice);
        // cudaMemcpy(dev_newt, newt, sizeof(double), cudaMemcpyHostToDevice);
        // cudaMemcpy(dev_inv_c, inv_c, sizeof(double), cudaMemcpyHostToDevice);
        // cudaMemcpy(dev_c_dt, c_dt, sizeof(double), cudaMemcpyHostToDevice);
        // cudaMemcpy(dev_c, c, sizeof(double), cudaMemcpyHostToDevice);
        // cudaMemcpy(dev_c_ds, c_ds, sizeof(double), cudaMemcpyHostToDevice);
        
        
        // child_EFIELD<<<dimGrid, dimBlock>>>(dev_elec,dev_x0,dev_OMEG,dev_newt,dev_inv_c,dev_c_dt,dev_sine,dev_sine1,dev_x,dev_c,dev_c_ds);
        
        // cudaMemcpy(sine, dev_sine, sizeof(dev_sine), cudaMemcpyDeviceToHost);
        // cudaMemcpy(sine1, dev_sine1, sizeof(dev_sine1), cudaMemcpyDeviceToHost);  

        // cudaFree(dev_elec);
        // cudaFree(dev_x0);
        // cudaFree(dev_OMEG);
        // cudaFree(dev_newt);
        // cudaFree(dev_inv_c);
        // cudaFree(dev_c_dt);
        // cudaFree(dev_sine);
        // cudaFree(dev_sine1);
        // cudaFree(dev_x);
        // cudaFree(dev_c);
        // cudaFree(dev_c_ds);
        //! interpolate corners, along x along y for E-,H-,electron density for child mesh from parent
        //! inter_value decides the fraction of data contributed only from updated parent data (provided on the factor th update of child)   
        interpolatecorners(1.0);
        interpolatex(1.0);
        interpolatey(1.0);
       	
       
        t=t+dt;
        newt=t;
        
        KRMS=1;
        if(fmod(c_n,c_nperdt)==0) {
                   
            KRMS=2;
        }
        child_RMS(KRMS);

      //   cudaMalloc((void**)&dev_child_elec, sizeof(child_elec));
      //   cudaMalloc((void**)&dev_z1, sizeof(double));
      //   cudaMalloc((void**)&dev_z2, sizeof(double));
      //   cudaMalloc((void**)&dev_inv_nperdt, sizeof(c_inv_nperdt));
      //   cudaMalloc((void**)&dev_ext, sizeof(c_ext));
      //   cudaMalloc((void**)&dev_eyt, sizeof(c_eyt));
      //   cudaMalloc((void**)&dev_erms2, sizeof(c_erms2));
      //   cudaMalloc((void**)&dev_KRMS, sizeof(KRMS));

      //   cudaMemcpy(dev_child_elec, child_elec, sizeof(child_elec), cudaMemcpyHostToDevice);
      //   cudaMemcpy(dev_inv_nperdt, &c_inv_nperdt, sizeof(c_inv_nperdt), cudaMemcpyHostToDevice);
      //   cudaMemcpy(dev_ext, c_ext, sizeof(c_ext), cudaMemcpyHostToDevice);
      //   cudaMemcpy(dev_eyt, c_eyt, sizeof(c_eyt), cudaMemcpyHostToDevice);
      //   //cudaMemcpy(dev_ERMSp, c_ERMSp, sizeof(c_ERMSp), cudaMemcpyHostToDevice);
      //   cudaMemcpy(dev_erms2, c_erms2, sizeof(c_erms2), cudaMemcpyHostToDevice);
      //   cudaMemcpy(dev_KRMS, &KRMS, sizeof(KRMS), cudaMemcpyHostToDevice);

      //   child_RMS<<<(ceil(child_elec->m/32),ceil(child_elec->n/32)),(32,32)>>>(dev_child_elec, dev_z1,dev_z2,dev_inv_nperdt,dev_ext,dev_eyt,dev_erms2,dev_KRMS);
        
      //   cudaMemcpy(c_erms2, dev_erms2, sizeof(dev_erms2), cudaMemcpyDeviceToHost);
      //   if(KRMS==2)
      //   {
      //       cudaMemcpy(child_elec->mesh, dev_child_elec->mesh, sizeof(dev_child_elec->mesh), cudaMemcpyDeviceToHost);
      //   }
       
      //   cudaFree(dev_child_elec);
      //  // cudaFree(dev_z1);
      //   //cudaFree(dev_z2);
      //   cudaFree(dev_inv_nperdt);
      //   cudaFree(dev_ext);
      //   cudaFree(dev_eyt);
      //   cudaFree(dev_erms2);
      //   cudaFree(dev_KRMS);
        

        if(fmod(c_n,c_nperdt*nmaxwell)==0){
            if(icpling!=0){
             child_ELEC_DENS();     //!  each child period update the child electron density

            }
            else{
                DTAC=1.0/FREQ;
                TIMD=(double)(c_n)*inv_nperdt*nmaxwell*DTAC/factor;
            }
            c_KELEC++;
        }

	        c2p();         //! only after the child data updated between a parent update copy the fine (child) data back to parent
                          //! to avoid mismatch in accuracy . copy cpb(child-parent) boundary 
                 
                    
    }while(1);

        //! prints the calculation time for the each of the subroutines (remove the rank as it is for MPI code)
        fptr2 = fopen("mes_result.csv","a");
        fprintf(fptr2,"Processor.Rank ");            
        fprintf(fptr2,"Cal.EFIELD ");            
        fprintf(fptr2,"Cal.Child_EFIELD ");            
        fprintf(fptr2,"Cal.HFIELD ");           
        fprintf(fptr2,"Cal.Child_HFIELD ");        
        fprintf(fptr2,"Cal.ELEC_DENS " );        
        fprintf(fptr2,"Cal.CHILD_ELEC_DENS ");        
        fprintf(fptr2,"Cal.RMS ");           
        fprintf(fptr2,"Cal.CHILD_RMS ");
        fprintf(fptr2,"Cal.c2p ");
        fprintf(fptr2,"Cal.interpolatex ");
        fprintf(fptr2,"Cal.interpolatexinitial ");
        fprintf(fptr2,"Cal.interpolatexnew ");
        fprintf(fptr2,"Cal.interpolateynew ");
        fprintf(fptr2,"Cal.anim ");
        fprintf(fptr2,"Total.time ");
        fprintf(fptr2,"Program.time\n");                
        fclose(fptr2);
        
       fptr2 = fopen("mes_result.csv","a");
        fprintf(fptr2,"%d ",rank);            
        fprintf(fptr2,"%f ",t_cal_efield);            
        fprintf(fptr2,"%f ",t_cal_child_efield);            
        fprintf(fptr2,"%f ",t_cal_hfield);           
        fprintf(fptr2,"%f ",t_cal_child_hfield);           
        fprintf(fptr2,"%f ",t_cal_elec_dens);            
        fprintf(fptr2,"%f ",t_cal_child_elec_dens);            
        fprintf(fptr2,"%f ",t_cal_rms);           
        fprintf(fptr2,"%f ",t_cal_child_rms);           
        fprintf(fptr2,"%f ",t_cal_c2p);           
        fprintf(fptr2,"%f ",t_cal_interpolatex); 
        fprintf(fptr2,"%f ",t_cal_interpolatexinitial);
        fprintf(fptr2,"%f ",t_cal_interpolatexnew);
        fprintf(fptr2,"%f ",t_cal_interpolateynew);
        fprintf(fptr2,"%f ",t_cal_anim);            
        fprintf(fptr2,"%f ",t1);                       
        gettimeofday(&program_end,NULL);    
        t1 = ((program_end.tv_sec - program_start.tv_sec) + ((program_end.tv_usec - program_start.tv_usec)/1000000.0));
        fprintf(fptr2,"%f\n",t1);                       
        fclose(fptr2);      

      printf("%e %d %d\n", dt,nx,ny); 
      printf("OUTPUT 2\n");
      gettimeofday(&total_end, NULL);
      printf("EFIELD time: %f s\n", t_efield_hfield);
      printf("ELEC_DENS time: %f s\n", t_elec_dens);
      printf("RMS time: %f s\n", t_rms);
      printf("Zero time: %f s\n", t_zero);
      printf("Anim time: %f s\n", t_anim);
      printf("Total time: %f s\n", ((total_end.tv_sec - total_start.tv_sec) + ((total_end.tv_usec - total_start.tv_usec)/1000000.0)));
      
  	free_all();
    return 0;
}

void EFIELD()
{    
     gettimeofday(&begin,NULL);
    for(i=0; i<root_elec->m;i++)
    {
    	x=i*ds;          
        sine=0.0;
        sine1=0.0;
        if(x<=(x0+c*t))
        { 
            sine = sin(OMEG*(t-(x-x0)*inv_c));  //!forward wave
        }
        if(x<=(x0+c*(t+dt))) 
        {   
            sine1 = sin(OMEG*(t+dt-(x-x0)*inv_c));  //!forward wave
        }

        

        for(j=0;j<root_elec->n;j++)
        {
            
            eyi[i][j] =   E0*(sine);  //! incident field () at t
            eyi1[i][j] =  E0*(sine1); //! incident field () at t+dt
         
        }

        for(j=0;j<root_elec->n;j++)
        {
            omp2x=(root_den->mesh[i][j]+root_den->mesh[i+1][j])*const3;
            betax=omp2x*const4;
            const5x=1.0/(1.0+betax);
            const6x=1.0-betax;
            if(j>0)
            {
            	extk=ext[i][j];
                exs_old[i][j] = exs[i][j];
             //! Scattered E-field x update depends on previous E-data, previous H-data (top/bottom), density data (space:left/ right),velocity (same location)     
                exs[i][j]=const5x*( exs[i][j]*(const6x)+qe*(root_den->mesh[i][j]+root_den->mesh[i+1][j])*vx[i][j]*const7
                                     -(exi[i][j]+exi1[i][j])*betax+(root_mag->mesh[i][j]-root_mag->mesh[i][j-1])*dteds); 

                vx_old[i][j] = vx[i][j];
				ext[i][j]=exs[i][j]+exi1[i][j];     //! Total field = inci + scattered 
                vx[i][j]=vx[i][j]*alpha - qmdt*(ext[i][j]+extk)*const8;
        
               
            }
            if(i>0)
            {

                omp2y=(root_den->mesh[i][j]+root_den->mesh[i][j+1])*const3;
                betay=omp2y*const4;
                const5y=1.0/(1.0+betay);
                const6y=1.0-betay;

                eytk=eyt[i][j];
                eys_old[i][j] = eys[i][j];      //! reqd for interpolation
                //!Scattered E-field y update depends on previous E-data, previous H-data (left/right), density data (space:top/ bottom),velocity (same location)
                eys[i][j]=const5y*(eys[i][j]*(const6y)+qe*(root_den->mesh[i][j]+root_den->mesh[i][j+1])*vy[i][j]*const7
                                 -(eyi[i][j]+eyi1[i][j])*betay-(root_mag->mesh[i][j]-root_mag->mesh[i-1][j])*dteds);


                eyt[i][j]=eys[i][j]+eyi1[i][j];    //! Total field = inci + scattered 
                vy_old[i][j] = vy[i][j];        //! reqd for interpolation
                vy[i][j]=vy[i][j]*alpha - qmdt*(eyt[i][j]+eytk)*const8;     //! velocity update
            }
        }
        MR_MUR(i);
    }
    gettimeofday(&end,NULL);    
    t_cal_efield += ((end.tv_sec - begin.tv_sec) + ((end.tv_usec - begin.tv_usec)/1000000.0));
}

  //int dev_child_elec,dev_x0,dev_OMEG,dev_newt,dev_inv_c,dev_c_dt,dev_sine,dev_sine1,dev_x 

// __global__ EFIELD(struct node * elec,double x0,double OMEG,double newt,double inv_c,double c_dt,double sine,double sine1,double x,double c)
// {
//     int i = blockIdx.x*blockDim.x+threadIdx.x;
//     if(i<elec->m)
//     {
//         x=i*ds;          
//          sine=0.0;
//          sine1=0.0;
//          if(x<=(x0+c*t))
//          { 
//              sine = sin(OMEG*(t-(x-x0)*inv_c));  //!forward wave
//          }
//          if(x<=(x0+c*(t+dt))) 
//          {   
//              sine1 = sin(OMEG*(t+dt-(x-x0)*inv_c));  //!forward wave
//          }
//     }

// }

// __global__ child_EFIELD(struct node * child_elec,double x0,double OMEG,double newt,double inv_c,double c_dt,double sine,double sine1,double x,double c,double c_ds)
// {
//     int i = blockIdx.x*blockDim.x+threadIdx.x;
//     if(i<elec->m)
//     {
//         x=(i+factor*child_elec->locx)*c_ds;          
//         sine=0.0;
//         sine1=0.0;
// 	        if(x<=(x0+c*(newt)))
// 	        { 
// 	            sine = sin(OMEG*(newt-(x-x0)*inv_c));    //!forward wave
	           
// 	        }
// 	        if(x<=(x0+c*(newt+c_dt))) 
// 	        {   
	        	
// 	            sine1 = sin(OMEG*(newt+c_dt-(x-x0)*inv_c));    //!forward wave
	            
// 	        }
	    
//     }
// }
void child_EFIELD()
{
	
   gettimeofday(&begin,NULL);
	for(i=0; i<child_elec->m;i++)
    {
    	x=(i+factor*child_elec->locx)*c_ds;          
        sine=0.0;
        sine1=0.0;
	        if(x<=(x0+c*(newt)))
	        { 
	            sine = sin(OMEG*(newt-(x-x0)*inv_c));    //!forward wave
	           
	        }
	        if(x<=(x0+c*(newt+c_dt))) 
	        {   
	        	
	            sine1 = sin(OMEG*(newt+c_dt-(x-x0)*inv_c));    //!forward wave
	            
	        }
	    
    	    
          for(j=0;j<child_elec->n;j++)
    	    {
      	    
              c_eyi[i][j] =   E0*(sine);    //! incident field ( forward direction) at newt
    	        
    	        c_eyi1[i][j] =  E0*(sine1);    //! incident field ( forward direction) at newt+c_dt
    	    }

        for(j=0;j<child_elec->n;j++)
        {
            omp2x=(child_den->mesh[i][j]+child_den->mesh[i+1][j])*const3;
            betax=omp2x*c_const4;
            const5x=1.0/(1.0+betax);
            const6x=1.0-betax;
            if(j>0)
            {
                extk=c_ext[i][j];
                c_exs[i][j]=const5x*( c_exs[i][j]*(const6x)+qe*(child_den->mesh[i][j]+child_den->mesh[i+1][j])*c_vx[i][j]*c_const7 
                                     -(c_exi[i][j]+c_exi1[i][j])*betax+(child_mag->mesh[i][j]-child_mag->mesh[i][j-1])*c_dteds); //const7 is divided by 2 because of the dt factor.

                c_exold[i][j]=c_exs[i][j];    //! reqd when expanding the child mesh to copy the previous mesh data
                
                c_ext[i][j]=c_exs[i][j]+c_exi1[i][j];   //! Total field (child)= inci + scattered 


                c_vx[i][j]=c_vx[i][j]*c_alpha - c_qmdt*(c_ext[i][j]+extk)*c_const8;
               
                c_vxold[i][j]=c_vx[i][j];     //! reqd when expanding the child mesh to copy the previous mesh data
            }
            if(i>0)
            {

                omp2y=(child_den->mesh[i][j]+child_den->mesh[i][j+1])*const3;
                betay=omp2y*c_const4;
                const5y=1.0/(1.0+betay);
                const6y=1.0-betay;

                eytk=c_eyt[i][j];

                double v = c_eys[i][j]*(const6y);

                c_eys[i][j]=const5y*(c_eys[i][j]*(const6y)+qe*(child_den->mesh[i][j]+child_den->mesh[i][j+1])*c_vy[i][j]*c_const7
                                 -(c_eyi[i][j]+c_eyi1[i][j])*betay-(child_mag->mesh[i][j]-child_mag->mesh[i-1][j])*c_dteds);
               	            
                c_eyold[i][j]=c_eys[i][j];  //! reqd when expanding the child mesh to copy the previous mesh data

                c_eyt[i][j]=c_eys[i][j]+c_eyi1[i][j];  //! Total field (child)= inci + scattered 

                c_vy[i][j]=c_vy[i][j]*c_alpha - c_qmdt*(c_eyt[i][j]+eytk)*c_const8;
                
                c_vyold[i][j]=c_vy[i][j];   //! reqd when expanding the child mesh to copy the previous mesh data
            }
        }
    }
    gettimeofday(&end,NULL);    
    t_cal_child_efield += ((end.tv_sec - begin.tv_sec) + ((end.tv_usec - begin.tv_usec)/1000000.0));
}

 void HFIELD()
 {
     gettimeofday(&begin,NULL); 
     for(i=0; i<root_mag->m; i++)
     {
     	for(j=0; j<root_mag->n; j++)
     	{
     		hzi[i][j] = root_mag->mesh[i][j]; //! if we give H-excitation rather E-excitation 
         //! H-field (parent) depends on E-field x:  top/bottom and E-field y: left/right previous (before E-field updates then H-field)
    		root_mag->mesh[i][j]+= (-(eys[i+1][j]-eys[i][j])+(exs[i][j+1]-exs[i][j]))*dtmds;
     	}
     }
     gettimeofday(&end,NULL);    
     t_cal_hfield += ((end.tv_sec - begin.tv_sec) + ((end.tv_usec - begin.tv_usec)/1000000.0));
 }   


 void child_HFIELD()
 {
    gettimeofday(&begin,NULL);
 	for(i=0; i<child_mag->m; i++)
     {
     	for(j=0; j<child_mag->n; j++)
     	{
     		c_hzi[i][j] = child_mag->mesh[i][j];
          //! H-field (parent) depends on E-field x:  top/bottom and E-field y: left/right previous (before E-field updates then H-field)
     		child_mag->mesh[i][j]+= (-(c_eys[i+1][j]-c_eys[i][j])+(c_exs[i][j+1]-c_exs[i][j]))*c_dtmds;
        
         c_hzold[i][j]=child_mag->mesh[i][j];  
     	}
     }
     gettimeofday(&end,NULL);    
     t_cal_child_hfield += ((end.tv_sec - begin.tv_sec) + ((end.tv_usec - begin.tv_usec)/1000000.0));
 }


void RMS(int k)
{
    //printf("In RMS\n");
    gettimeofday(&begin,NULL);
    double z1,z2;
    
    int i;
    for(i=1;i<root_elec->m;i++)
    {
        for(j=0;j<root_elec->n;j++)
        {
            z1=(ext[i][j]*ext[i][j]+ext[i-1][j]*ext[i-1][j])*.5f;   //! avg of the two scattered field is required (E_eff) for the density update
            z2=(eyt[i][j]*eyt[i][j]+eyt[i][j-1]*eyt[i][j-1])*.5f;   //! avg of the two scattered field is required (E_eff) for the density update
            ERMSp[i][j] = erms2[i][j];
            erms2[i][j]=erms2[i][j]+(z1+z2)*inv_nperdt;     //! time updates and averages (parent)
        }
    }

    if(k==2)
    {
        int j;
        for(i=0;i<root_elec->m;i++)
        {
            for(j=0;j<root_elec->n;j++)
            {
            	if (erms2[i][j]<0)
            	{
            		printf("Alert!!\n");
            	}
                root_elec->mesh[i][j] = sqrt(erms2[i][j]);  //! completes a period and then squre root the time avg data (parent)
                erms2[i][j]=0.0f;
            } 
        }
    }
    gettimeofday(&end,NULL);
    t_cal_rms += ((end.tv_sec - begin.tv_sec) + ((end.tv_usec - begin.tv_usec)/1000000.0));
}

void child_RMS(int k)
{
	// printf("In RMS\n");
   gettimeofday(&begin,NULL); 
	double z1,z2;

    int i;
    for(i=1;i<child_elec->m;i++)
    {
        for(j=1;j<child_elec->n;j++)
        {
            z1=(c_ext[i][j]*c_ext[i][j]+c_ext[i-1][j]*c_ext[i-1][j])*.5f; //! avg of the two scattered field is required (E_eff) for the density update
            z2=(c_eyt[i][j]*c_eyt[i][j]+c_eyt[i][j-1]*c_eyt[i][j-1])*.5f; //! avg of the two scattered field is required (E_eff) for the density update
            c_erms2[i][j]=c_erms2[i][j]+(z1+z2)*c_inv_nperdt;  //! time updates and averages (child)
            
        }
    }

    if(k==2)
    {
        int j;
        
        for(i=0;i<child_elec->m;i++)
        {
            for(j=0;j<child_elec->n;j++)
            {
            	if (c_erms2[i][j]<0)
            	{
            		printf("Alert!!\n");
            	}
                child_elec->mesh[i][j] = sqrt(c_erms2[i][j]);  //! completes a period and then squre root the time avg data (child)
                c_erms2[i][j]=0.0f;
            } 
        }
    }
    gettimeofday(&end,NULL);    
    t_cal_child_rms += ((end.tv_sec - begin.tv_sec) + ((end.tv_usec - begin.tv_usec)/1000000.0));
}

void MR_MUR(int row)
{
    //! idea for the MUR boundary is to provide the E-field value at the four boundaries (top (ny)/bottom(0) and left (0)/right(nx))
    //! the values present at the +1 from left, -1 from right and -1from top and +1 from bottom from two time updates
   //! To avoid discontinuity due to truncation of the infinite boundary and absorb the scattered waves
    double csym;
    csym=1.0;
    if(IABSOR==2)
        csym=0.0;

    int column=0;
    if(row==1)
    {   
        for(column=0;column<ny;column++){
            eys[0][column]=eys1[1][column]+c1*(eys[1][column]-eys1[0][column]);
            eys1[0][column]=eys[0][column];  
            eys1[1][column]=eys[1][column];      
        } 
    }
    if(row==nx-1)
    {
        for(column=0;column<ny;column++)
        {
	        eys[nx][column]=eys1[nx-1][column]+c1*(eys[nx-1][column]-eys1[nx][column]);
	        eys1[nx][column]=eys[nx][column];
	        eys1[nx-1][column]=eys[nx-1][column];
        }   
    }
    exs[row][0]=exs1[row][1] +csym*c1*(exs[row][1]-exs1[row][0]);
    exs1[row][0]=exs[row][0];
    exs1[row][1]=exs[row][1];
    exs[row][ny]=exs1[row][ny-1] +csym*c1*(exs[row][ny-1]-exs1[row][ny]);
    exs1[row][ny]=exs[row][ny];
    exs1[row][ny-1]=exs[row][ny-1];
}

void child_MR_MUR_0(int row)
{
    double csym;
    csym=1.0;
    if(IABSOR==2)
        csym=0.0;

    int column=0;
    if(row==1)
    {   
        column=0;
        c_eys[0][column]=(c_eys1[1][column]+c_c1*(c_eys[1][column]-c_eys1[0][column]));
        c_eys1[0][column]=c_eys[0][column];  
        c_eys1[1][column]=c_eys[1][column];      
    }
    if(row==child_elec->m-1)
    {
        column=0;
        c_eys[child_elec->m][column]=(c_eys1[child_elec->m-1][column]+c_c1*(c_eys[child_elec->m-1][column]-c_eys1[child_elec->m][column]));
        c_eys1[child_elec->m][column]=c_eys[child_elec->m][column];
        c_eys1[child_elec->m-1][column]=c_eys[child_elec->m-1][column];
    }
    c_exs[row][0]=(c_exs1[row][1] +csym*c_c1*(c_exs[row][1]-c_exs1[row][0]));
    c_exs1[row][0]=c_exs[row][0];
    c_exs1[row][1]=c_exs[row][1];
}

void child_MR_MUR_1(int row)
{
    double csym;
    csym=1.0;
    if(IABSOR==2)
        csym=0.0;

    int column=0;
    if(row==child_elec->m-1)
    {
        for(column=0;column<child_elec->n;column++)
        {
	        c_eys[child_elec->m][column]=(c_eys1[child_elec->m-1][column]+c_c1*(c_eys[child_elec->m-1][column]-c_eys1[child_elec->m][column]));
	        c_eys1[child_elec->m][column]=c_eys[child_elec->m][column];
	        c_eys1[child_elec->m-1][column]=c_eys[child_elec->m-1][column];
        }   
	    c_exs[row][0]=(c_exs1[row][1] +csym*c_c1*(c_exs[row][1]-c_exs1[row][0]));
	    c_exs1[row][0]=c_exs[row][0];
	    c_exs1[row][1]=c_exs[row][1];
	    c_exs[row][child_elec->n]=(c_exs1[row][child_elec->n-1] +csym*c_c1*(c_exs[row][child_elec->n-1]-c_exs1[row][child_elec->n]));
	    c_exs1[row][child_elec->n]=c_exs[row][child_elec->n];
	    c_exs1[row][child_elec->n-1]=c_exs[row][child_elec->n-1];
	}
}

void child_MR_MUR_2(int row)
{
    double csym;
    csym=1.0;
    if(IABSOR==2)
        csym=0.0;

    int column=0;
    if(row==1)
    {   
        column=child_elec->n-1;
        c_eys[0][column]=(c_eys1[1][column]+c_c1*(c_eys[1][column]-c_eys1[0][column]));
        c_eys1[0][column]=c_eys[0][column];  
        c_eys1[1][column]=c_eys[1][column];      
    }
    if(row==child_elec->m-1)
    {
        column=child_elec->n-1;
        c_eys[child_elec->m][column]=(c_eys1[child_elec->m-1][column]+c_c1*(c_eys[child_elec->m-1][column]-c_eys1[child_elec->m][column]));
        c_eys1[child_elec->m][column]=c_eys[child_elec->m][column];
        c_eys1[child_elec->m-1][column]=c_eys[child_elec->m-1][column];
    }
    c_exs[row][child_elec->n]=(c_exs1[row][child_elec->n-1] +csym*c_c1*(c_exs[row][child_elec->n-1]-c_exs1[row][child_elec->n]));
    c_exs1[row][child_elec->n]=c_exs[row][child_elec->n];
    c_exs1[row][child_elec->n-1]=c_exs[row][child_elec->n-1];
}

void child_MR_MUR_3(int row)
{
    double csym;
    csym=1.0;
    if(IABSOR==2)
        csym=0.0;

    int column=0;
    if(row==1)
    {   
        for(column=0;column<child_elec->n;column++){
            c_eys[0][column]=(c_eys1[1][column]+c_c1*(c_eys[1][column]-c_eys1[0][column]));
            c_eys1[0][column]=c_eys[0][column];  
            c_eys1[1][column]=c_eys[1][column];      
        } 
    }
    
    if(row ==1)
    {
	    c_exs[row][0]=(c_exs1[row][1] +csym*c_c1*(c_exs[row][1]-c_exs1[row][0]));
	    c_exs1[row][0]=c_exs[row][0];
	    c_exs1[row][1]=c_exs[row][1];
	    c_exs[row][child_elec->n]=(c_exs1[row][child_elec->n-1] +csym*c_c1*(c_exs[row][child_elec->n-1]-c_exs1[row][child_elec->n]));
	    c_exs1[row][child_elec->n]=c_exs[row][child_elec->n];
	    c_exs1[row][child_elec->n-1]=c_exs[row][child_elec->n-1];
    }
}

double FIONIZ(double EE,double PR)
{
    //printf("In FIONIZ\n");
    double fioniz,ARG,VD;
    amu=QSM/FNUM*PRESSURE/PR;
    VD=amu*EE*PR;
    
    if(EE>2.0e4)
    {
	    ARG=BB1/EE;
	    fioniz=AA1*PR*exp(-ARG)*VD; 
    } 
    else if(EE<5.0e3)
    {
	    ARG=BB3*(1.0/EE-1.0/EE3);
	    fioniz=AA3*PR*(exp(-ARG)-1.0)*VD;
    } 
    else
    {
	    ARG=BB2/EE;
	    fioniz=AA2*PR*exp(-ARG)*VD;
    }
    return fioniz;
}

void ELEC_DENS()
{
	// printf("In parent elec dens %d\n",n);
    // printf("In ELEC_DENS\n");
    /*
    routine to calculate electron DENsity
    free diffusion + ambipolar dIFfusion + ioniZation
    */
    gettimeofday(&begin,NULL); 
    
    int ie,ied,je,jed,iii;
    double fioniz,aad,taumij;
    double coref,ee,cf,da,dac,fnui,fnua;
    double omgc2,rec1;
    double d0,dimax,ecm;
    double frqij,dtacmax,tcycle;

    //!=============================================================
    dnma=0.0;   //! max density at grids
    dimax=0.0;

    ACCEL=naccel; 
    if(ACCEL<=1) 
    ACCEL=1.0;

    coref=FNUM/sqrt(nu2omeg2);
    amu=EMOB;

    //c---------------------------------------------------------
    int i;
    
    
    for(i=1;i<nx;i++)
    {
        for(j=1;j<ny;j++)
        { 
            denp[i][j]=root_den->mesh[i][j];
            ee=root_elec->mesh[i][j];
            ee=ee/PRESSURE*coref;
            
            //c------------ calculation of ionization frequency---------
            frqio[i][j]=FIONIZ(ee,PRESSURE);    //ionization frequency
           
            //c---------------------------------------------------------
            
            cf=1.0;
            if(root_den->mesh[i][j]<=0.1) 
                cf=0.0;
            //          FNUIMAX=AMAX1(FNUIMAX,cf*frqio(I,J))
            frqij=frqio[i][j];

            //c----- calculation of diffusion coefficient--------------

            EDIF=dife;
            if(frqij<0)
                frqij=0.0;        

            if(cene==0.0){    //use formula for ETEM calc. Edif = diffusion of electron (free)
                ecm=ee*0.01;    
                ETEM=0.1+2.0/3.0*pow((0.0021*ecm*(91.0+ecm)),0.3333);
                EDIF=EMOB*ETEM;
            }

            if(cene<0.0)    //! to calculate the effective diffusion coefficient Edif, using a fixed ETEM=2 eV
            {
                taumij=eps0/(qe*(root_den->mesh[i][j]+1.0)*amu);
                aad=frqij*taumij;
                EDIF=(aad*dife+difa)/(aad+1.0);
            }
        
            //c---------------------------------------------------------              
            DIFFUSION[i][j]=EDIF;
         
            dimax=Max(dimax,DIFFUSION[i][j]);
            
        }
   	}
    //c---------------------------------------------------------              

    //c------- time step for fluid equation------------
    dtacmax=0.20*ds*ds/dimax;
    tcycle=1.0/FREQ;
    DTAC=Min(dtacmax, ndifmax*tcycle*ACCEL);

    for(i=1;i<root_den->m;i++)
    {
        for(j=1;j<root_den->n;j++)
        {
            da=DIFFUSION[i][j];
            dac=da*DTAC/(ds*ds);
            
            rec1=RECOMB*root_den->mesh[i][j]*DTAC;
            
            //c------------ ionization and attachment frequency-------------
            fnua=0.0;
            fnui=frqio[i][j];
            
            if(fnui<0.0)
            {
                fnua=-fnui;
                fnui=0.0;
            }
           
                if(isnan(denp[i][j]*exp(fnui*DTAC)))
            	{
            		printf("entered here%d\t%d\t%e\t%d\n",i,j,root_elec->mesh[i][j],n);
                	exit(0);
                }

            //c---------- Density equation updates------------------------

            //! parent mesh electron density update requires the previous densities (left, right, top and bottom and the center one)      
            root_den->mesh[i][j]=denp[i][j]*exp(fnui*DTAC)+dac*(denp[i+1][j]+denp[i-1][j]+denp[i][j+1]+denp[i][j-1]-4.0*denp[i][j]);
            root_den->mesh[i][j]=root_den->mesh[i][j]/(1.0+rec1+fnua*DTAC);

            if(root_den->mesh[i][j]<=1e-15)
                root_den->mesh[i][j]=root_den->mesh[i][j]*0.0;

            dnma=Max(dnma,root_den->mesh[i][j]);
        }
    
	}
	
    //c---------------Actual time calculation -------------------------
    TIMD=TIMD+DTAC;
    //c------------------------------------------------------------
    
    omgc2=qe*qe*inv_cmasse/eps0*dnma;
    PARC=omgc2/(pow(OMEG,2)+pow(FNUM,2));
  
    gettimeofday(&end,NULL);    
    t_cal_elec_dens += ((end.tv_sec - begin.tv_sec) + ((end.tv_usec - begin.tv_usec)/1000000.0));
}

void child_ELEC_DENS()
{
	// printf("In child elec dens %d\n",n);
  gettimeofday(&begin,NULL);
	int ie,ied,je,jed,iii;
    double fioniz,aad,taumij;
    double coref,ee,cf,da,dac,fnui,fnua;
    double omgc2,rec1;
    double d0,dimax,ecm;
    double frqij,dtacmax,tcycle;

    //!=============================================================
    dnma=0.0;   //! max density at grids
    dimax=0.0;

    ACCEL=naccel; 
    if(ACCEL<=1) 
    ACCEL=1.0;

    coref=FNUM/sqrt(nu2omeg2);
    amu=EMOB;

    //c---------------------------------------------------------
    int i;
    
    
    for(i=1;i<child_den->m;i++)
    {
        for(j=1;j<child_den->n;j++)
        { 
        	c_denp[i][j]=child_den->mesh[i][j];
            ee=child_elec->mesh[i][j];
            ee=ee/PRESSURE*coref;
            
            //c------------ calculation of ionization frequency---------
            c_frqio[i][j]=FIONIZ(ee,PRESSURE);   //! the ionization frequency
            //c---------------------------------------------------------
            
            cf=1.0;
            if(child_den->mesh[i][j]<=0.1) 
                cf=0.0;
            //          FNUIMAX=AMAX1(FNUIMAX,cf*frqio(I,J))
            frqij=c_frqio[i][j];

            //c----- calculation of diffusion coefficient--------------

            EDIF=dife;
            if(frqij<0)
                frqij=0.0;        

            if(cene==0.0){
                ecm=ee*0.01;    
                ETEM=0.1+2.0/3.0*pow((0.0021*ecm*(91.0+ecm)),1/3);
                EDIF=EMOB*ETEM;
            }

            if(cene<0.0)
            {
                taumij=eps0/(qe*(child_den->mesh[i][j]+1.0)*amu);
                aad=frqij*taumij;
                EDIF=(aad*dife+difa)/(aad+1.0);
            }
        
            //c---------------------------------------------------------              
            c_DIFFUSION[i][j]=EDIF;
         
            dimax=Max(dimax,c_DIFFUSION[i][j]);
            
        }
   	}
    //c---------------------------------------------------------              

    //c------- time step for fluid equation------------
    dtacmax=0.20*c_ds*c_ds/dimax;
    tcycle=1.0/FREQ;
    DTAC=Min(dtacmax, ndifmax*tcycle*ACCEL);

    for(i=1;i<child_den->m;i++)
    {
        for(j=1;j< child_den->n;j++)
        {
            da=c_DIFFUSION[i][j];
            dac=da*DTAC/(c_ds*c_ds);
            
            rec1=RECOMB*child_den->mesh[i][j]*DTAC;
            
            //c------------ ionization and attachment frequency-------------
            fnua=0.0;
            fnui=c_frqio[i][j];
            if(fnui<0.0)
            {
                fnua=-fnui;
                fnui=0.0;
            }

            if(isnan(c_denp[i][j]*exp(fnui*DTAC)))
            {
            	printf("Child, n = %d\n", n);
            	exit(0);
            }

            //c---------- Density equation updates------------------------
            //! child mesh electron density update requires the previous densities (left, right, top and bottom and the center one)      
            child_den->mesh[i][j]=c_denp[i][j]*exp(fnui*DTAC)+dac*(c_denp[i+1][j]+c_denp[i-1][j]+c_denp[i][j+1]+c_denp[i][j-1]-4.0*c_denp[i][j]);
            


            child_den->mesh[i][j]=child_den->mesh[i][j]/(1.0+rec1+fnua*DTAC);
            
            c_denpold[i][j]=child_den->mesh[i][j];      //! reqd when expanding the child mesh to copy the previous mesh data
            
            
            if(child_den->mesh[i][j]<=1e-15)
                child_den->mesh[i][j]=child_den->mesh[i][j]*0.0;

            dnma=Max(dnma,child_den->mesh[i][j]);
        }
    
	}

    //c---------------Actual time calculation -------------------------
    TIMD=TIMD+DTAC;
    //c------------------------------------------------------------
    
    omgc2=qe*qe*inv_cmasse/eps0*dnma;
    PARC=omgc2/(pow(OMEG,2)+pow(FNUM,2));
    
    gettimeofday(&end,NULL);    
    t_cal_child_elec_dens += ((end.tv_sec - begin.tv_sec) + ((end.tv_usec - begin.tv_usec)/1000000.0));
}

void SETUP()
{
    // printf("In SETUP\n");
    static int k;
    static double ardiy,ardix,dinig;
    static double xd0,yd0,xxi,yyj;
    static double c_ardiy,c_ardix,c_dinig;
    static double c_xd0,c_yd0,c_xxi,c_yyj;
    static int xtotal,ytotal;
    //static int xend,yend,xstar,ystar;
//!---------------------------divisible by 4(child data)-----    
    //xend=(9*nx/10);yend=(9*ny/10);xstar=(1*nx/2);ystar=(1*ny/10);
    xend=(xe*nx);yend=(ye*ny);xstar=(xs*nx);ystar=(ys*ny);
    /*  
     xend = 0.65*nx;
    xstar = 0.35*nx;
    yend = 0.9*ny;
    ystar = 0.1*ny;
  */
  /*
     xend = 0.8*nx;
    xstar = 0.2*nx;
    yend = 0.9*ny;
    ystar = 0.1*ny;
  */
   /* 
    xend = 0.95*nx;
    xstar = 0.05*nx;
    yend = 0.9*ny;
    ystar = 0.1*ny;
    */
    
   printf("Before,xend=%d\tyend=%d\txstar=%d\tystar=%d\t\n",xend,yend,xstar,ystar);
	
    xtotal=(xend-xstar)*factor;
    ytotal=(yend-ystar)*factor;
    printf("after,xend=%d\tyend=%d\txstar=%d\tystar=%d\t\n",xend,yend,xstar,ystar);
//!--------------------------------------------------------------
	
		
  root_elec = newnode(nx, ny, 0, 0, 0);
	child_elec = newnode(xtotal, ytotal, xstar, ystar, 1);
	
  root_elec->mesh = root_mesh_elec;
	child_elec->mesh = child_mesh_elec;

	root_elec->children[0] = child_elec;
	root_elec->parent = NULL;

	child_elec->parent = root_elec;


	root_mag = newnode(nx, ny, 0, 0, 0);
	child_mag = newnode(xtotal, ytotal, xstar, ystar, 1);

	root_mag->mesh = root_mesh_mag;
	child_mag->mesh = child_mesh_mag;

	root_mag->children[0] = child_mag;
	root_mag->parent = NULL;

	child_mag->parent = root_mag;


	root_den = newnode(nx, ny, 0, 0, 0);
	child_den = newnode(xtotal, ytotal, xstar, ystar, 1);

	root_den->mesh = root_mesh_den;
	child_den->mesh = child_mesh_den;

	root_den->children[0] = child_den;
	root_den->parent = NULL;

	child_den->parent = root_den;


    dt=1.0/(double)(nperdt)/FREQ;
    ds=c/(double)(nlamb)/FREQ;

    c_dt=dt/factor;
    c_ds=ds/factor;

    OMEG=2.0*pi*FREQ;
    QSM=qe/cmasse;
    e2_epsm=qe*qe/eps0/cmasse;
    
    /*
    !  AIR DAT
    !  FNUM=electron-neutral coll frequency
    !  RECOMB=electron-ion recombination coefficient
    !  EMOB=electron mobility
    !  ETEM=electron temperature
    !  DIFE=electron diffusion coeff
    */

    FNUM=5.3*pow(10,9)*PRESSURE;
    RECOMB=crec*1.0*pow(10,-13);
    EMOB=QSM/FNUM;
    ETEM=2.0*abs(cene);
    dife=EMOB*ETEM;
    difa=dife/100.0;
    nu2omeg2=pow(OMEG,2)+pow(FNUM,2);
    nu_omeg=FNUM/OMEG;

    //!================Initial density location =======================

    for(K=1;K<=nini;K++)
    {
	      imid[K]=0;
	      jmid[K]=0;
        c_imid[K]=0;
	      c_jmid[K]=0;
	   
	      if(xmid[K]>0.0) 
        {
	        
        imid[K] = xpos*nx;
        
        c_imid[K] = xpos*(factor*(nx));
        
                         printf("xmid[%d]=%f\n",K,xmid[K] ); 
                         printf("imid[%d]=%d\n",K,imid[K] );
                         printf("c_imid[%d]=%d\n",K,c_imid[K] );
	      } 
	      if(ymid[K]>0.0) 
        {  
	        
          jmid[K] = ypos*ny;
	         
          c_jmid[K] = ypos*(factor*ny);
                  
                 printf("ymid[%d]=%f\n",K,ymid[K] );
                 printf("jmid[%d]=%d\n",K,jmid[K] );
                 printf("c_jmid[%d]=%d\n",K,c_jmid[K] );
	      } 
	   if(xmid[K]==0.0) 
	   {     imid[K] = nx/2.0;
	        c_imid[K] = factor/2.0*nx;
	    }
      if(ymid[K]==0.0) 
	     {
          jmid[K] =ny/2.0;
	        c_jmid[K] = factor/2.0*ny;
      }
    }

    //!=============================================================
    TEMP0=300.0;
    DENG0=PRESSURE/760.0*101300.0/akb/TEMP0;

    radius = nx/5*ds;
    printf("radius: %f\n", radius);
    printf("ds: %f\t c_ds: %f\n", ds,c_ds);
    int i;
    
   

	for (i = 0; i <= root_elec->m; i+=1)
	{
		for (j = 0; j <= root_elec->n; j+=1)
		{
			root_den->mesh[i][j] = 0.0;
			root_elec->mesh[i][j] = E0/sqrt(2.0);
		}
	}

	for (i = 0; i <= child_elec->m; i+=1)
	{
		for (j = 0; j <=child_elec->n; j+=1)
		{
			child_den->mesh[i][j] = 0;
			child_elec->mesh[i][j] = E0/sqrt(2.0);
		}
	}


    for(K=1;K<=nini;K++)
    {
	    xd0=ds*imid[K];
        // printf("%.20f\t%.20f\t%d\n",xd0,ds,imid[K]);

	    yd0=ds*jmid[K];
	    c_xd0=c_ds*c_imid[K];
	    c_yd0=c_ds*c_jmid[K];
	    //printf("%f %f\n",xd0,yd0 );
	    if(xmid[K]<0) 
	        xd0=-xmid[K];
	    if(ymid[K]<0) 
	        yd0=-ymid[K];

	    if(xmid[K]<0) 
	        c_xd0=-xmid[K];
	    if(ymid[K]<0) 
	        c_yd0=-ymid[K];

	    //!================Initial density, Gaussian, defined =======================
	    //!make DEN and DENP =0
	    int j;
        
	    for(i=0;i<=nx;i++)
	    {
	        xxi=ds*i;
	        ardix=0.0;
	        if(sgdx0[K]>0)
	            ardix=(-pow((xxi-xd0),2))/2.0/sgdx0[K]/sgdx0[K];

            for(j=0;j<=ny;j++){    
                yyj=ds*j;
                ardiy=0.0;
                if(sgdy0[K]>0) 
                   ardiy=-pow((yyj-yd0),2)/2.0/sgdy0[K]/sgdy0[K];
                dinig=DINI[K]*exp(ardix+ardiy);
                //if(dinig<=DINI[K]*1.0*exp(-2)) 
                  if(dinig<=1.0e13)
                    dinig=0;
                // den[i][j]=den[i][j]+dinig;
                root_den->mesh[i][j] = root_den->mesh[i][j] + dinig;

            }
	    }
		int ic=xpos*nx,jc=ypos*ny;
		printf("root_den->mesh[%d][%d]=%e \n",ic,jc, root_den->mesh[ic][jc]);
			
	    // exit(0);

	    for(i=0;i<=child_den->m;i+=1)
	    {
	    	c_xxi=c_ds*(i+factor*child_den->locx);
	    	c_ardix=0.0;
	    	if(sgdx0[K]>0)
	            c_ardix=(-pow((c_xxi-c_xd0),2))/2.0/sgdx0[K]/sgdx0[K];
	    	for(j=0;j<=child_den->n;j+=1)
	    	{
	    		// printf("%d\n",child_den->n );
	    		c_yyj=c_ds*(j+factor*child_den->locy);
                c_ardiy=0.0;
                if(sgdy0[K]>0){
                	// printf("Happens!!\n");
                   c_ardiy=-pow((c_yyj-c_yd0),2)/2.0/sgdy0[K]/sgdy0[K];
                }
                dinig=DINI[K]*exp(c_ardix+c_ardiy);
                //if(dinig<=DINI[K]*1.0*exp(-2)) 
                  if(dinig<=1.0e13)
                    dinig=0;
                child_den->mesh[i][j] = child_den->mesh[i][j] + dinig;
	    	}
	    }
    }
    dte=dt/eps0;
    dtm=dt/xmu0;
    dteds=dte/ds;
    dtmds=dtm/ds;

    c_dte=c_dt/(eps0);
    c_dtm=c_dt/(xmu0);
    c_dteds=c_dte/c_ds;
    c_dtmds=c_dtm/c_ds;
  
    //mur constants
    c1=(c*dt-ds)/(c*dt+ds);
    c2=2.0*ds/(c*dt+ds);
    c3=(c*dt*c*dt)/(2.0*ds*(c*dt+ds));

    c_c1 = (c*c_dt-c_ds)/(c*c_dt+c_ds);
    
           panim();
           panimE();
           canim();
		       canimE();
    
    
}

void SETUP2()
{

    static int k,contr=-1,contr2=-1;
    static double ardiy,ardix,dinig;
    static double xd0,yd0,xxi,yyj;
    static double c_ardiy,c_ardix,c_dinig;
    static double c_xd0,c_yd0,c_xxi,c_yyj;
    
    static int xtotal2,ytotal2;
    

  if(n==0){
     
   //!   The code should be run for factor =2 /4 not more refinement for coarse nlamb= 256/128 minimum for accuracy,
  //! for  faster experiments to check whether parallelization working , nlamb=64 can be used not below, use E0=5e6 instead of 5.5e6
  //!  to properly resolve the field and density
  //! Code only works for power of 2 (refinement factor only 2,4) and nlamb = 64, 128, 256 ... (parent ) child (512 min for accuracy)
  if(nlamb<64){
  printf("minimum cells required for accuracy is nlamb=128 factor=4 use nlamb=64 for small experiments\n");
  exit(0);
  }
   
     	
          //!------------Activate for MR ---------------- //! for static mesh refinement (MR)
   if(option==0)
   { xend=(xe*nx);yend=(ye*ny);xstar=(xs*nx);ystar=(ys*ny);}
  
  
   //!--------------------Activate for dynamic MR --------------- 
   if(option==1)
   {
      //! The expanding box grows towards left in x, top and bottom along y, 
      //!the right side expansion along x is restricted as no filament propagation

       inibxs= (int)(bxsize*nlamb);  //!inibxs: initial box size x
       inibys= (int)(bysize*nlamb);  //!inibys: initial box size y
       
       xstar=(int)(xpos*nx)-(int)(inibxs*0.5); //!inibxs: initial box size x, xstar: start cell of the child mesh location on parent 
       
       xend= (int)(xpos*nx)+(int)(inibxs*0.5); //!inibxs: initial box size x  xend: end cell of the child mesh location on parent 
       
       ystar=(int)(ypos*ny)-(int)(inibys*0.5); //!inibys: initial box size y  ystar: start cell of the child mesh location on parent 
       
       yend= (int)(ypos*ny)+(int)(inibys*0.5); //!inibys: initial box size y  yend: end cell of the child mesh location on parent
       
       xstarfix=xstar;    //! use for calculating fixed cells to increament along  x (towards left : source side)
       ystarfix=ystar;    //! use for calculating fixed cells to increament along y (towards bottom)
       yendfix=yend;      //! use for calculating fixed cells to increament along y (towards top)
       xstarol=xstar;   //! previous start cell number of initial box along x (left)
       ystarol=ystar;   //! previous start cell number of initial box along y (bottom)
       yendol=yend;     //! previous start cell number of initial box along y  (top)
   
       xinc=(int)(incrx*xstarfix);
       yincs=(int)(incry*ystarfix);
       //yincn=(int)(incry*yendfix);
       
       printf("initial size of box x =%d\n",inibxs);
       printf("initial size of box y =%d\n",inibys);
       printf("Before,xend=%d\tyend=%d\txstar=%d\tystar=%d\t\n",xend,yend,xstar,ystar);
       printf("xstarol=%d ystarol=%d yendol=%d\n",xstarol,ystarol,yendol); 
       printf("xstarfix=%d ystarfix=%d yendfix=%d\n",xstarfix,ystarfix,yendfix); 
       printf("No.expboxcells xs = %d ys = %d yen = %d\n",xinc,yincs,yincs);
       //!----------------------------------------------------------------------- 
    
    //! The cells calculation along x and y for two conditions (a: minimum 25 cells required x and y (each side) expansion)
    //!                       (b: minimum nx>=128 (coarse and fine (factor =4) :512 min. cells per lambda accuracy))
    
    if(xinc>=30){
     xinc=(int)(incrx*xstarfix); 
     //xinc=40;  //activate if bigger domain required
     }
     else{
    
        if(nx>=128)
         {xinc=30;}
         else
         {xinc=20;}
    }
     
     if(yincs>=30){
     yincs=(int)(incry*ystarfix);
     //yincs=40;    //!activate if bigger domain required
     }
     else{
     
     if(ny>=128)
       {yincs=30;}
     else
       {yincs=20;}
     }
    
    }
    xtotal2=(xend-xstar)*factor;      //! Total cells in the refine mesh x-direction
    ytotal2=(yend-ystar)*factor;      //! Total cells in the refine mesh y-direction
    
    printf("after,xend=%d\tyend=%d\txstar=%d\tystar=%d\t\n",xend,yend,xstar,ystar);
//!--------------------------------------------------------------
	
		//!========== Mesh generation (initial) =====================

    //!---------E-field----------
    root_elec = newnode(nx, ny, 0, 0, 0);
	
    child_elec = newnode(xtotal2, ytotal2, xstar, ystar, 1);
  	
  	root_elec->mesh = root_mesh_elec;
  	child_elec->mesh = child_mesh_elec;
  
  	root_elec->children[0] = child_elec;
  	root_elec->parent = NULL;
  
  	child_elec->parent = root_elec;

   //!---------H-field----------
  	root_mag = newnode(nx, ny, 0, 0, 0);
  	child_mag = newnode(xtotal2, ytotal2, xstar, ystar, 1);
  
  	root_mag->mesh = root_mesh_mag;
  	child_mag->mesh = child_mesh_mag;
  
  	root_mag->children[0] = child_mag;
  	root_mag->parent = NULL;
  
  	child_mag->parent = root_mag;

     //!---------electron density----------
  	root_den = newnode(nx, ny, 0, 0, 0);
  	child_den = newnode(xtotal2, ytotal2, xstar, ystar, 1);
  
  	root_den->mesh = root_mesh_den;
  	child_den->mesh = child_mesh_den;
  
  	root_den->children[0] = child_den;
  	root_den->parent = NULL;
  
  	child_den->parent = root_den;
   
  //! ============== parent and child time and grid step (space), constants=======================================
    dt=1.0/(double)(nperdt)/FREQ;
    ds=c/(double)(nlamb)/FREQ;

    c_dt=dt/factor;
    c_ds=ds/factor;

    OMEG=2.0*pi*FREQ;
    QSM=qe/cmasse;
    e2_epsm=qe*qe/eps0/cmasse;
   
   
    FNUM=5.3*pow(10,9)*PRESSURE;
    RECOMB=crec*1.0*pow(10,-13);
    EMOB=QSM/FNUM;
    ETEM=2.0*abs(cene);
    dife=EMOB*ETEM;
    difa=dife/100.0;
    nu2omeg2=pow(OMEG,2)+pow(FNUM,2);
    nu_omeg=FNUM/OMEG;
    
    TEMP0=300.0;
    DENG0=PRESSURE/760.0*101300.0/akb/TEMP0;

    // int i;
     
    //  for (i = 0; i <= root_elec->m; i+=1)
  	// {
  	// 	for (j = 0; j <= root_elec->n; j+=1)
  	// 	{
  	// 		root_den->mesh[i][j] = 0.0;
  	// 		root_elec->mesh[i][j] = E0/sqrt(2.0);
  	// 	}
  	// }
  
  	// for (i = 0; i <= child_elec->m; i+=1)
  	// {
  	// 	for (j = 0; j <=child_elec->n; j+=1)
  	// 	{
  	// 		child_den->mesh[i][j] = 0;
  	// 		child_elec->mesh[i][j] = E0/sqrt(2.0);
  	// 	}
  	// }
      cudaMalloc((void**)&dev_root_elec, sizeof(root_elec));
      cudaMalloc((void**)&dev_den, sizeof(root_den));
      cudaMalloc((void**)&dev_E0, sizeof(double));
      cudaMemcpy(dev_root_elec, root_elec, sizeof(root_elec), cudaMemcpyHostToDevice);
      cudaMemcpy(dev_E0, &E0, sizeof(double), cudaMemcpyHostToDevice);
      setup_init<<<(ceil(root_elec->m/32),ceil(root_elec->n/32)),(32,32)>>>(dev_root_elec,dev_den,dev_E0);    
    cudaMemcpy(root_den, dev_den, sizeof(dev_den), cudaMemcpyDeviceToHost);
    cudaMemcpy(root_elec, dev_root_elec, sizeof(dev_den), cudaMemcpyDeviceToHost);
    cudaFree(dev_root_elec);
    cudaFree(dev_den);
    cudaFree(dev_E0);

    cudaMalloc((void**)&dev_root_elec, sizeof(child_elec));
      cudaMalloc((void**)&dev_den, sizeof(child_den));
      cudaMalloc((void**)&dev_E0, sizeof(double));
      cudaMemcpy(dev_root_elec, child_elec, sizeof(child_elec), cudaMemcpyHostToDevice);
      cudaMemcpy(dev_E0, &E0, sizeof(double), cudaMemcpyHostToDevice);
      setup_init<<<(ceil(child_elec->m/32),ceil(child_elec->n/32)),(32,32)>>>(dev_root_elec,dev_den,dev_E0);    
    cudaMemcpy(child_den, dev_den, sizeof(dev_den), cudaMemcpyDeviceToHost);
    cudaMemcpy(child_elec, dev_root_elec, sizeof(dev_den), cudaMemcpyDeviceToHost);
    cudaFree(dev_root_elec);
    cudaFree(dev_den);
    cudaFree(dev_E0);

            
      //!========================= centering the quantities in the new mesh refinement region =============
      
      //!================Initial density location =======================

    for(K=1;K<=nini;K++)
    {
	      imid[K]=0;
	      jmid[K]=0;
        c_imid[K]=0;
	      c_jmid[K]=0;
	   
	      if(xmid[K]>0.0) 
        {
	        
        imid[K] = xpos*(nx);
        
        c_imid[K] = xpos*(factor*nx);         
                        
	      } 
	      if(ymid[K]>0.0) 
        {  
	        
          jmid[K] = ypos*ny;
	         
          c_jmid[K] = ypos*(factor*ny);          
                 
	      } 
    
    }

    //!=============================================================
          
       for(K=1;K<=nini;K++)
    {
	    xd0=ds*imid[K];
      yd0=ds*jmid[K];
	    c_xd0=c_ds*c_imid[K];
	    c_yd0=c_ds*c_jmid[K];
	    
      
	    //!================Initial density, Gaussian, defined =======================
	    
	   int j;

	    for(i=0;i<=nx;i++)
	    {
	        xxi=ds*i;
	        ardix=0.0;
	        if(sgdx0[1]>0)
	           ardix=(-pow((xxi-xd0),2))/2.0/sgdx0[1]/sgdx0[1];

            for(j=0;j<=ny;j++){    
                yyj=ds*j;
                ardiy=0.0;
                if(sgdy0[K]>0) 
                  ardiy=-pow((yyj-yd0),2)/2.0/sgdy0[K]/sgdy0[K];
                 dinig=DINI[K]*exp(ardix+ardiy);
                 if(dinig<=1.0e13)
                   dinig=0;
                 
                root_den->mesh[i][j] = root_den->mesh[i][j]+ dinig;
                denp[i][j]=root_den->mesh[i][j];
            }
	    }

        

        // cudaMalloc((void**)&dev_den, sizeof(root_den));
        // cudaMalloc((void**)&dev_ny, sizeof(int));
        // cudaMalloc((void**)&dev_nx, sizeof(int));
        // cudaMalloc((void**)&dev_xxi, sizeof(double));
        // cudaMalloc((void**)&dev_ds, sizeof(double));
        // cudaMalloc((void**)&dev_ardix, sizeof(double));
        // cudaMalloc((void**)&dev_yyj, sizeof(double));
        // cudaMalloc((void**)&dev_ardiy, sizeof(double));
        // cudaMalloc((void**)&dev_yd0, sizeof(double));
        // cudaMalloc((void**)&dev_dinig, sizeof(double));
        // cudaMalloc((void**)&dev_sgdx0, sizeof(sgdx0));
        // cudaMalloc((void**)&dev_sgdy0, sizeof(dev_sgdy0));
        // cudaMalloc((void**)&dev_DINI, sizeof(dev_DINI));
        // cudaMalloc((void**)&dev_K, sizeof(K));

        // cudaMemcpy(root_den, dev_den, sizeof(dev_den), cudaMemcpyDeviceToHost);
        // cudaMemcpy(dev_ny, &ny, sizeof(ny), cudaMemcpyDeviceToHost);
        // cudaMemcpy(dev_nx, &nx, sizeof(nx), cudaMemcpyDeviceToHost);
        // cudaMemcpy(dev_xxi, &xxi, sizeof(xxi), cudaMemcpyDeviceToHost);
        // cudaMemcpy(dev_ds, &ds, sizeof(ds), cudaMemcpyDeviceToHost);
        // cudaMemcpy(dev_ardix, &ardix, sizeof(ardix), cudaMemcpyDeviceToHost);
        // cudaMemcpy(dev_yyj, &yyj, sizeof(yyj), cudaMemcpyDeviceToHost);
        // cudaMemcpy(dev_ardiy, &ardiy, sizeof(ardiy), cudaMemcpyDeviceToHost);
        // cudaMemcpy(dev_xd0, &xd0, sizeof(xd0), cudaMemcpyDeviceToHost);
        // cudaMemcpy(dev_yd0, &yd0, sizeof(yd0), cudaMemcpyDeviceToHost);
        // cudaMemcpy(dev_dinig, &dinig, sizeof(dinig), cudaMemcpyDeviceToHost);
        // cudaMemcpy(dev_sgdx0, sgdx0, sizeof(sgdx0), cudaMemcpyDeviceToHost);
        // cudaMemcpy(dev_sgdy0, sgdy0, sizeof(sgdy0), cudaMemcpyDeviceToHost);
        // cudaMemcpy(dev_DINI, DINI, sizeof(DINI), cudaMemcpyDeviceToHost);
        // cudaMemcpy(dev_K, &K, sizeof(K), cudaMemcpyDeviceToHost);

        
        // setup_init1<<<(ceil(nx/32),ceil(ny/32)),(32,32)>>>(dev_den, dev_ny, dev_nx, dev_xxi, dev_ds, dev_ardix,dev_yyj,dev_ardiy,dev_xd0,dev_yd0,dev_dinig,dev_sgdx0,dev_sgdy0,dev_DINI,dev_K);

        // cudaMemcpy(root_den, dev_den, sizeof(dev_den), cudaMemcpyDeviceToHost);
        // cudaMemcpy(denp, dev_den->mesh, sizeof(dev_den->mesh), cudaMemcpyDeviceToHost);

        // cudaFree(dev_den);
        // cudaFree(dev_ny);
        // cudaFree(dev_nx);
        // cudaFree(dev_xxi);
        // cudaFree(dev_ds);
        // cudaFree(dev_ardix);
        // cudaFree(dev_yyj);
        // cudaFree(dev_ardiy);
        // cudaFree(dev_xd0);
        // cudaFree(dev_yd0);
        // cudaFree(dev_dinig);
        // cudaFree(dev_sgdx0);
        // cudaFree(dev_sgdy0);
        // cudaFree(dev_DINI);
  
          printf("Entered here \n");
          interpolatexinitial(1.0);    //! to interpolate the initial gaussian density parent to child location
          printf("Entered here \n");
          
          printf("left here \n");
    }
      
    
    dte=dt/eps0;
    dtm=dt/xmu0;
    dteds=dte/ds;
    dtmds=dtm/ds;

    c_dte=c_dt/(eps0);
    c_dtm=c_dt/(xmu0);
    c_dteds=c_dte/c_ds;
    c_dtmds=c_dtm/c_ds;
  
    //mur constants
    c1=(c*dt-ds)/(c*dt+ds);
    c2=2.0*ds/(c*dt+ds);
    c3=(c*dt*c*dt)/(2.0*ds*(c*dt+ds));

    c_c1 = (c*c_dt-c_ds)/(c*c_dt+c_ds);
    
    //!======== Initialize the child electric field , velcocity, density and magnetic field in the new refinement region ========

          
  
    //!=============================================================
    
      
           printf("Entered here \n");
           panim();
           panimE();
           canim();
		       canimE(); 
 
  }
 if(n>0){  
     
    if(option==1){
    if(chc==1)
   {
      //! cells expansion in x , must not go below the iend or xs*nx
     if(xstarol>xinc)
     {
     xstarnew=xstarol-xinc; 
    
     if(xstarnew>=(int)(xs*nx))
     { xstar=xstarnew;}
     else
      {xstar=(int)(xs*nx);}
     
     }  
     else
     {
       xstar=(int)(xs*nx);
       xstarnew=xstar;
     }
   
   }
   
   if(chc==2)
   {
      //! cells expansion in y , must not go below the ystar (final) or ys*ny
     if(ystarol>yincs)
     {
       ystarnew=ystarol-yincs; 
       
       if(ystarnew>=(int)(ys*ny))
       { ystar=ystarnew;}
       else
       { ystar=(int)(ys*ny);} 
     }
     else
     {
        ystar=(int)(ys*ny);
        ystarnew=ystar;
     }
   }
  
   if(chc==2)
   {
      //! cells expansion in y , must not go above the ystar (final) or ye*ny
     if(yendol<=(int)(ye*ny))
     {
       yendnew=yendol+yincs; 
       
       if(yendnew<=(int)(ye*ny))
        {yend=yendnew;}
       else
        {yend=(int)(ye*ny);}
     }
     else
     {
        yend=(int)(ye*ny);
        yendnew=yend;
     }
   
   }
   
    xtotal2=(xend-xstar)*factor;
    ytotal2=(yend-ystar)*factor;
    printf("after,xend=%d\tyend=%d\txstar=%d\tystar=%d\t\n",xend,yend,xstar,ystar);
//!--------------------------------------------------------------
	
	 //!========== Mesh generation (expanding) =====================
    //! every steps follows same technique: mesh generation, define the parent and child relation , here only one child 
    //! no further division exists
    //!---------E-field----------	
    root_elec = newnode(nx, ny, 0, 0, 0);
	
    child_elec = newnode(xtotal2, ytotal2, xstar, ystar, 1);
  	
  	root_elec->mesh = root_mesh_elec;
  	child_elec->mesh = child_mesh_elec;
  
  	root_elec->children[0] = child_elec;
  	root_elec->parent = NULL;
  
  	child_elec->parent = root_elec;

    //!---------H-field---------- 
  	root_mag = newnode(nx, ny, 0, 0, 0);
  	child_mag = newnode(xtotal2, ytotal2, xstar, ystar, 1);
  
  	root_mag->mesh = root_mesh_mag;
  	child_mag->mesh = child_mesh_mag;
  
  	root_mag->children[0] = child_mag;
  	root_mag->parent = NULL;
  
  	child_mag->parent = root_mag;

    //!---------electron density---------- 
  	root_den = newnode(nx, ny, 0, 0, 0);
  	child_den = newnode(xtotal2, ytotal2, xstar, ystar, 1);
  
  	root_den->mesh = root_mesh_den;
  	child_den->mesh = child_mesh_den;
  
  	root_den->children[0] = child_den;
  	root_den->parent = NULL;
  
  	child_den->parent = root_den;
  
   //!=======================================
    printf("Entered here \n");
    //! chc==1/2 indicate the interpolation correspond to x or y expansion
    if(chc==1)
    {
      interpolatexnew(1.0);   //! interpolates all previous mesh data along x
      printf("Entered here \n");
    }
    if(chc==2)
    {
      interpolateynew(1.0);   //! interpolates all previous mesh data along y
      printf("Entered here \n");
    }
    
    printf("left here \n");

    if(chc==1)
   {
       if(xstarnew>(int)(xs*nx))
       xstarol=xstarnew;
       
   }
   
   if(chc==2)
   {
        if(ystarnew>(int)(ys*ny))
        { ystarol=ystarnew;}
        
        if(yendnew<(int)(ye*ny))
        { yendol=yendnew;}
        
    } 
     
    
  }
 
 }                  

}


//! to perform interpolation of parent data f(x,y): f represents either E-,H-, electron density etc., on child at 2d positions 
//! bilinear interpolation
double interpolate2d(double a,double b,double c,double d, double posx, double posy)
{
	double val;
	val = (1-posx)*(1-posy)*a + (posx)*(1-posy)*b + (posx)*posy*c + (1-posx)*posy*d;
	return val;
}

//! to perform interpolation of parent data f(x) or f(y): f represents either E-,H-, electron density etc., on child at 1d positions 
//! linear interpolation
double interpolate1d(double a, double b, double posx)
{
	return ((1-posx)*a + b*posx);
}

//! perform interpolation at the corners (density), at center (H-field), edge (E-field) only on the cpb (child parent boundary)
void interpolatecorners(double t1)
{
	double tau = t1;
	int basex, basey;
	double distx, disty;
	//int factor=8;
	if(child_elec->locx>0)
	{
		basex = child_elec->locx-1;
		basey = child_elec->locy;
		distx =	0.5+ 0.5/factor;
    	c_exs[0][0] = (1-tau)*(interpolate1d(exs_old[basex][basey],
			exs_old[basex+1][basey],distx))+(tau)*(interpolate1d(
			exs[basex][basey],exs[basex+1][basey],distx));
    					
    	c_vx[0][0] = (1-tau)*(interpolate1d(vx_old[basex][basey],
		     vx_old[basex+1][basey],distx))+(tau)*(interpolate1d(
		     vx[basex][basey],vx[basex+1][basey],distx));
	}

	if(child_elec->locy>0)
	{
		basex = child_elec->locx;
		basey = child_elec->locy-1;
		disty = 0.5+0.5/factor;
		c_eys[0][0] = (1-tau)*interpolate1d(eys_old[basex][basey],
			       eys_old[basex][basey+1],disty) + 
			      (tau)*interpolate1d(eys[basex][basey],
			       eys[basex][basey+1],disty);  
		c_vy[0][0] = (1-tau)*interpolate1d(vy_old[basex][basey],
			     vy_old[basex][basey+1],disty) + (tau)*interpolate1d(			      vy[basex][basey],vy[basex][basey+1],disty);
	}

	if(child_elec->locx>0&&child_elec->locy>0)
	{
		basex = child_elec->locx-1;
		basey = child_elec->locy-1;
		distx = 0.5 + 0.5/factor;
		disty = 0.5 + 0.5/factor;
		child_mag->mesh[0][0] = (1-tau)*(interpolate2d(hzi[basex][basey],					  hzi[basex+1][basey],
					  hzi[basex+1][basey+1],
					 hzi[basex][basey+1],distx,disty))+
					(tau)*(interpolate2d(
					 root_mag->mesh[basex][basey],
					root_mag->mesh[basex+1][basey],
					root_mag->mesh[basex+1][basey+1],
					root_mag->mesh[basex][basey+1],
					distx,disty));
	}

	child_den->mesh[0][0] = (1-tau)*denp[child_den->locx][child_den->locy] + 			     tau*root_den->mesh[child_den->locx][child_den->locy];

	basex = child_elec->locx+(child_elec->m)/factor;
	basey = child_elec->locy-1;
	disty = 0.5 + 0.5/factor;
	c_eys[child_elec->m][0] = (1-tau)*interpolate1d(eys_old[basex][basey],
				  eys_old[basex][basey+1],disty) + 
				  tau*interpolate1d(eys[basex][basey],
				  eys[basex][basey+1],disty);

	basex = child_elec->locx-1;
	basey = child_elec->locy+(child_elec->n)/factor;
	distx = 0.5 + 0.5/factor;
	c_exs[0][child_elec->n] = (1-tau)*interpolate1d(exs_old[basex][basey],
				exs_old[basex+1][basey],distx) + 
				tau*interpolate1d(exs[basex][basey],
				exs[basex+1][basey],distx);

	c_eyt[0][0] = c_eys[0][0] + c_eyi1[0][0];
	c_ext[0][0] = c_exs[0][0] + c_exi1[0][0];

c_ext[0][child_elec->n]= c_exs[0][child_elec->n] + c_exi1[0][child_elec->n];
c_eyt[child_elec->m][0] = c_exs[child_elec->m][0] + c_exi1[child_elec->m][0];

}

void interpolatecornersnew(double t1)
{
	double tau = t1;
	int basex, basey,lx,ly;
	double distx, disty;
	
	if(child_elec->locx>0)
	{
		basex = child_elec->locx-1;
		basey = child_elec->locy;
		distx =	0.5+ 0.5/factor;
    	c_exs[0][0] = (1-tau)*(interpolate1d(exs_old[basex][basey],
			exs_old[basex+1][basey],distx))+(tau)*(interpolate1d(
			exs[basex][basey],exs[basex+1][basey],distx));
    					
    	c_vx[0][0] = (1-tau)*(interpolate1d(vx_old[basex][basey],
		     vx_old[basex+1][basey],distx))+(tau)*(interpolate1d(
		     vx[basex][basey],vx[basex+1][basey],distx));
	}

	if(child_elec->locy>0)
	{
		basex = child_elec->locx;
		basey = child_elec->locy-1;
		disty = 0.5+0.5/factor;
		c_eys[0][0] = (1-tau)*interpolate1d(eys_old[basex][basey],
			       eys_old[basex][basey+1],disty) + 
			      (tau)*interpolate1d(eys[basex][basey],
			       eys[basex][basey+1],disty);  
		c_vy[0][0] = (1-tau)*interpolate1d(vy_old[basex][basey],
			     vy_old[basex][basey+1],disty) + (tau)*interpolate1d(			      vy[basex][basey],vy[basex][basey+1],disty);
	}

	if(child_elec->locx>0&&child_elec->locy>0)
	{
		basex = child_elec->locx-1;
		basey = child_elec->locy-1;
		distx = 0.5 + 0.5/factor;
		disty = 0.5 + 0.5/factor;
		child_mag->mesh[0][0] = (1-tau)*(interpolate2d(hzi[basex][basey],					  hzi[basex+1][basey],
					  hzi[basex+1][basey+1],
					 hzi[basex][basey+1],distx,disty))+
					(tau)*(interpolate2d(
					 root_mag->mesh[basex][basey],
					root_mag->mesh[basex+1][basey],
					root_mag->mesh[basex+1][basey+1],
					root_mag->mesh[basex][basey+1],
					distx,disty));
	
  }

  i=(child_den->m);
  j=(child_den->n);
	child_den->mesh[0][0] = (1-tau)*denp[child_den->locx][child_den->locy] + tau*root_den->mesh[child_den->locx][child_den->locy];
  
  child_den->mesh[0][j] = (1-tau)*denp[child_den->locx][child_den->locy+j/factor] + tau*root_den->mesh[child_den->locx][child_den->locy+j/factor];
  
  child_den->mesh[i][j] = (1-tau)*denp[child_den->locx+i/factor][child_den->locy+j/factor] + tau*root_den->mesh[child_den->locx+i/factor][child_den->locy+j/factor];
  
  child_den->mesh[i][0] = (1-tau)*denp[child_den->locx+i/factor][child_den->locy] + tau*root_den->mesh[child_den->locx+i/factor][child_den->locy];
  
 
	basex = child_elec->locx+(child_elec->m)/factor;
	basey = child_elec->locy-1;
	disty = 0.5 + 0.5/factor;
	c_eys[child_elec->m][0] = (1-tau)*interpolate1d(eys_old[basex][basey],
				  eys_old[basex][basey+1],disty) + 
				  tau*interpolate1d(eys[basex][basey],
				  eys[basex][basey+1],disty);
  

	basex = child_elec->locx-1;
	basey = child_elec->locy+(child_elec->n)/factor;
	distx = 0.5 + 0.5/factor;
	c_exs[0][child_elec->n] = (1-tau)*interpolate1d(exs_old[basex][basey],
				exs_old[basex+1][basey],distx) + 
				tau*interpolate1d(exs[basex][basey],
				exs[basex+1][basey],distx);

	c_eyt[0][0] = c_eys[0][0] + c_eyi1[0][0];
	c_ext[0][0] = c_exs[0][0] + c_exi1[0][0];

c_ext[0][child_elec->n]= c_exs[0][child_elec->n] + c_exi1[0][child_elec->n];
c_eyt[child_elec->m][0] = c_eys[child_elec->m][0] + c_eyi1[child_elec->m][0];

}

void interpolatex(double t1)
{
    gettimeofday(&begin,NULL);
    double tau = t1;
    int basex, basey;
    double distx, disty;
    //int factor=8;
	// printf("interpolatex\n");
/*********************interpolation******************************/
for(j=1;j<child_elec->n;j++)
{

    double ita,ep;
    

    	if(child_elec->locx>0)
    	{
    		basex = child_elec->locx-1;
    		basey = child_elec->locy+j/factor;
		distx = 0.5 + 0.5/factor;
		disty = (double)(j%factor)/factor;
    		
	    	c_exs[0][j] = (1-tau)*(interpolate2d(exs_old[basex][basey],
			      exs_old[basex+1][basey],exs_old[basex+1][basey+1],
			      exs_old[basex][basey+1],distx,disty))+(tau)*(
			      interpolate2d(exs[basex][basey],exs[basex+1][basey]			       ,exs[basex+1][basey+1],exs[basex][basey+1],distx,
			      disty));
	    					
	    	c_vx[0][j] = (1-tau)*(interpolate2d(vx_old[basex][basey],
			     vx_old[basex+1][basey],vx_old[basex+1][basey+1],
			     vx_old[basex][basey+1],distx,disty))+(tau)*(
			     interpolate2d(vx[basex][basey],vx[basex+1][basey],
			     vx[basex+1][basey+1],vx[basex][basey+1],distx,
			     disty));
    	}

    	if(child_elec->locx+(child_elec->m-1)/factor+1<root_elec->m)
    	{
			basex = child_elec->locx+(child_elec->m-1)/factor;
			basey = child_elec->locy+j/factor;
			distx = 0.5 - 0.5/factor;
			disty = (double)(j%factor)/factor;
			c_exs[child_elec->m-1][j] = (1-tau)*(interpolate2d(
						    exs_old[basex][basey],
						    exs_old[basex+1][basey],
						    exs_old[basex+1][basey+1],
						    exs_old[basex][basey+1],
						    distx,disty))+(tau)*(
						    interpolate2d(
						    exs[basex][basey],
						    exs[basex+1][basey],
					            exs[basex+1][basey+1],
						    exs[basex][basey+1],distx,
						    disty));
		
			c_vx[child_elec->m-1][j] = (1-tau)*(interpolate2d(
						   vx_old[basex][basey],
						   vx_old[basex+1][basey],
						   vx_old[basex+1][basey+1],
						   vx_old[basex][basey+1],distx,
						   disty))+(tau)*(
						   interpolate2d(vx[basex][basey]						    ,vx[basex+1][basey],
						   vx[basex+1][basey+1],
						   vx[basex][basey+1],distx,
						   disty));
		}
    	
	if(j%factor<factor/2)
     {	
		basex = child_elec->locx;
		basey = child_elec->locy+j/factor-1;
		disty = 0.5 + (double)(j%factor)/(factor) + 0.5/(factor);

		c_eys[0][j] = (1-tau)*interpolate1d(eys_old[basex][basey],
			      eys_old[basex][basey+1],disty) + 
			      (tau)*interpolate1d(eys[basex][basey],
			       eys[basex][basey+1],disty);
		c_vy[0][j] = (1-tau)*interpolate1d(vy_old[basex][basey],
			     vy_old[basex][basey+1],disty) + 
			     (tau)*interpolate1d(vy[basex][basey],
			     vy[basex][basey+1],disty);
     }
	else
    {   
 	basex = child_elec->locx;
    	basey = child_elec->locy+j/factor;
	disty = (double)(j%factor)/(factor) + 0.5/(factor)-0.5;

    	c_eys[0][j] = (1-tau)*interpolate1d(eys_old[basex][basey],
		      eys_old[basex][basey+1],disty) + (tau)*interpolate1d(
                      eys[basex][basey],eys[basex][basey+1],disty);
    	c_vy[0][j] = (1-tau)*interpolate1d(vy_old[basex][basey],
		     vy_old[basex][basey+1],disty) + (tau)*interpolate1d(
		     vy[basex][basey],vy[basex][basey+1],disty);
    	
    }	
    
    if (j%factor<factor/2)
    {
		basex = child_elec->locx + (child_elec->m)/factor;
		basey = child_elec->locy+j/factor-1;
		disty = 0.5 + (double)(j%factor)/(factor) + 0.5/(factor);
		
	c_eys[child_elec->m][j] = (1-tau)*(interpolate1d(eys_old[basex][basey],
				  eys_old[basex][basey+1],disty))+(tau)*(
				  interpolate1d(eys[basex][basey],
				  eys[basex][basey+1],disty));

	c_vy[child_elec->m][j] = (1-tau)*(interpolate1d(vy_old[basex][basey],
				vy_old[basex][basey+1],disty))+(tau)*(
				interpolate1d(vy[basex][basey],
				vy[basex][basey+1],disty));
   }
  else
  {
	basex = child_elec->locx + (child_elec->m)/factor;
    	basey = child_elec->locy+j/factor;
	disty = (double)(j%factor)/(factor) + 0.5/(factor)-0.5;

    	c_eys[child_elec->m][j] = (1-tau)*(interpolate1d(eys_old[basex][basey],
				eys_old[basex][basey+1],disty))+(tau)*(
				interpolate1d(eys[basex][basey],
				eys[basex][basey+1],disty));
	
	c_vy[child_elec->m][j] = (1-tau)*(interpolate1d(vy_old[basex][basey],
				vy_old[basex][basey+1],disty))+(tau)*(
				interpolate1d(vy[basex][basey],
				vy[basex][basey+1],disty));

  }

    if(j%factor<factor/2 && child_elec->locx>0)
    {
	basex = child_elec->locx-1;
	basey = child_elec->locy+j/factor-1;
	distx = 0.5 + 0.5/factor;
	disty = 0.5 + 0.5/factor + (double)(j%factor)/factor;
	
	child_mag->mesh[0][j] = (1-tau)*(interpolate2d(hzi[basex][basey],
				hzi[basex+1][basey],hzi[basex+1][basey+1],
				hzi[basex][basey+1],distx,disty))+(tau)*(
				interpolate2d(root_mag->mesh[basex][basey],
				root_mag->mesh[basex+1][basey],
				root_mag->mesh[basex+1][basey+1],
				root_mag->mesh[basex][basey+1],distx,disty));
    }

    else if (child_elec->locx>0)
    {
		basex = child_elec->locx-1;
		basey = child_elec->locy+j/factor;
		distx = 0.5 + 0.5/factor;
		disty = 0.5/factor + (double)(j%factor)/factor-0.5;

	child_mag->mesh[0][j] = (1-tau)*(interpolate2d(hzi[basex][basey],
				hzi[basex+1][basey],hzi[basex+1][basey+1],
				hzi[basex][basey+1],distx,disty))+(tau)*(
				interpolate2d(root_mag->mesh[basex][basey],
				root_mag->mesh[basex+1][basey],
				root_mag->mesh[basex+1][basey+1],
				root_mag->mesh[basex][basey+1],distx,disty));
    } 
	
	if(j%factor<factor/2)
        {
                basex = child_elec->locx+(child_elec->m-1)/factor;
                basey = child_elec->locy+j/factor-1;
                distx = 0.5 - 0.5/factor;
                disty = 0.5 + (double)(j%factor)/factor + 0.5/factor;

         child_mag->mesh[child_elec->m-1][j] = (1-tau)*(interpolate2d(
						hzi[basex][basey],
						hzi[basex+1][basey],
						hzi[basex+1][basey+1],
						hzi[basex][basey+1],distx,disty))                                                +(tau)*(interpolate2d(
						root_mag->mesh[basex][basey],
						root_mag->mesh[basex+1][basey],
						root_mag->mesh[basex+1][basey+1],						root_mag->mesh[basex][basey+1],
						distx,disty));
        }

	else
        {
                basex = child_elec->locx+(child_elec->m-1)/factor;
                basey = child_elec->locy+j/factor;
                distx = 0.5 - 0.5/factor;
                disty = (double)(j%factor)/factor + 0.5/factor - 0.5;

         child_mag->mesh[child_elec->m-1][j] = (1-tau)*(interpolate2d(
						hzi[basex][basey],
						hzi[basex+1][basey],
						hzi[basex+1][basey+1],
						hzi[basex][basey+1],distx,disty))
					       +(tau)*(interpolate2d(
						root_mag->mesh[basex][basey],
						root_mag->mesh[basex+1][basey],
						root_mag->mesh[basex+1][basey+1],						 root_mag->mesh[basex][basey+1],
						distx,disty));
        }

	basex = child_elec->locx;
        basey = child_elec->locy+j/factor;
        disty = (double)(j%factor)/factor;
 
        child_den->mesh[0][j] = (1-tau)*(interpolate1d(denp[basex][basey],
				denp[basex][basey+1],disty))+(tau)*(
				interpolate1d(root_den->mesh[basex][basey],
				root_den->mesh[basex][basey+1],disty));

        basex = child_elec->locx+(child_elec->m)/factor;
        basey = child_elec->locy+j/factor;
        disty = (double)(j%factor)/factor;
 
       child_den->mesh[child_elec->m][j] =(1-tau)*(interpolate1d(
					  denp[basex][basey],
					  denp[basex][basey+1],disty))+(tau)*(
					  interpolate1d(
					  root_den->mesh[basex][basey],
					  root_den->mesh[basex][basey+1],
					  disty));

}
    gettimeofday(&end,NULL);
    t_cal_interpolatex += ((end.tv_sec - begin.tv_sec) + ((end.tv_usec - begin.tv_usec)/1000000.0));
}

void interpolatexinitial(double t1)
{
    gettimeofday(&begin,NULL);
    double tau = t1;
    int basex, basey,l,g;
    double distx, disty;
    
/*********************interpolation******************************/
for(i=0;i<=child_elec->m;i++)
{
  for(j=0;j<=child_elec->n;j++)
  {

    double ita,ep;
    double fac=0.25;
       
              //!=========== e-density ============
    
              basex = child_elec->locx+i/factor;
              basey = child_elec->locy+j/factor;
              distx = (double)(i%factor)/factor;
              disty = (double)(j%factor)/factor;
            
                        
              child_den->mesh[i][j] = (1-tau)*(interpolate2d(denp[basex][basey],
                    denp[basex+1][basey],denp[basex+1][basey+1],
                    denp[basex][basey+1],distx,disty))+(tau)*(
                    interpolate2d(root_den->mesh[basex][basey],
                    root_den->mesh[basex+1][basey],
                    root_den->mesh[basex+1][basey+1],
                    root_den->mesh[basex][basey+1],distx,disty));
                    
                    
                 //!=========== E-field + vel x ============
           
        if((child_elec->locx>0) && (i<child_elec->m))
        {
         if(i%factor<factor/2)
         {
          basex = child_elec->locx+i/factor-1;
          basey = child_elec->locy+j/factor;
          //distx = 0.5 + 0.5/factor;
          distx = 0.5 + 0.5/factor + (double)(i%factor)/factor;
          disty = (double)(j%factor)/factor;
          
          c_exs[i][j] = (1-tau)*(interpolate2d(exs_old[basex][basey],
              exs_old[basex+1][basey],exs_old[basex+1][basey+1],
              exs_old[basex][basey+1],distx,disty))+(tau)*(
              interpolate2d(exs[basex][basey],exs[basex+1][basey]            ,exs[basex+1][basey+1],exs[basex][basey+1],distx,
              disty));
                  
          c_vx[i][j] = (1-tau)*(interpolate2d(vx_old[basex][basey],
             vx_old[basex+1][basey],vx_old[basex+1][basey+1],
             vx_old[basex][basey+1],distx,disty))+(tau)*(
             interpolate2d(vx[basex][basey],vx[basex+1][basey],
             vx[basex+1][basey+1],vx[basex][basey+1],distx,
             disty));
        
               
        }
        else
        {   
          if(child_elec->locx+i/factor+1<(root_elec->m))
          {
          basex = child_elec->locx+i/factor;
          basey = child_elec->locy+j/factor;
          
          distx = 0.5/factor+ (double)(i%factor)/factor - 0.5;
          disty = (double)(j%factor)/factor;
          
          c_exs[i][j] = (1-tau)*(interpolate2d(exs_old[basex][basey],
              exs_old[basex+1][basey],exs_old[basex+1][basey+1],
              exs_old[basex][basey+1],distx,disty))+(tau)*(
              interpolate2d(exs[basex][basey],exs[basex+1][basey]            ,exs[basex+1][basey+1],exs[basex][basey+1],distx,
              disty));
                  
          c_vx[i][j] = (1-tau)*(interpolate2d(vx_old[basex][basey],
             vx_old[basex+1][basey],vx_old[basex+1][basey+1],
             vx_old[basex][basey+1],distx,disty))+(tau)*(
             interpolate2d(vx[basex][basey],vx[basex+1][basey],
             vx[basex+1][basey+1],vx[basex][basey+1],distx,
             disty));
          }
        }
        
       }
       
        //!=========== E-field + vel y ============
       
     if((child_elec->locy>0) && (j<child_elec->n))
     {
      if(j%factor<factor/2)
      {
                basex = child_elec->locx+(i)/factor;
                basey = child_elec->locy+j/factor-1;
                distx = (double)(i%factor)/factor;
                
                disty = 0.5 + 0.5/factor + (double)(j%factor)/factor;
                c_eys[i][j] = (1-tau)*(interpolate2d(eys_old[basex][basey],
                  eys_old[basex+1][basey],eys_old[basex+1][child_elec->locy],
                  eys_old[basex][child_elec->locy],distx,disty))+(tau)*(interpolate2d(eys[basex][basey],eys[basex+1][basey],
                  eys[basex+1][basey+1],eys[basex][basey+1],distx,disty));

                c_vy[i][j] = (1-tau)*(interpolate2d(vy_old[basex][basey],vy_old[basex+1][basey],
                  vy_old[basex+1][basey+1],vy_old[basex][basey+1],distx,disty))+(tau)*(interpolate2d(vy[basex][basey],
                  vy[basex+1][basey],vy[basex+1][basey+1],vy[basex][basey+1],distx,disty));
      }
      
      else
      {
          if(child_elec->locy+j/factor+1<(root_elec->n))
          {
                basex = child_elec->locx+(i)/factor;
                basey = child_elec->locy+j/factor;
                distx = (double)(i%factor)/factor;
                //disty = 0.5 + 0.5/factor;
                disty = 0.5/factor+ (double)(j%factor)/factor - 0.5;
                c_eys[i][j] = (1-tau)*(interpolate2d(eys_old[basex][basey],
                  eys_old[basex+1][basey],eys_old[basex+1][child_elec->locy],
                  eys_old[basex][child_elec->locy],distx,disty))+(tau)*(interpolate2d(eys[basex][basey],
                  eys[basex+1][basey],eys[basex+1][basey+1],eys[basex][basey+1],distx,disty));

                c_vy[i][j] = (1-tau)*(interpolate2d(vy_old[basex][basey],
                  vy_old[basex+1][basey],vy_old[basex+1][basey+1],
                  vy_old[basex][basey+1],distx,disty))+(tau)*(interpolate2d(vy[basex][basey],
                  vy[basex+1][basey],vy[basex+1][basey+1],vy[basex][basey+1],distx,disty));
      
         }
       }
           
      }        
                    
        //!====================== H-field ================
       
       if((child_elec->locy>0) && (j<child_elec->n))
          {
            if(j%factor<factor/2)
            {
              
              if (i%factor<factor/2)
              {  
        
        
                basex = child_elec->locx+(i)/factor-1;
                basey = child_elec->locy-1+j/factor;
                distx = 0.5 + 0.5/factor + (double)(i%factor)/factor;
                disty = 0.5 + 0.5/factor + (double)(j%factor)/factor;
                child_mag->mesh[i][j] = (1-tau)*(interpolate2d(hzi[basex][basey],hzi[basex+1][basey],
                  hzi[basex+1][basey+1],hzi[basex][basey+1],distx,disty))+(tau)*(interpolate2d(root_mag->mesh[basex][basey],
                  root_mag->mesh[basex+1][basey],root_mag->mesh[basex+1][basey+1],root_mag->mesh[basex][basey+1],distx,disty));
                                                                
              }
            
            else
            {
                
              basex = child_elec->locx+i/factor;
              basey = child_elec->locy+j/factor-1;
              distx = 0.5/factor + (double)(i%factor)/factor - 0.5;
              disty = 0.5 + 0.5/factor + (double)(j%factor)/factor;
              child_mag->mesh[i][j] = (1-tau)*(interpolate2d(hzi[basex][basey],hzi[basex+1][basey],
                hzi[basex+1][basey+1],hzi[basex][basey+1],distx,disty))+(tau)*(interpolate2d(root_mag->mesh[basex][basey],
                root_mag->mesh[basex+1][basey],root_mag->mesh[basex+1][basey+1],root_mag->mesh[basex][basey+1],distx,disty));
                
             }  
            
            }    
                
       else    
        {
            if (i%factor<factor/2)
             {  
                basex = child_elec->locx+(i)/factor-1;
                basey = child_elec->locy+j/factor;
                distx = 0.5 + 0.5/factor + (double)(i%factor)/factor; 
                disty = 0.5/factor+ (double)(j%factor)/factor - 0.5; 
                child_mag->mesh[i][j] = (1-tau)*(interpolate2d(hzi[basex][basey],hzi[basex+1][basey],
                  hzi[basex+1][basey+1],hzi[basex][basey+1],distx,disty))+(tau)*(interpolate2d(root_mag->mesh[basex][basey],
                  root_mag->mesh[basex+1][basey],root_mag->mesh[basex+1][basey+1],root_mag->mesh[basex][basey+1],distx,disty));
            
              }
            else
            {
              basex = child_elec->locx+i/factor;
              basey = child_elec->locy+j/factor;
              distx = 0.5/factor + (double)(i%factor)/factor - 0.5;
              disty = 0.5/factor+ (double)(j%factor)/factor - 0.5; 
               child_mag->mesh[i][j] = (1-tau)*(interpolate2d(hzi[basex][basey],hzi[basex+1][basey],
                 hzi[basex+1][basey+1],hzi[basex][basey+1],distx,disty))+(tau)*(interpolate2d(root_mag->mesh[basex][basey],
                 root_mag->mesh[basex+1][basey],root_mag->mesh[basex+1][basey+1],root_mag->mesh[basex][basey+1],distx,disty));
            }
        
          }        
      }       
      
      //!========================================================             
                    
       
    }
}
    gettimeofday(&end,NULL);
    t_cal_interpolatexinitial += ((end.tv_sec - begin.tv_sec) + ((end.tv_usec - begin.tv_usec)/1000000.0));                 
                    

}

void interpolatexnew(double t1)
{
    gettimeofday(&begin,NULL);
    double tau = t1;
    int basex, basey,l,g;
    double distx, disty;
    
/*********************interpolation along x(expanding box)******************************/
for(i=0;i<=child_elec->m;i++)
{
  for(j=0;j<=child_elec->n;j++)
  {

    double ita,ep;
    double fac=0.25;
       
    
              basex = child_elec->locx+i/factor;
            	basey = child_elec->locy+j/factor;
              distx = (double)(i%factor)/factor;
              disty = (double)(j%factor)/factor;
            
                      	
            	child_den->mesh[i][j] = (1-tau)*(interpolate2d(denp[basex][basey],
            				denp[basex+1][basey],denp[basex+1][basey+1],
            				denp[basex][basey+1],distx,disty))+(tau)*(
            				interpolate2d(root_den->mesh[basex][basey],
            				root_den->mesh[basex+1][basey],
            				root_den->mesh[basex+1][basey+1],
            				root_den->mesh[basex][basey+1],distx,disty));
           
          
            if((child_den->locx+i/factor)>=xstarol && i<=child_elec->m)
               {
                     	child_den->mesh[i][j] = c_denpold[i-(int)((xstarol-xstar)*factor)][j];
                     
               }
            
                   
             if((child_den->locx+i/factor)==xstarol)
             {
                    basex = child_elec->locx+i/factor;
                  	basey = child_elec->locy+j/factor;
                    distx = (double)(i%factor)/factor;
                    disty = (double)(j%factor)/factor;
             
                   	child_den->mesh[i][j] = (1-fac)*(c_denpold[i-(int)((xstarol-xstar)*factor)][j])+(fac)*(
                    interpolate2d(root_den->mesh[basex][basey],
            				root_den->mesh[basex+1][basey],
            				root_den->mesh[basex+1][basey+1],
            				root_den->mesh[basex][basey+1],distx,disty));
                  
             }
           
            //!=========== E-field + vel x ============
           
        if((child_elec->locx>0) && (i<child_elec->m))
        {
         if(i%factor<factor/2)
         {
          basex = child_elec->locx+i/factor-1;
      		basey = child_elec->locy+j/factor;
      		
          distx = 0.5 + 0.5/factor + (double)(i%factor)/factor;
      		disty = (double)(j%factor)/factor;
      		
  	    	c_exs[i][j] = (1-tau)*(interpolate2d(exs_old[basex][basey],
  			      exs_old[basex+1][basey],exs_old[basex+1][basey+1],
  			      exs_old[basex][basey+1],distx,disty))+(tau)*(
  			      interpolate2d(exs[basex][basey],exs[basex+1][basey]			       ,exs[basex+1][basey+1],exs[basex][basey+1],distx,
  			      disty));
  	    					
  	    	c_vx[i][j] = (1-tau)*(interpolate2d(vx_old[basex][basey],
  			     vx_old[basex+1][basey],vx_old[basex+1][basey+1],
  			     vx_old[basex][basey+1],distx,disty))+(tau)*(
  			     interpolate2d(vx[basex][basey],vx[basex+1][basey],
  			     vx[basex+1][basey+1],vx[basex][basey+1],distx,
  			     disty));
        
               
        }
        else
        {   
          if(child_elec->locx+i/factor+1<(root_elec->m))
          {
          basex = child_elec->locx+i/factor;
      		basey = child_elec->locy+j/factor;
      		
          distx = 0.5/factor+ (double)(i%factor)/factor - 0.5;
      		disty = (double)(j%factor)/factor;
      		
  	    	c_exs[i][j] = (1-tau)*(interpolate2d(exs_old[basex][basey],
  			      exs_old[basex+1][basey],exs_old[basex+1][basey+1],
  			      exs_old[basex][basey+1],distx,disty))+(tau)*(
  			      interpolate2d(exs[basex][basey],exs[basex+1][basey]			       ,exs[basex+1][basey+1],exs[basex][basey+1],distx,
  			      disty));
  	    					
  	    	c_vx[i][j] = (1-tau)*(interpolate2d(vx_old[basex][basey],
  			     vx_old[basex+1][basey],vx_old[basex+1][basey+1],
  			     vx_old[basex][basey+1],distx,disty))+(tau)*(
  			     interpolate2d(vx[basex][basey],vx[basex+1][basey],
  			     vx[basex+1][basey+1],vx[basex][basey+1],distx,
  			     disty));
          }
        }
        
    	 }
          
         if((child_elec->locx+i/factor)>=xstarol && i<child_elec->m)
               {
                     	c_exs[i][j] = c_exold[i-(int)((xstarol-xstar)*factor)][j];
                      
                      c_vx[i][j] = c_vxold[i-(int)((xstarol-xstar)*factor)][j];
               } 
          
             
       //!=========== E-field + vel y ============
       
     if((child_elec->locy>0) && (j<child_elec->n))
     {
      if(j%factor<factor/2)
      {
		            basex = child_elec->locx+(i)/factor;
                basey = child_elec->locy+j/factor-1;
                distx = (double)(i%factor)/factor;
                
                disty = 0.5 + 0.5/factor + (double)(j%factor)/factor;
                c_eys[i][j] = (1-tau)*(interpolate2d(eys_old[basex][basey],
                  eys_old[basex+1][basey],eys_old[basex+1][child_elec->locy],
                  eys_old[basex][child_elec->locy],distx,disty))+(tau)*(interpolate2d(eys[basex][basey],eys[basex+1][basey],
                  eys[basex+1][basey+1],eys[basex][basey+1],distx,disty));

                c_vy[i][j] = (1-tau)*(interpolate2d(vy_old[basex][basey],vy_old[basex+1][basey],
                  vy_old[basex+1][basey+1],vy_old[basex][basey+1],distx,disty))+(tau)*(interpolate2d(vy[basex][basey],
                  vy[basex+1][basey],vy[basex+1][basey+1],vy[basex][basey+1],distx,disty));
      }
      
      else
      {
          if(child_elec->locy+j/factor+1<(root_elec->n))
          {
		            basex = child_elec->locx+(i)/factor;
                basey = child_elec->locy+j/factor;
                distx = (double)(i%factor)/factor;
                
                disty = 0.5/factor+ (double)(j%factor)/factor - 0.5;
                c_eys[i][j] = (1-tau)*(interpolate2d(eys_old[basex][basey],
                  eys_old[basex+1][basey],eys_old[basex+1][child_elec->locy],
                  eys_old[basex][child_elec->locy],distx,disty))+(tau)*(interpolate2d(eys[basex][basey],
                  eys[basex+1][basey],eys[basex+1][basey+1],eys[basex][basey+1],distx,disty));

                c_vy[i][j] = (1-tau)*(interpolate2d(vy_old[basex][basey],
                  vy_old[basex+1][basey],vy_old[basex+1][basey+1],
                  vy_old[basex][basey+1],distx,disty))+(tau)*(interpolate2d(vy[basex][basey],
                  vy[basex+1][basey],vy[basex+1][basey+1],vy[basex][basey+1],distx,disty));
      
         }
       }
           
      }
         
          if((child_elec->locx+i/factor)>=xstarol && i<child_elec->m)
               {
                     	c_eys[i][j] = c_eyold[i-(int)((xstarol-xstar)*factor)][j];
                      
                      c_vy[i][j] = c_vyold[i-(int)((xstarol-xstar)*factor)][j];
               } 
        
     
       //!====================== H-field ================
       
       if((child_elec->locy>0) && (j<child_elec->n))
          {
            if(j%factor<factor/2)
            {
              
              if (i%factor<factor/2)
              {  
        
        
                basex = child_elec->locx+(i)/factor-1;
                basey = child_elec->locy-1+j/factor;
                distx = 0.5 + 0.5/factor + (double)(i%factor)/factor;
                disty = 0.5 + 0.5/factor + (double)(j%factor)/factor;
                child_mag->mesh[i][j] = (1-tau)*(interpolate2d(hzi[basex][basey],hzi[basex+1][basey],
                  hzi[basex+1][basey+1],hzi[basex][basey+1],distx,disty))+(tau)*(interpolate2d(root_mag->mesh[basex][basey],
                  root_mag->mesh[basex+1][basey],root_mag->mesh[basex+1][basey+1],root_mag->mesh[basex][basey+1],distx,disty));
                                                                
              }
            
            else
            {
                
              basex = child_elec->locx+i/factor;
              basey = child_elec->locy+j/factor-1;
              distx = 0.5/factor + (double)(i%factor)/factor - 0.5;
              disty = 0.5 + 0.5/factor + (double)(j%factor)/factor;
              child_mag->mesh[i][j] = (1-tau)*(interpolate2d(hzi[basex][basey],hzi[basex+1][basey],
                hzi[basex+1][basey+1],hzi[basex][basey+1],distx,disty))+(tau)*(interpolate2d(root_mag->mesh[basex][basey],
                root_mag->mesh[basex+1][basey],root_mag->mesh[basex+1][basey+1],root_mag->mesh[basex][basey+1],distx,disty));
                
             }  
            
            }    
                
       else    
        {
            if (i%factor<factor/2)
             {  
                basex = child_elec->locx+(i)/factor-1;
                basey = child_elec->locy+j/factor;
                distx = 0.5 + 0.5/factor + (double)(i%factor)/factor; 
                disty = 0.5/factor+ (double)(j%factor)/factor - 0.5; 
                child_mag->mesh[i][j] = (1-tau)*(interpolate2d(hzi[basex][basey],hzi[basex+1][basey],
                  hzi[basex+1][basey+1],hzi[basex][basey+1],distx,disty))+(tau)*(interpolate2d(root_mag->mesh[basex][basey],
                  root_mag->mesh[basex+1][basey],root_mag->mesh[basex+1][basey+1],root_mag->mesh[basex][basey+1],distx,disty));
            
              }
            else
            {
              basex = child_elec->locx+i/factor;
              basey = child_elec->locy+j/factor;
              distx = 0.5/factor + (double)(i%factor)/factor - 0.5;
              disty = 0.5/factor+ (double)(j%factor)/factor - 0.5; 
               child_mag->mesh[i][j] = (1-tau)*(interpolate2d(hzi[basex][basey],hzi[basex+1][basey],
                 hzi[basex+1][basey+1],hzi[basex][basey+1],distx,disty))+(tau)*(interpolate2d(root_mag->mesh[basex][basey],
                 root_mag->mesh[basex+1][basey],root_mag->mesh[basex+1][basey+1],root_mag->mesh[basex][basey+1],distx,disty));
            }
        
          }        
      }          
         
         if((child_elec->locx+i/factor)>=xstarol && i<child_elec->m)
               {
                     	child_mag->mesh[i][j] = c_hzold[i-(int)((xstarol-xstar)*factor)][j];
                      
               } 
   
       
      //!===================================     
         
      
    
    }
}
    gettimeofday(&end,NULL);
    t_cal_interpolatexnew += ((end.tv_sec - begin.tv_sec) + ((end.tv_usec - begin.tv_usec)/1000000.0));
}

void interpolateynew(double t1)
{
    gettimeofday(&begin,NULL);
    double tau = t1;
    int basex, basey,l,g;
    double distx, disty;
    
/*********************interpolation along y(expanding box)******************************/
for(i=0;i<=child_elec->m;i++)
{
  for(j=0;j<=child_elec->n;j++)
  {
    double ita,ep;
    double fac=0.25;
       
    
              basex = child_elec->locx+i/factor;
            	basey = child_elec->locy+j/factor;
              distx = (double)(i%factor)/factor;
              disty = (double)(j%factor)/factor;
            
                      	
            	child_den->mesh[i][j] = (1-tau)*(interpolate2d(denp[basex][basey],
            				denp[basex+1][basey],denp[basex+1][basey+1],
            				denp[basex][basey+1],distx,disty))+(tau)*(
            				interpolate2d(root_den->mesh[basex][basey],
            				root_den->mesh[basex+1][basey],
            				root_den->mesh[basex+1][basey+1],
            				root_den->mesh[basex][basey+1],distx,disty));
           
          
              if((child_den->locy+j/factor)>=(ystarol) && (child_den->locy+j/factor)<=(yendol) )
              
               {
                        child_den->mesh[i][j] = c_denpold[i][j-(int)((ystarol-ystar)*factor)];
                    
               }
            
                   
            
             if((child_den->locy+j/factor)==ystarol || (child_den->locy+j/factor)==(yendol))
             {
                    basex = child_elec->locx+i/factor;
                  	basey = child_elec->locy+j/factor;
                    distx = (double)(i%factor)/factor;
                    disty = (double)(j%factor)/factor;
             
                   child_den->mesh[i][j] = (1-fac)*(c_denpold[i][j-(int)((ystarol-ystar)*factor)])+(fac)*(
                    interpolate2d(root_den->mesh[basex][basey],
            				root_den->mesh[basex+1][basey],
            				root_den->mesh[basex+1][basey+1],
            				root_den->mesh[basex][basey+1],distx,disty));
                 
             }
           
            //!=========== E-field + vel x ============
           
        if((child_elec->locx>0) && (i<child_elec->m))
        {
         if(i%factor<factor/2)
         {
          basex = child_elec->locx+i/factor-1;
      		basey = child_elec->locy+j/factor;
      		
          distx = 0.5 + 0.5/factor + (double)(i%factor)/factor;
      		disty = (double)(j%factor)/factor;
      		
  	    	c_exs[i][j] = (1-tau)*(interpolate2d(exs_old[basex][basey],
  			      exs_old[basex+1][basey],exs_old[basex+1][basey+1],
  			      exs_old[basex][basey+1],distx,disty))+(tau)*(
  			      interpolate2d(exs[basex][basey],exs[basex+1][basey]			       ,exs[basex+1][basey+1],exs[basex][basey+1],distx,
  			      disty));
  	    					
  	    	c_vx[i][j] = (1-tau)*(interpolate2d(vx_old[basex][basey],
  			     vx_old[basex+1][basey],vx_old[basex+1][basey+1],
  			     vx_old[basex][basey+1],distx,disty))+(tau)*(
  			     interpolate2d(vx[basex][basey],vx[basex+1][basey],
  			     vx[basex+1][basey+1],vx[basex][basey+1],distx,
  			     disty));
        
               
        }
        else
        {   
          if(child_elec->locx+i/factor+1<(root_elec->m))
          {
          basex = child_elec->locx+i/factor;
      		basey = child_elec->locy+j/factor;
      		
          distx = 0.5/factor+ (double)(i%factor)/factor - 0.5;
      		disty = (double)(j%factor)/factor;
      		
  	    	c_exs[i][j] = (1-tau)*(interpolate2d(exs_old[basex][basey],
  			      exs_old[basex+1][basey],exs_old[basex+1][basey+1],
  			      exs_old[basex][basey+1],distx,disty))+(tau)*(
  			      interpolate2d(exs[basex][basey],exs[basex+1][basey]			       ,exs[basex+1][basey+1],exs[basex][basey+1],distx,
  			      disty));
  	    					
  	    	c_vx[i][j] = (1-tau)*(interpolate2d(vx_old[basex][basey],
  			     vx_old[basex+1][basey],vx_old[basex+1][basey+1],
  			     vx_old[basex][basey+1],distx,disty))+(tau)*(
  			     interpolate2d(vx[basex][basey],vx[basex+1][basey],
  			     vx[basex+1][basey+1],vx[basex][basey+1],distx,
  			     disty));
          }
        }
        
    	 }
          
        
              if((child_den->locy+j/factor)>=(ystarol) && (child_den->locy+j/factor)<=(yendol) )
               {
                     	c_exs[i][j] = c_exold[i][j-(int)((ystarol-ystar)*factor)];
                      
                      c_vx[i][j] = c_vxold[i][j-(int)((ystarol-ystar)*factor)];
               } 
          
             
       //!=========== E-field + vel y ============
       
     if((child_elec->locy>0) && (j<child_elec->n))
     {
      if(j%factor<factor/2)
      {
		            basex = child_elec->locx+(i)/factor;
                basey = child_elec->locy+j/factor-1;
                distx = (double)(i%factor)/factor;
                
                disty = 0.5 + 0.5/factor + (double)(j%factor)/factor;
                c_eys[i][j] = (1-tau)*(interpolate2d(eys_old[basex][basey],
                  eys_old[basex+1][basey],eys_old[basex+1][child_elec->locy],
                  eys_old[basex][child_elec->locy],distx,disty))+(tau)*(interpolate2d(eys[basex][basey],eys[basex+1][basey],
                  eys[basex+1][basey+1],eys[basex][basey+1],distx,disty));

                c_vy[i][j] = (1-tau)*(interpolate2d(vy_old[basex][basey],vy_old[basex+1][basey],
                  vy_old[basex+1][basey+1],vy_old[basex][basey+1],distx,disty))+(tau)*(interpolate2d(vy[basex][basey],
                  vy[basex+1][basey],vy[basex+1][basey+1],vy[basex][basey+1],distx,disty));
      }
      
      else
      {
          if(child_elec->locy+j/factor+1<(root_elec->n))
          {
		            basex = child_elec->locx+(i)/factor;
                basey = child_elec->locy+j/factor;
                distx = (double)(i%factor)/factor;
                
                disty = 0.5/factor+ (double)(j%factor)/factor - 0.5;
                c_eys[i][j] = (1-tau)*(interpolate2d(eys_old[basex][basey],
                  eys_old[basex+1][basey],eys_old[basex+1][child_elec->locy],
                  eys_old[basex][child_elec->locy],distx,disty))+(tau)*(interpolate2d(eys[basex][basey],
                  eys[basex+1][basey],eys[basex+1][basey+1],eys[basex][basey+1],distx,disty));

                c_vy[i][j] = (1-tau)*(interpolate2d(vy_old[basex][basey],
                  vy_old[basex+1][basey],vy_old[basex+1][basey+1],
                  vy_old[basex][basey+1],distx,disty))+(tau)*(interpolate2d(vy[basex][basey],
                  vy[basex+1][basey],vy[basex+1][basey+1],vy[basex][basey+1],distx,disty));
      
         }
       }
           
      }
         
         
               if((child_den->locy+j/factor)>=(ystarol) && (child_den->locy+j/factor)<=(yendol) )
               {
                     	c_eys[i][j] = c_eyold[i][j-(int)((ystarol-ystar)*factor)];
                      
                      c_vy[i][j] = c_vyold[i][j-(int)((ystarol-ystar)*factor)];
               } 
        
     
       //!====================== H-field ================
       
       if((child_elec->locy>0) && (j<child_elec->n))
          {
            if(j%factor<factor/2)
            {
              
              if (i%factor<factor/2)
              {  
        
        
                basex = child_elec->locx+(i)/factor-1;
                basey = child_elec->locy-1+j/factor;
                distx = 0.5 + 0.5/factor + (double)(i%factor)/factor;
                disty = 0.5 + 0.5/factor + (double)(j%factor)/factor;
                child_mag->mesh[i][j] = (1-tau)*(interpolate2d(hzi[basex][basey],hzi[basex+1][basey],
                  hzi[basex+1][basey+1],hzi[basex][basey+1],distx,disty))+(tau)*(interpolate2d(root_mag->mesh[basex][basey],
                  root_mag->mesh[basex+1][basey],root_mag->mesh[basex+1][basey+1],root_mag->mesh[basex][basey+1],distx,disty));
                                                                
              }
            
            else
            {
                
              basex = child_elec->locx+i/factor;
              basey = child_elec->locy+j/factor-1;
              distx = 0.5/factor + (double)(i%factor)/factor - 0.5;
              disty = 0.5 + 0.5/factor + (double)(j%factor)/factor;
              child_mag->mesh[i][j] = (1-tau)*(interpolate2d(hzi[basex][basey],hzi[basex+1][basey],
                hzi[basex+1][basey+1],hzi[basex][basey+1],distx,disty))+(tau)*(interpolate2d(root_mag->mesh[basex][basey],
                root_mag->mesh[basex+1][basey],root_mag->mesh[basex+1][basey+1],root_mag->mesh[basex][basey+1],distx,disty));
                
             }  
            
            }    
                
       else    
        {
            if (i%factor<factor/2)
             {  
                basex = child_elec->locx+(i)/factor-1;
                basey = child_elec->locy+j/factor;
                distx = 0.5 + 0.5/factor + (double)(i%factor)/factor; 
                disty = 0.5/factor+ (double)(j%factor)/factor - 0.5; 
                child_mag->mesh[i][j] = (1-tau)*(interpolate2d(hzi[basex][basey],hzi[basex+1][basey],
                  hzi[basex+1][basey+1],hzi[basex][basey+1],distx,disty))+(tau)*(interpolate2d(root_mag->mesh[basex][basey],
                  root_mag->mesh[basex+1][basey],root_mag->mesh[basex+1][basey+1],root_mag->mesh[basex][basey+1],distx,disty));
            
              }
            else
            {
              basex = child_elec->locx+i/factor;
              basey = child_elec->locy+j/factor;
              distx = 0.5/factor + (double)(i%factor)/factor - 0.5;
              disty = 0.5/factor+ (double)(j%factor)/factor - 0.5; 
               child_mag->mesh[i][j] = (1-tau)*(interpolate2d(hzi[basex][basey],hzi[basex+1][basey],
                 hzi[basex+1][basey+1],hzi[basex][basey+1],distx,disty))+(tau)*(interpolate2d(root_mag->mesh[basex][basey],
                 root_mag->mesh[basex+1][basey],root_mag->mesh[basex+1][basey+1],root_mag->mesh[basex][basey+1],distx,disty));
            }
        
          }        
      }          
         
        
               if((child_den->locy+j/factor)>=(ystarol) && (child_den->locy+j/factor)<=(yendol) ) 
               {
                     	child_mag->mesh[i][j] = c_hzold[i][j-(int)((ystarol-ystar)*factor)];
                      
               } 
   
       
      //!===================================     
       
     
   
       
    }
}
    gettimeofday(&end,NULL);
    t_cal_interpolateynew += ((end.tv_sec - begin.tv_sec) + ((end.tv_usec - begin.tv_usec)/1000000.0));
}

void interpolatexall(double t1)
{
    gettimeofday(&begin,NULL);
    double tau = t1;
    int basex, basey;
    double distx, disty;
    
/*********************interpolation******************************/
for(i=0;i<=child_elec->m;i++)
{
  for(j=1;j<child_elec->n;j++)
  {
//	int factor=2;
    double ita,ep;

 
    
    if((child_elec->locx>0) && (i<child_elec->m))
    {
       if(i%factor<factor/2)
      {
          basex = child_elec->locx+i/factor-1;
      		basey = child_elec->locy+j/factor;
      		
          distx = 0.5 + 0.5/factor + (double)(i%factor)/factor;
      		disty = (double)(j%factor)/factor;
      		
  	    	c_exs[i][j] = (1-tau)*(interpolate2d(exs_old[basex][basey],
  			      exs_old[basex+1][basey],exs_old[basex+1][basey+1],
  			      exs_old[basex][basey+1],distx,disty))+(tau)*(
  			      interpolate2d(exs[basex][basey],exs[basex+1][basey]			       ,exs[basex+1][basey+1],exs[basex][basey+1],distx,
  			      disty));
  	    					
  	    	c_vx[i][j] = (1-tau)*(interpolate2d(vx_old[basex][basey],
  			     vx_old[basex+1][basey],vx_old[basex+1][basey+1],
  			     vx_old[basex][basey+1],distx,disty))+(tau)*(
  			     interpolate2d(vx[basex][basey],vx[basex+1][basey],
  			     vx[basex+1][basey+1],vx[basex][basey+1],distx,
  			     disty));
        }
      else
      {   
          if(child_elec->locx+i/factor+1<(root_elec->m))
          {
          basex = child_elec->locx+i/factor;
      		basey = child_elec->locy+j/factor;
      		
          distx = 0.5/factor+ (double)(i%factor)/factor - 0.5;
      		disty = (double)(j%factor)/factor;
      		
  	    	c_exs[i][j] = (1-tau)*(interpolate2d(exs_old[basex][basey],
  			      exs_old[basex+1][basey],exs_old[basex+1][basey+1],
  			      exs_old[basex][basey+1],distx,disty))+(tau)*(
  			      interpolate2d(exs[basex][basey],exs[basex+1][basey]			       ,exs[basex+1][basey+1],exs[basex][basey+1],distx,
  			      disty));
  	    					
  	    	c_vx[i][j] = (1-tau)*(interpolate2d(vx_old[basex][basey],
  			     vx_old[basex+1][basey],vx_old[basex+1][basey+1],
  			     vx_old[basex][basey+1],distx,disty))+(tau)*(
  			     interpolate2d(vx[basex][basey],vx[basex+1][basey],
  			     vx[basex+1][basey+1],vx[basex][basey+1],distx,
  			     disty));
          }
      }
        
    	}
  
   
   
	if(j%factor<factor/2)
     {	
        
        basex = child_elec->locx;
    		basey = child_elec->locy+j/factor-1;
    		disty = 0.5 + (double)(j%factor)/(factor) + 0.5/(factor);
    
    	c_eys[0][j] = (1-tau)*interpolate1d(eys_old[basex][basey],
			      eys_old[basex][basey+1],disty) + 
			      (tau)*interpolate1d(eys[basex][basey],
			       eys[basex][basey+1],disty);
		  c_vy[0][j] = (1-tau)*interpolate1d(vy_old[basex][basey],
			     vy_old[basex][basey+1],disty) + 
			     (tau)*interpolate1d(vy[basex][basey],
			     vy[basex][basey+1],disty);
        
     }
	else
    {   
      
        basex = child_elec->locx;
    	  basey = child_elec->locy+j/factor;
	      disty = (double)(j%factor)/(factor) + 0.5/(factor)-0.5;

      	c_eys[0][j] = (1-tau)*interpolate1d(eys_old[basex][basey],
		      eys_old[basex][basey+1],disty) + (tau)*interpolate1d(
                      eys[basex][basey],eys[basex][basey+1],disty);
     	  c_vy[0][j] = (1-tau)*interpolate1d(vy_old[basex][basey],
		     vy_old[basex][basey+1],disty) + (tau)*interpolate1d(
		     vy[basex][basey],vy[basex][basey+1],disty);
              
        
    	
    }	
    
    if (j%factor<factor/2)
    {
      
		basex = child_elec->locx + (child_elec->m)/factor;
		basey = child_elec->locy+j/factor-1;
		disty = 0.5 + (double)(j%factor)/(factor) + 0.5/(factor);
		
  	c_eys[child_elec->m][j] = (1-tau)*(interpolate1d(eys_old[basex][basey],
  				  eys_old[basex][basey+1],disty))+(tau)*(
  				  interpolate1d(eys[basex][basey],
  				  eys[basex][basey+1],disty));
  
  	c_vy[child_elec->m][j] = (1-tau)*(interpolate1d(vy_old[basex][basey],
  				vy_old[basex][basey+1],disty))+(tau)*(
  				interpolate1d(vy[basex][basey],
  				vy[basex][basey+1],disty));
      
   }
  else
  {
  
     
      	basex = child_elec->locx + (child_elec->m)/factor;
        basey = child_elec->locy+j/factor;
      	disty = (double)(j%factor)/(factor) + 0.5/(factor)-0.5;

    		c_eys[child_elec->m][j] = (1-tau)*(interpolate1d(eys_old[basex][basey],
				eys_old[basex][basey+1],disty))+(tau)*(
				interpolate1d(eys[basex][basey],
				eys[basex][basey+1],disty));
	
      	c_vy[child_elec->m][j] = (1-tau)*(interpolate1d(vy_old[basex][basey],
      				vy_old[basex][basey+1],disty))+(tau)*(
      				interpolate1d(vy[basex][basey],
      				vy[basex][basey+1],disty));
       
      

  }

    if(j%factor<factor/2 && child_elec->locx>0)
    {
       
        	basex = child_elec->locx-1;
        	basey = child_elec->locy+j/factor-1;
        	distx = 0.5 + 0.5/factor;
        	disty = 0.5 + 0.5/factor + (double)(j%factor)/factor;
        	
          child_mag->mesh[0][j] = (1-tau)*(interpolate2d(hzi[basex][basey],
  				hzi[basex+1][basey],hzi[basex+1][basey+1],
  				hzi[basex][basey+1],distx,disty))+(tau)*(
  				interpolate2d(root_mag->mesh[basex][basey],
  				root_mag->mesh[basex+1][basey],
  				root_mag->mesh[basex+1][basey+1],
  				root_mag->mesh[basex][basey+1],distx,disty));
      
    }

    else if (child_elec->locx>0)
    {
     
      		basex = child_elec->locx-1;
      		basey = child_elec->locy+j/factor;
      		distx = 0.5 + 0.5/factor;
      		disty = 0.5/factor + (double)(j%factor)/factor-0.5;
      
      	child_mag->mesh[0][j] = (1-tau)*(interpolate2d(hzi[basex][basey],
				hzi[basex+1][basey],hzi[basex+1][basey+1],
				hzi[basex][basey+1],distx,disty))+(tau)*(
				interpolate2d(root_mag->mesh[basex][basey],
				root_mag->mesh[basex+1][basey],
				root_mag->mesh[basex+1][basey+1],
				root_mag->mesh[basex][basey+1],distx,disty));
     
    } 
	
   //!================ density ============
     
        	basex = child_elec->locx+i/factor;
        	basey = child_elec->locy+j/factor;
          distx = (double)(i%factor)/factor;
          disty = (double)(j%factor)/factor;
        
      
        	
        	child_den->mesh[i][j] = (1-tau)*(interpolate2d(denp[basex][basey],
        				denp[basex+1][basey],denp[basex+1][basey+1],
        				denp[basex][basey+1],distx,disty))+(tau)*(
        				interpolate2d(root_den->mesh[basex][basey],
        				root_den->mesh[basex+1][basey],
        				root_den->mesh[basex+1][basey+1],
        				root_den->mesh[basex][basey+1],distx,disty));
    
       
      //!===================================     
       
   
    
    }
}
    gettimeofday(&end,NULL);
    t_cal_interpolatex += ((end.tv_sec - begin.tv_sec) + ((end.tv_usec - begin.tv_usec)/1000000.0));
}

void interpolateyall(double t1)
{
     double tau = t1;
     int basex, basey;
     double distx,disty;
     //int factor=8;

for(j=0;j<=child_elec->n;j++ ) 
{
  for(i=1;i<child_elec->m;i++)
  {
	//int factor=2;
     double ita,ep;
  
    if((child_elec->locy>0) && (j<child_elec->n))
    {
     // if(j<(child_elec->n-1))
      if(j%factor<factor/2)
      {
		            basex = child_elec->locx+(i)/factor;
                basey = child_elec->locy+j/factor-1;
                distx = (double)(i%factor)/factor;
                //disty = 0.5 + 0.5/factor;
                disty = 0.5 + 0.5/factor + (double)(j%factor)/factor;
                c_eys[i][j] = (1-tau)*(interpolate2d(eys_old[basex][basey],eys_old[basex+1][basey],eys_old[basex+1][child_elec->locy],eys_old[basex][child_elec->locy],distx,disty))
                                                +(tau)*(interpolate2d(eys[basex][basey],eys[basex+1][basey],eys[basex+1][basey+1],eys[basex][basey+1],distx,disty));

                c_vy[i][j] = (1-tau)*(interpolate2d(vy_old[basex][basey],vy_old[basex+1][basey],vy_old[basex+1][basey+1],vy_old[basex][basey+1],distx,disty))
                                                +(tau)*(interpolate2d(vy[basex][basey],vy[basex+1][basey],vy[basex+1][basey+1],vy[basex][basey+1],distx,disty));
      }
      
      else
      {
          if(child_elec->locy+j/factor+1<(root_elec->n))
          {
		            basex = child_elec->locx+(i)/factor;
                basey = child_elec->locy+j/factor;
                distx = (double)(i%factor)/factor;
                //disty = 0.5 + 0.5/factor;
                disty = 0.5/factor+ (double)(j%factor)/factor - 0.5;
                c_eys[i][j] = (1-tau)*(interpolate2d(eys_old[basex][basey],eys_old[basex+1][basey],eys_old[basex+1][child_elec->locy],eys_old[basex][child_elec->locy],distx,disty))
                                                +(tau)*(interpolate2d(eys[basex][basey],eys[basex+1][basey],eys[basex+1][basey+1],eys[basex][basey+1],distx,disty));

                c_vy[i][j] = (1-tau)*(interpolate2d(vy_old[basex][basey],vy_old[basex+1][basey],vy_old[basex+1][basey+1],vy_old[basex][basey+1],distx,disty))
                                                +(tau)*(interpolate2d(vy[basex][basey],vy[basex+1][basey],vy[basex+1][basey+1],vy[basex][basey+1],distx,disty));
      
         }
      }
      
      
    }
    
	if(i%factor<factor/2)
        {
         
                basex = child_elec->locx+(i)/factor-1;
                basey = child_elec->locy;
                distx = 0.5 + 0.5/factor + (double)(i%factor)/factor;
                c_exs[i][0] = (1-tau)*(interpolate1d(exs_old[basex][basey],exs_old[basex+1][basey],distx))
                                                +(tau)*(interpolate1d(exs[basex][basey],exs[basex+1][basey],distx));

                c_vx[i][0] = (1-tau)*(interpolate1d(vx_old[basex][basey],vx_old[basex+1][basey],distx))
                                                +(tau)*(interpolate1d(vx[basex][basey],vx[basex+1][basey],distx));
                                                
           
        }
	
	else
        {
           
                basex = child_elec->locx+(i)/factor;
                basey = child_elec->locy;
                distx = 0.5/factor + (double)(i%factor)/factor - 0.5;
                c_exs[i][0] = (1-tau)*(interpolate1d(exs_old[basex][basey],exs_old[basex+1][basey],distx))
                                                +(tau)*(interpolate1d(exs[basex][basey],exs[basex+1][basey],distx));

                c_vx[i][0] = (1-tau)*(interpolate1d(vx_old[basex][basey],vx_old[basex+1][basey],distx))
                                                +(tau)*(interpolate1d(vx[basex][basey],vx[basex+1][basey],distx));
          
        }
	
	if (i%factor<factor/2)
        {
          
                basex = child_elec->locx+(i)/factor-1;
                basey = child_elec->locy+(child_elec->n)/factor;
                distx = 0.5 + 0.5/factor + (double)(i%factor)/factor;
                c_exs[i][child_elec->n] = (1-tau)*(interpolate1d(exs_old[basex][basey],exs_old[basex+1][basey],distx))
                                                        +(tau)*(interpolate1d(exs[basex][basey],exs[basex+1][basey],distx));

                c_vx[i][child_elec->n] = (1-tau)*(interpolate1d(vx_old[basex][basey],vx_old[basex+1][basey],distx))
                                                        +(tau)*(interpolate1d(vx[basex][basey],vx[basex+1][basey],distx));
          
        }
	
	else
        {
          
                basex = child_elec->locx+(i)/factor;
                basey = child_elec->locy+(child_elec->n)/factor;
                distx = 0.5/factor + (double)(i%factor)/factor - 0.5;

                c_exs[i][child_elec->n] = (1-tau)*(interpolate1d(exs_old[basex][basey],exs_old[basex+1][basey],distx))
                                                        +(tau)*(interpolate1d(exs[basex][basey],exs[basex+1][basey],distx));

                c_vx[i][child_elec->n] = (1-tau)*(interpolate1d(vx_old[basex][basey],vx_old[basex+1][basey],distx))
                                                        +(tau)*(interpolate1d(vx[basex][basey],vx[basex+1][basey],distx));
          
        }
	
       	basex = child_elec->locx+(i)/factor;
        basey = child_elec->locy;
        distx = (double)(i%factor)/factor;
        child_den->mesh[i][0] = (1-tau)*(interpolate1d(denp[basex][basey],denp[basex+1][basey],distx))
                                                        +(tau)*(interpolate1d(root_den->mesh[basex][basey],root_den->mesh[basex+1][basey],distx));

        basex = child_elec->locx+(i)/factor;
        basey = child_elec->locy+(child_elec->n)/factor;
        distx = (double)(i%factor)/factor;
        child_den->mesh[i][child_elec->n] =(1-tau)*(interpolate1d(denp[basex][basey],denp[basex+1][basey],distx))
                                                                                +(tau)*(interpolate1d(root_den->mesh[basex][basey],root_den->mesh[basex+1][basey],distx));
 
 
   
       
          if((child_elec->locy>0) && (j<child_elec->n))
          {
            if(j%factor<factor/2)
            {
              
              if (i%factor<factor/2)
              {  
        
        
                basex = child_elec->locx+(i)/factor-1;
                basey = child_elec->locy-1+j/factor;
                distx = 0.5 + 0.5/factor + (double)(i%factor)/factor;
               // disty = 0.5 + 0.5/factor;
               disty = 0.5 + 0.5/factor + (double)(j%factor)/factor;
                child_mag->mesh[i][j] = (1-tau)*(interpolate2d(hzi[basex][basey],hzi[basex+1][basey],hzi[basex+1][basey+1],hzi[basex][basey+1],distx,disty))
                                                                +(tau)*(interpolate2d(root_mag->mesh[basex][basey],root_mag->mesh[basex+1][basey],root_mag->mesh[basex+1][basey+1],root_mag->mesh[basex][basey+1],distx,disty));
                                                                
              }
            
            else
            {
                
              basex = child_elec->locx+i/factor;
              basey = child_elec->locy+j/factor-1;
              distx = 0.5/factor + (double)(i%factor)/factor - 0.5;
              //disty = 0.5 - 0.5/factor;
              disty = 0.5 + 0.5/factor + (double)(j%factor)/factor;
              child_mag->mesh[i][j] = (1-tau)*(interpolate2d(hzi[basex][basey],hzi[basex+1][basey],hzi[basex+1][basey+1],hzi[basex][basey+1],distx,disty))

               +(tau)*(interpolate2d(root_mag->mesh[basex][basey],root_mag->mesh[basex+1][basey],root_mag->mesh[basex+1][basey+1],root_mag->mesh[basex][basey+1],distx,disty));
                
             }  
            
            }    
                
       else    
        {
            if (i%factor<factor/2)
             {  
                basex = child_elec->locx+(i)/factor-1;
                basey = child_elec->locy+j/factor;
                distx = 0.5 + 0.5/factor + (double)(i%factor)/factor; 
                disty = 0.5/factor+ (double)(j%factor)/factor - 0.5; 
                child_mag->mesh[i][j] = (1-tau)*(interpolate2d(hzi[basex][basey],hzi[basex+1][basey],hzi[basex+1][basey+1],hzi[basex][basey+1],distx,disty))
                                                                +(tau)*(interpolate2d(root_mag->mesh[basex][basey],root_mag->mesh[basex+1][basey],root_mag->mesh[basex+1][basey+1],root_mag->mesh[basex][basey+1],distx,disty));
            
              }
            else
            {
              basex = child_elec->locx+i/factor;
              basey = child_elec->locy+j/factor;
              distx = 0.5/factor + (double)(i%factor)/factor - 0.5;
              //disty = 0.5 - 0.5/factor;
               disty = 0.5/factor+ (double)(j%factor)/factor - 0.5; 
               child_mag->mesh[i][j] = (1-tau)*(interpolate2d(hzi[basex][basey],hzi[basex+1][basey],hzi[basex+1][basey+1],hzi[basex][basey+1],distx,disty))

               +(tau)*(interpolate2d(root_mag->mesh[basex][basey],root_mag->mesh[basex+1][basey],root_mag->mesh[basex+1][basey+1],root_mag->mesh[basex][basey+1],distx,disty));
            }
        
          }        
      }              
                
                
          
         if(j<(child_elec->n))
         {
           	c_eyt[i][j] = c_eys[i][j]+c_eyi1[i][j];
            //c_eyt[i][child_elec->n-1] = c_eys[i][child_elec->n-1]+c_eyi1[i][child_elec->n-1];
         }
        
         if(j<=(child_elec->n))
         { 
            c_ext[i][j] = c_exs[i][j]+c_exi1[i][j];
          
         }
      
      
   }	
          

}



	// //========================================================================================================================================

}

void interpolatey(double t1)
{
     double tau = t1;
     int basex, basey;
     double distx,disty;
     
for(i=1;i<child_elec->m;i++)
{
	
    double ita,ep;
   if(child_elec->locy>0)
   {
		 basex = child_elec->locx+(i)/factor;
                basey = child_elec->locy-1;
                distx = (double)(i%factor)/factor;
                disty = 0.5 + 0.5/factor;
                c_eys[i][0] = (1-tau)*(interpolate2d(eys_old[basex][basey],eys_old[basex+1][basey],eys_old[basex+1][child_elec->locy],eys_old[basex][child_elec->locy],distx,disty))
                                                +(tau)*(interpolate2d(eys[basex][basey],eys[basex+1][basey],eys[basex+1][basey+1],eys[basex][basey+1],distx,disty));

                c_vy[i][0] = (1-tau)*(interpolate2d(vy_old[basex][basey],vy_old[basex+1][basey],vy_old[basex+1][basey+1],vy_old[basex][basey+1],distx,disty))
                                                +(tau)*(interpolate2d(vy[basex][basey],vy[basex+1][basey],vy[basex+1][basey+1],vy[basex][basey+1],distx,disty));
   }
   
        basex = child_elec->locx+(i)/factor;
        basey = child_elec->locy+(child_elec->n-1)/factor;
        distx = (double)(i%factor)/factor;
        disty = 0.5 - 0.5/factor;
        c_eys[i][child_elec->n-1] = (1-tau)*(interpolate2d(eys_old[basex][basey],eys_old[basex+1][basey],eys_old[basex+1][basey+1],eys_old[basex][basey+1],distx,disty))
                                                                +(tau)*(interpolate2d(eys[basex][basey],eys[basex+1][basey],eys[basex+1][basey+1],eys[basex][basey+1],distx,disty));

        c_vy[i][child_elec->n-1] = (1-tau)*(interpolate2d(vy_old[basex][basey],vy_old[basex+1][basey],vy_old[basex+1][basey+1],vy_old[basex][basey+1],distx,disty))
                                                                +(tau)*(interpolate2d(vy[basex][basey],vy[basex+1][basey],vy[basex+1][basey+1],vy[basex][basey+1],distx,disty));
   
	if(i%factor<factor/2)
        {
                basex = child_elec->locx+(i)/factor-1;
                basey = child_elec->locy;
                distx = 0.5 + 0.5/factor + (double)(i%factor)/factor;
                c_exs[i][0] = (1-tau)*(interpolate1d(exs_old[basex][basey],exs_old[basex+1][basey],distx))
                                                +(tau)*(interpolate1d(exs[basex][basey],exs[basex+1][basey],distx));

                c_vx[i][0] = (1-tau)*(interpolate1d(vx_old[basex][basey],vx_old[basex+1][basey],distx))
                                                +(tau)*(interpolate1d(vx[basex][basey],vx[basex+1][basey],distx));
        }
	
	else
        {
                basex = child_elec->locx+(i)/factor;
                basey = child_elec->locy;
                distx = 0.5/factor + (double)(i%factor)/factor - 0.5;
                c_exs[i][0] = (1-tau)*(interpolate1d(exs_old[basex][basey],exs_old[basex+1][basey],distx))
                                                +(tau)*(interpolate1d(exs[basex][basey],exs[basex+1][basey],distx));

                c_vx[i][0] = (1-tau)*(interpolate1d(vx_old[basex][basey],vx_old[basex+1][basey],distx))
                                                +(tau)*(interpolate1d(vx[basex][basey],vx[basex+1][basey],distx));
        }
	
	if (i%factor<factor/2)
        {
                basex = child_elec->locx+(i)/factor;
                basey = child_elec->locy+(child_elec->n)/factor;
                distx = 0.5 + 0.5/factor + (double)(i%factor)/factor;
                c_exs[i][child_elec->n] = (1-tau)*(interpolate1d(exs_old[basex][basey],exs_old[basex+1][basey],distx))
                                                        +(tau)*(interpolate1d(exs[basex][basey],exs[basex+1][basey],distx));

                c_vx[i][child_elec->n] = (1-tau)*(interpolate1d(vx_old[basex][basey],vx_old[basex+1][basey],distx))
                                                        +(tau)*(interpolate1d(vx[basex][basey],vx[basex+1][basey],distx));
        }
	
	else
        {
                basex = child_elec->locx+(i)/factor;
                basey = child_elec->locy+(child_elec->n)/factor;
                distx = 0.5/factor + (double)(i%factor)/factor - 0.5;

                c_exs[i][child_elec->n] = (1-tau)*(interpolate1d(exs_old[basex][basey],exs_old[basex+1][basey],distx))
                                                        +(tau)*(interpolate1d(exs[basex][basey],exs[basex+1][basey],distx));

                c_vx[i][child_elec->n] = (1-tau)*(interpolate1d(vx_old[basex][basey],vx_old[basex+1][basey],distx))
                                                        +(tau)*(interpolate1d(vx[basex][basey],vx[basex+1][basey],distx));
        }
	
	basex = child_elec->locx+(i)/factor;
        basey = child_elec->locy;
        distx = (double)(i%factor)/factor;
        child_den->mesh[i][0] = (1-tau)*(interpolate1d(denp[basex][basey],denp[basex+1][basey],distx))
                                                        +(tau)*(interpolate1d(root_den->mesh[basex][basey],root_den->mesh[basex+1][basey],distx));

        basex = child_elec->locx+(i)/factor;
        basey = child_elec->locy+(child_elec->n)/factor;
        distx = (double)(i%factor)/factor;
        child_den->mesh[i][child_elec->n] =(1-tau)*(interpolate1d(denp[basex][basey],denp[basex+1][basey],distx))
                                                                                +(tau)*(interpolate1d(root_den->mesh[basex][basey],root_den->mesh[basex+1][basey],distx));
	
	if (i%factor<factor/2 && child_mag->locy>0)
        {
                basex = child_elec->locx+(i)/factor-1;
                basey = child_elec->locy-1;
                distx = 0.5 + 0.5/factor + (double)(i%factor)/factor;
                disty = 0.5 + 0.5/factor;
                child_mag->mesh[i][0] = (1-tau)*(interpolate2d(hzi[basex][basey],hzi[basex+1][basey],hzi[basex+1][basey+1],hzi[basex][basey+1],distx,disty))
                                                                +(tau)*(interpolate2d(root_mag->mesh[basex][basey],root_mag->mesh[basex+1][basey],root_mag->mesh[basex+1][basey+1],root_mag->mesh[basex][basey+1],distx,disty));
        }
	
	else if(child_mag->locy>0)
        {
                basex = child_elec->locx+(i)/factor;
                basey = child_elec->locy-1;
                distx = 0.5/factor + (double)(i%factor)/factor - 0.5;
                disty = 0.5 + 0.5/factor;
                child_mag->mesh[i][0] = (1-tau)*(interpolate2d(hzi[basex][basey],hzi[basex+1][basey],hzi[basex+1][basey+1],hzi[basex][basey+1],distx,disty))
                                                                +(tau)*(interpolate2d(root_mag->mesh[basex][basey],root_mag->mesh[basex+1][basey],root_mag->mesh[basex+1][basey+1],root_mag->mesh[basex][basey+1],distx,disty));
        }
	
	if (i%factor<factor/2)
        {
                basex = child_elec->locx+i/factor-1;
        basey = child_elec->locy+(child_elec->n-1)/factor;
        distx = 0.5 + 0.5/factor + (double)(i%factor)/factor;
        disty = 0.5 - 0.5/factor;
                child_mag->mesh[i][child_elec->n-1] = (1-tau)*(interpolate2d(hzi[basex][basey],hzi[basex+1][basey],hzi[basex+1][basey+1],hzi[basex][basey+1],distx,disty))
                                                                                                +(tau)*(interpolate2d(root_mag->mesh[basex][basey],root_mag->mesh[basex+1][basey],root_mag->mesh[basex+1][basey+1],root_mag->mesh[basex][basey+1],distx,disty));
        }
	
	else
        {
                basex = child_elec->locx+i/factor;
        basey = child_elec->locy+(child_elec->n-1)/factor;
        distx = 0.5/factor + (double)(i%factor)/factor - 0.5;
        disty = 0.5 - 0.5/factor;
                child_mag->mesh[i][child_elec->n-1] = (1-tau)*(interpolate2d(hzi[basex][basey],hzi[basex+1][basey],hzi[basex+1][basey+1],hzi[basex][basey+1],distx,disty))

               +(tau)*(interpolate2d(root_mag->mesh[basex][basey],root_mag->mesh[basex+1][basey],root_mag->mesh[basex+1][basey+1],root_mag->mesh[basex][basey+1],distx,disty));
        }
	
	c_eyt[i][0] = c_eys[i][0]+c_eyi1[i][0];
    c_eyt[i][child_elec->n-1] = c_eys[i][child_elec->n-1]+c_eyi1[i][child_elec->n-1];
    c_ext[i][0] = c_exs[i][0]+c_exi1[i][0];
    c_ext[i][child_elec->n] = c_exs[i][child_elec->n]+c_exi1[i][child_elec->n];
	
}



	// //========================================================================================================================================

}

void c2p()
{	
    FILE *fpr;
    int l=1;
    int ic,jc;
    int cm1=0,cn1=0,cm2=0,cn2=0,cm3=0,cn3=0,cm4=0,cn4=0;
    int centr;
    gettimeofday(&begin,NULL);
    
    centr=(factor)/2;
    
     for(i=factor;i<child_elec->m-factor;i+=factor)
        {
        
              
                for(j=factor;j<child_elec->n;j+=factor)
                {
                   
                   
                   exs[child_elec->locx+i/factor][child_elec->locy+j/factor]
                 = (c_exs[i+centr-1][j]+c_exs[i+centr][j])/2;
                 
                 vx[child_elec->locx+i/factor][child_elec->locy+j/factor]
                 = (c_vx[i+centr-1][j]+c_vx[i+centr][j])/2;
                
                ext[child_elec->locx+i/factor][child_elec->locy+j/factor]
                 = exs[child_elec->locx+i/factor][child_elec->locy+j/factor]
                + exi1[child_elec->locx+i/factor][child_elec->locy+j/factor];
                 

                }
    
              
        }
        
	
	for(i=factor;i<child_elec->m;i+=factor)
        {
            
                for(j=factor;j<child_elec->n-(factor);j+=(factor))
                {    
                    
                    
                   eys[child_elec->locx+i/factor][child_elec->locy+j/factor]
                   = (c_eys[i][j+centr-1]+c_eys[i][j+centr])/2;
                   
                  eyt[child_elec->locx+i/factor][child_elec->locy+j/factor]
                   = eys[child_elec->locx+i/factor][child_elec->locy+j/factor]
                   + eyi1[child_elec->locx+i/factor][child_elec->locy+j/factor];
               
                   vy[child_elec->locx+i/factor][child_elec->locy+j/factor]
                   = (c_vy[i][j+centr-1]+c_vy[i][j+centr])/2;
                 
                    
                 }
                 
        }
	
	for(i=factor;i<child_elec->m;i+=factor)
        {
                for(j=factor;j<child_elec->n;j+=factor)
                {
                
                                
                 
                 root_den->mesh[child_den->locx+(i)/factor][child_den->locy+(j)/factor]
              = child_den->mesh[i][j];
                
                   
            
            
                }
        }
	
	for(i=factor;i<child_elec->m-(factor);i+=(factor))
        {
              
        
                for(j=factor;j<child_elec->n-(factor);j+=(factor))
                {
                    
              
		
		root_mag->mesh[child_mag->locx+i/factor][child_mag->locy+j/factor]
              = (child_mag->mesh[i+centr-1][j+centr-1]
                +child_mag->mesh[i+centr][j+centr-1]
                +child_mag->mesh[i+centr-1][j+centr]
                +child_mag->mesh[i+centr][j+centr])/4;
                    
                   
                }
                
        } 
	 
    gettimeofday(&end,NULL);
    t_cal_c2p += ((end.tv_sec - begin.tv_sec) + ((end.tv_usec - begin.tv_usec)/1000000.0));

}
