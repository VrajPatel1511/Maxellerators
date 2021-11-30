#include "xyfdtd.h"
#include <cuda.h>

int tag=1;
double c=2.99795e8,pi=3.14159265,eps0=8.854e-12,xmu0=1.55063706e-6,qe=1.602176487e-19,cmasse=9.10938215e-31,akb=1.3806503e-23;
double radius=0;
int rank =0, sizar=3000;
struct node *newnode(int m, int n ,int locx, int locy, int level)
{
	struct node *temp = malloc(sizeof(struct node));
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

	ERMSp = malloc(sizeof(double *) * sizar);
	if (ERMSp){
		for (i = 0; i < sizar; i++){
			ERMSp[i] = malloc(sizeof(double) * sizar);
		}
	}
	erms2 = malloc(sizeof(double *) * sizar);
	if (erms2){
		for (i = 0; i < sizar; i++){
			erms2[i] = malloc(sizeof(double) * sizar);
		}
	}
	// den = malloc(sizeof(double *) * sizar);
	// if (den){
	// 	for (i = 0; i < sizar; i++){
	// 		den[i] = malloc(sizeof(double) * sizar);
	// 	}
	// }
	
	exi = malloc(sizeof(double *) * sizar);
	if (exi){
		for (i = 0; i < sizar; i++){
			exi[i] = malloc(sizeof(double) * sizar);
		}
	}
	eyi = malloc(sizeof(double *) * sizar);
	if (eyi){
		for (i = 0; i < sizar; i++){
			eyi[i] = malloc(sizeof(double) * sizar);
		}
	}
	exi1 = malloc(sizeof(double *) * sizar);
	if (exi1){
		for (i = 0; i < sizar; i++){
			exi1[i] = malloc(sizeof(double) * sizar);
		}
	}
	eyi1 = malloc(sizeof(double *) * sizar);
	if (eyi1){
		for (i = 0; i < sizar; i++){
			eyi1[i] = malloc(sizeof(double) * sizar);
		}
	}
	exs = malloc(sizeof(double *) * sizar);
	if (exs){
		for (i = 0; i < sizar; i++){
			exs[i] = malloc(sizeof(double) * sizar);
		}
	}
	eys = malloc(sizeof(double *) * sizar);
	if (eys){
		for (i = 0; i < sizar; i++){
			eys[i] = malloc(sizeof(double) * sizar);
		}
	}
	hzi = malloc(sizeof(double *) * sizar);
	if (hzi){
		for (i = 0; i < sizar; i++){
			hzi[i] = malloc(sizeof(double) * sizar);
		}
	}
	vx = malloc(sizeof(double *) * sizar);
	if (vx){
		for (i = 0; i < sizar; i++){
			vx[i] = malloc(sizeof(double) * sizar);
		}
	}
	vy = malloc(sizeof(double *) * sizar);
	if (vy){
		for (i = 0; i < sizar; i++){
			vy[i] = malloc(sizeof(double) * sizar);
		}
	}
	ext = malloc(sizeof(double *) * sizar);
	if (ext){
		for (i = 0; i < sizar; i++){
			ext[i] = malloc(sizeof(double) * sizar);
		}
	}
	eyt = malloc(sizeof(double *) * sizar);
	if (eyt){
		for (i = 0; i < sizar; i++){
			eyt[i] = malloc(sizeof(double) * sizar);
		}
	}
	exs1 = malloc(sizeof(double *) * sizar);
	if (exs1){
		for (i = 0; i < sizar; i++){
			exs1[i] = malloc(sizeof(double) * sizar);
		}
	}
	eys1 = malloc(sizeof(double *) * sizar);
	if (eys1){
		for (i = 0; i < sizar; i++){
			eys1[i] = malloc(sizeof(double) * sizar);
		}
	}
	xmid = malloc(sizeof(double) * sizar);
	ymid = malloc(sizeof(double) * sizar);
	sgdx0 = malloc(sizeof(double) * sizar);
	sgdy0 = malloc(sizeof(double) * sizar);
	DINI = malloc(sizeof(double) * sizar);
	DIFFUSION = malloc(sizeof(double *) * sizar);
	if (DIFFUSION){
		for (i = 0; i < sizar; i++){
			DIFFUSION[i] = malloc(sizeof(double) * sizar);
		}
	}
	frqio = malloc(sizeof(double *) * sizar);
	if (frqio){
		for (i = 0; i < sizar; i++){
			frqio[i] = malloc(sizeof(double) * sizar);
		}
	}
	denp = malloc(sizeof(double *) * sizar);
	if (denp){
		for (i = 0; i < sizar; i++){
			denp[i] = malloc(sizeof(double) * sizar);
		}
	}
	
	imid = malloc(sizeof(int) * sizar);
	jmid = malloc(sizeof(int) * sizar);

	exs_old = malloc(sizeof(double *) * sizar);
	if (exs_old){
		for (i = 0; i < sizar; i++){
			exs_old[i] = malloc(sizeof(double) * sizar);
		}
	}
	eys_old = malloc(sizeof(double *) * sizar);
	if (eys_old){
		for (i = 0; i < sizar; i++){
			eys_old[i] = malloc(sizeof(double) * sizar);
		}
	}
	vx_old = malloc(sizeof(double *) * sizar);
	if (vx_old){
		for (i = 0; i < sizar; i++){
			vx_old[i] = malloc(sizeof(double) * sizar);
		}
	}
	vy_old = malloc(sizeof(double *) * sizar);
	if (vy_old){
		for (i = 0; i < sizar; i++){
			vy_old[i] = malloc(sizeof(double) * sizar);
		}
	}

	root_mesh_elec = malloc(sizeof(double)*sizar);
	child_mesh_elec = malloc(sizeof(double)*sizar);
	if(root_mesh_elec && child_mesh_elec)
	{
		for (i = 0; i < sizar; ++i)
		{
			root_mesh_elec[i] = calloc(sizar, sizeof(double));
			child_mesh_elec[i] = calloc(sizar, sizeof(double));	
		}
	}

	root_mesh_mag = malloc(sizeof(double)*sizar);
	child_mesh_mag = malloc(sizeof(double)*sizar);
	if(root_mesh_mag && child_mesh_mag)
	{
		for (i = 0; i < sizar; ++i)
		{
			root_mesh_mag[i] = calloc(sizar, sizeof(double));
			child_mesh_mag[i] = calloc(sizar, sizeof(double));	
		}
	}

	root_mesh_den = malloc(sizeof(double)*sizar);
	child_mesh_den = malloc(sizeof(double)*sizar);
	if(root_mesh_den && child_mesh_den)
	{
		for (i = 0; i < sizar; ++i)
		{
			root_mesh_den[i] = calloc(sizar, sizeof(double));
			child_mesh_den[i] = calloc(sizar, sizeof(double));	
		}
	}

	c_exs = malloc(sizeof(double *) * sizar);
	if (c_exs){
		for (i = 0; i < sizar; i++){
			c_exs[i] = calloc(sizar, sizeof(double));
		}
	}
	c_eys = malloc(sizeof(double *) * sizar);
	if (c_eys){
		for (i = 0; i < sizar; i++){
			c_eys[i] = calloc(sizar, sizeof(double));
		}
	}
	c_eyi = malloc(sizeof(double *) * sizar);
	if (c_eyi){
		for (i = 0; i < sizar; i++){
			c_eyi[i] = calloc(sizar, sizeof(double));
		}
	}
	c_eyi1 = malloc(sizeof(double *) * sizar);
	if (c_eyi1){
		for (i = 0; i < sizar; i++){
			c_eyi1[i] = calloc(sizar, sizeof(double));
		}
	}

	c_erms2 = malloc(sizeof(double *) * sizar);
	if (c_erms2){
		for (i = 0; i < sizar; i++){
			c_erms2[i] = malloc(sizeof(double) * sizar);
		}
	}
	// den = malloc(sizeof(double *) * sizar);
	// if (den){
	// 	for (i = 0; i < sizar; i++){
	// 		den[i] = malloc(sizeof(double) * sizar);
	// 	}
	// }
	c_vx = malloc(sizeof(double *) * sizar);
	if (c_vx){
		for (i = 0; i < sizar; i++){
			c_vx[i] = malloc(sizeof(double) * sizar);
		}
	}
	c_vy = malloc(sizeof(double *) * sizar);
	if (c_vy){
		for (i = 0; i < sizar; i++){
			c_vy[i] = malloc(sizeof(double) * sizar);
		}
	}
	c_exi = malloc(sizeof(double *) * sizar);
	if (c_exi){
		for (i = 0; i < sizar; i++){
			c_exi[i] = malloc(sizeof(double) * sizar);
		}
	}
	c_exi1 = malloc(sizeof(double *) * sizar);
	if (c_exi1){
		for (i = 0; i < sizar; i++){
			c_exi1[i] = malloc(sizeof(double) * sizar);
		}
	}
	c_ext = malloc(sizeof(double *) * sizar);
	if (c_ext){
		for (i = 0; i < sizar; i++){
			c_ext[i] = malloc(sizeof(double) * sizar);
		}
	}
	c_eyt = malloc(sizeof(double *) * sizar);
	if (c_eyt){
		for (i = 0; i < sizar; i++){
			c_eyt[i] = malloc(sizeof(double) * sizar);
		}
	}
	c_exs1 = malloc(sizeof(double *) * sizar);
	if (c_exs1){
		for (i = 0; i < sizar; i++){
			c_exs1[i] = malloc(sizeof(double) * sizar);
		}
	}
	c_eys1 = malloc(sizeof(double *) * sizar);
	if (c_eys1){
		for (i = 0; i < sizar; i++){
			c_eys1[i] = malloc(sizeof(double) * sizar);
		}
	}
	c_xmid = malloc(sizeof(double) * sizar);
	c_ymid = malloc(sizeof(double) * sizar);
	c_sgdx0 = malloc(sizeof(double) * sizar);
	c_sgdy0 = malloc(sizeof(double) * sizar);
	c_DINI = malloc(sizeof(double) * sizar);
	c_DIFFUSION = malloc(sizeof(double *) * sizar);
	if (c_DIFFUSION){
		for (i = 0; i < sizar; i++){
			c_DIFFUSION[i] = malloc(sizeof(double) * sizar);
		}
	}
	c_frqio = malloc(sizeof(double *) * sizar);
	if (c_frqio){
		for (i = 0; i < sizar; i++){
			c_frqio[i] = malloc(sizeof(double) * sizar);
		}
	}
	c_denp = malloc(sizeof(double *) * sizar);
	if (c_denp){
		for (i = 0; i < sizar; i++){
			c_denp[i] = malloc(sizeof(double) * sizar);
		}
	}

  c_denpold = malloc(sizeof(double *) * sizar);
	if (c_denpold){
		for (i = 0; i < sizar; i++){
			c_denpold[i] = malloc(sizeof(double) * sizar);
		}
	}
  
  c_exold = malloc(sizeof(double *) * sizar);
	if (c_exold){
		for (i = 0; i < sizar; i++){
			c_exold[i] = malloc(sizeof(double) * sizar);
		}
	}
  
   c_eyold = malloc(sizeof(double *) * sizar);
	if (c_eyold){
		for (i = 0; i < sizar; i++){
			c_eyold[i] = malloc(sizeof(double) * sizar);
		}
	}
 
  c_vxold = malloc(sizeof(double *) * sizar);
	if (c_vxold){
		for (i = 0; i < sizar; i++){
			c_vxold[i] = malloc(sizeof(double) * sizar);
		}
	} 
 
  c_vyold = malloc(sizeof(double *) * sizar);
	if (c_vyold){
		for (i = 0; i < sizar; i++){
			c_vyold[i] = malloc(sizeof(double) * sizar);
		}
	}
  
  c_hzold = malloc(sizeof(double *) * sizar);
	if (c_hzold){
		for (i = 0; i < sizar; i++){
			c_hzold[i] = malloc(sizeof(double) * sizar);
		}
	}
   
	c_hzi = malloc(sizeof(double *) * sizar);
	if (c_hzi){
		for (i = 0; i < sizar; i++){
			c_hzi[i] = malloc(sizeof(double) * sizar);
		}
	}
	
	c_imid = malloc(sizeof(int) * sizar);
	c_jmid = malloc(sizeof(int) * sizar);
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
