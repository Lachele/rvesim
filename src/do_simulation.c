#include "rvesim.h"

/**************** do_simulation *****************/ 

/* This function loops through all the line emissions (whether calculated in a 
   transition function or read from a file) and adds their contributions to 
   appropriate places in the simulation.  The distinction between an 
   instrument that scans or jumps from point to point is made here.  */
void do_simulation(){
/* counting variables */
int da=0,db=0,dc=0,dd=0,dn=0,dbin=0,dnn=0,dsw=0;
/* flags */
int ddoflag=0,ddoneflag=0;
/* variables for assigning intensities */
double dposn=0,dcslope=0,dnslope=0,ddumposn=0,ddumnint=0,ddumcint=0;
double ddend=0,dwidth=0,dintmax=0;
/* temporary simulations for writing separate atomic, molecular,
   etc, simulations to files */
sim dmnsim,dmcsim,dasim,dosim,dtotsim;
/* variables about files to create */
char dnfile[2000],dcfile[2000];
FILE *DNFILE,*DCFILE;

/* Make the directory for storing simulation output */
sprintf(dnfile,"mkdir %s_simulations",PREF);
system(dnfile); 
/* get number of bins and allocate memory for simulations */
dn=(int)(ceil((MAXV-MINV)/PSTEP));
if(ddebugflag<DEBUG){
fprintf(DBG,"Top of do_simulation.  dn=%d\n",dn);
fflush(DBG);
}
dtotsim.n=dmnsim.n=dmcsim.n=dasim.n=dosim.n=dn;
dmnsim.x=(double*)calloc(dn,sizeof(double));
dmnsim.y=(double*)calloc(dn,sizeof(double)); 
dmcsim.x=(double*)calloc(dn,sizeof(double));
dmcsim.y=(double*)calloc(dn,sizeof(double)); 
dasim.x=(double*)calloc(dn,sizeof(double));
dasim.y=(double*)calloc(dn,sizeof(double)); 
dosim.x=(double*)calloc(dn,sizeof(double));
dosim.y=(double*)calloc(dn,sizeof(double)); 
dtotsim.x=(double*)calloc(dn,sizeof(double)); 
dtotsim.y=(double*)calloc(dn,sizeof(double)); 

/******************* Diatomic emissions ******************/
for(da=0;da<NUMMOL;da++){/* Loop over all molecules */ 
	for(db=0;db<MOL[da].trans;db++){ /* then loop over transitions */
if(ddebugflag<DEBUG){
fprintf(DBG,"MOL[da].t[db].Nhi=%s, MOL[%d].t[%d].Nlo=%s",\
		MOL[da].t[db].Nhi,da,db,MOL[da].t[db].Nlo);
fflush(DBG);
}
		for(dc=0;dc<MOL[da].t[db].nS;dc++){/* loop through simsets */
if(ddebugflag<DEBUG){
fprintf(DBG,"da=%d, db=%d, dc=%d, MOL[da].t[db].fS[dc]=%d, ",\
		da,db,dc,MOL[da].t[db].fS[dc]);
fflush(DBG);
}
/* if inclusion flag says to include this simset */
		if(MOL[da].t[db].fS[dc]==0){
/* loop through transitions in simset */
if(ddebugflag<DEBUG){
fprintf(DBG,"MOL[%d].t[%d].S[%d].n=%d, ",\
		da,db,dc,MOL[da].t[db].S[dc].n); 
fflush(DBG);
}
			for(dd=0;dd<MOL[da].t[db].S[dc].n;dd++){ 
/* if the current frequency is not equal to zero */
			if(MOL[da].t[db].S[dc].f[dd]!=0){
/* see if line falls within requested simulation min/max */
if(UNITTYPE==0) dposn=1e7/MOL[da].t[db].S[dc].f[dd];
else dposn=MOL[da].t[db].S[dc].f[dd];
ddoflag=0;
if(dposn<(MINV-RESN)) ddoflag=1;
if(dposn>(MAXV+RESN)) ddoflag=1;
if(ddebugflag<DEBUG){
fprintf(DBG,"UNITTYPE=%d, dposn=%12.6e, ddoflag=%d\n",UNITTYPE,dposn,ddoflag);
fprintf(DBG,"MOL[%d].t[%d].S[%d].f[%d]=%12.6e\n",da,db,dc,dd,\
		MOL[da].t[db].S[dc].f[dd]);
fflush(DBG);
}
if(ddoflag==0){
/* find bin where line originates, also get native and cascade slopes */
	dbin=(int)(floor((dposn-MINV)/PSTEP));
	dcslope=MOL[da].pop*MOL[da].t[db].S[dc].ci[dd]/RESN;
	dnslope=MOL[da].pop*MOL[da].t[db].S[dc].ni[dd]/RESN; 
if(ddebugflag<DEBUG){
fprintf(DBG,"RESN=%12.6e, MOL[%d].pop=%12.6e\n",RESN,da,MOL[da].pop);
fprintf(DBG,"MOL[%d].t[%d].S[%d].ci[%d]=%12.6e\n",da,db,dc,dd,\
		MOL[da].t[db].S[dc].ci[dd]);
fprintf(DBG,"MOL[%d].t[%d].S[%d].ni[%d]=%12.6e\n",da,db,dc,dd,\
		MOL[da].t[db].S[dc].ni[dd]);
fprintf(DBG,"dposn=%12.6e, MINV=%12.6e, PSTEP=%12.6e\n",dposn,MINV,PSTEP);
fflush(DBG);
}
/* increment through them and calculate the intensities for each bin -- 
   integrate if SCAN and point if JUMP */
	ddoneflag=1;
	ddumposn=dposn;
	ddumnint=MOL[da].pop*MOL[da].t[db].S[dc].ni[dd]; 
	ddumcint=MOL[da].pop*MOL[da].t[db].S[dc].ci[dd]; 
	dnn=0;
	if(dbin<0) ddoneflag=0; 
if(ddebugflag<DEBUG){
fprintf(DBG,"dbin=%d, dcslope=%12.6e, dnslope=%12.6e, ddumnint=%12.6e, ",\
		dbin,dcslope,dnslope,ddumnint);
fprintf(DBG,"ddumcint=%12.6e,ddoneflag=%d\n",ddumcint,ddoneflag);
fflush(DBG);
}
	while(ddoneflag==1){ /* fill bins at/below dbin */ 
		dsw=0;
		ddend=MINV+(dbin-dnn)*PSTEP;
		if(ddumposn==dposn) dsw+=1;
		if(ddend<(dposn-RESN)) dsw+=2; 
		if(SCANTYPE==1) dsw+=4;
		if((dbin-dnn)>(dn-1)) dsw=8;
if(ddebugflag<DEBUG){
fprintf(DBG,"AT/BELOW  ddumposn=%12.6e, ddend=%12.6e, ",ddumposn,ddend);
fprintf(DBG,"SCANTYPE=%d, dsw=%d, ",SCANTYPE,dsw);
fflush(DBG);
}
		switch (dsw){
/* the following for equipment that scans between "points" */
			case 0:
				dmnsim.y[dbin-dnn]+=\
					(ddumnint-0.5*dnslope*PSTEP)*PSTEP;
				dmcsim.y[dbin-dnn]+=\
					(ddumcint-0.5*dcslope*PSTEP)*PSTEP;
				ddumnint-=dnslope*PSTEP; 
				ddumcint-=dcslope*PSTEP; 
				ddumposn-=PSTEP;
				if(dnn==dbin) ddoneflag=0;
				dnn++;
if(ddebugflag<DEBUG){
fprintf(DBG,"Case 0. dmnsim.y[%d]=%12.6e, dmcsim.y[%d]=%12.6e\n",\
	(dbin-dnn+1),dmnsim.y[dbin-dnn+1],(dbin-dnn+1),dmcsim.y[dbin-dnn+1]);
fprintf(DBG,"ddumnint=%12.6e, ddumcint=%12.6e, ddumposn=%12.6e, dnn=%d\n",\
		ddumnint,ddumcint,ddumposn,dnn);
fflush(DBG);
}
				break;
			case 1:
				dwidth=ddumposn-ddend; 
				dmnsim.y[dbin-dnn]+=\
					(ddumnint-0.5*dnslope*dwidth)*dwidth;
				dmcsim.y[dbin-dnn]+=\
					(ddumcint-0.5*dcslope*dwidth)*dwidth;
				ddumnint-=dnslope*dwidth; 
				ddumcint-=dcslope*dwidth; 
				ddumposn=ddend;
				if(dnn==dbin) ddoneflag=0;
				dnn++;
if(ddebugflag<DEBUG){
fprintf(DBG,"Case 1. dwidth=%12.6e, dmnsim.y[%d]=%12.6e, ",dwidth,\
	(dbin-dnn+1),dmnsim.y[dbin-dnn+1]);
fprintf(DBG,"dmcsim.y[%d]=%12.6e\n",(dbin-dnn+1),dmcsim.y[dbin-dnn+1]);
fprintf(DBG,"ddumnint=%12.6e, ddumcint=%12.6e, ddumposn=%12.6e, dnn=%d\n",\
		ddumnint,ddumcint,ddumposn,dnn);
fflush(DBG);
}
				break;
			case 2:
				if(dnslope!=0) dwidth=ddumnint/dnslope; 
				else{
					if(dcslope!=0) dwidth=ddumcint/dcslope; 
					else dwidth=0;
					}
				dmnsim.y[dbin-dnn]+=ddumnint*0.5*dwidth;
				dmcsim.y[dbin-dnn]+=ddumcint*0.5*dwidth;
				ddoneflag=0;
if(ddebugflag<DEBUG){
fprintf(DBG,"Case 2. dwidth=%12.6e, dmnsim.y[%d]=%12.6e, ",dwidth,\
	(dbin-dnn),dmnsim.y[dbin-dnn]);
fprintf(DBG,"dmcsim.y[%d]=%12.6e\n",(dbin-dnn),dmcsim.y[dbin-dnn]);
fprintf(DBG,"ddumnint=%12.6e, ddumcint=%12.6e, ddumposn=%12.6e, dnn=%d\n",\
		ddumnint,ddumcint,ddumposn,dnn);
fflush(DBG);
}
				break;
			case 3:
				dmnsim.y[dbin-dnn]+=ddumnint*0.5*PSTEP;
				dmcsim.y[dbin-dnn]+=ddumcint*0.5*PSTEP;
				ddoneflag=0;
if(ddebugflag<DEBUG){
fprintf(DBG,"Case 3. dmnsim.y[%d]=%12.6e, ",(dbin-dnn),dmnsim.y[dbin-dnn]);
fprintf(DBG,"dmcsim.y[%d]=%12.6e\n",(dbin-dnn),dmcsim.y[dbin-dnn]);
fprintf(DBG,"ddumnint=%12.6e, ddumcint=%12.6e, ddumposn=%12.6e, dnn=%d\n",\
		ddumnint,ddumcint,ddumposn,dnn);
fflush(DBG);
}
				break; 
/* the following for equipment that jumps from point to point */
			case 4:
				ddumnint-=dnslope*PSTEP; 
				ddumcint-=dcslope*PSTEP; 
				dmnsim.y[dbin-dnn]+=ddumnint;
				dmcsim.y[dbin-dnn]+=ddumcint;
				ddumposn-=PSTEP;
				if(dnn==dbin) ddoneflag=0;
				dnn++;
if(ddebugflag<DEBUG){
fprintf(DBG,"Case 4. dmnsim.y[%d]=%12.6e, dmcsim.y[%d]=%12.6e\n",\
	(dbin-dnn+1),dmnsim.y[dbin-dnn+1],(dbin-dnn+1),dmcsim.y[dbin-dnn+1]);
fprintf(DBG,"ddumnint=%12.6e, ddumcint=%12.6e, ddumposn=%12.6e, dnn=%d\n",\
		ddumnint,ddumcint,ddumposn,dnn);
fflush(DBG);
}
				break; 
			case 5:
				dwidth=ddumposn-ddend; 
				ddumnint-=dnslope*dwidth; 
				ddumcint-=dcslope*dwidth; 
				dmnsim.y[dbin-dnn]+=ddumnint;
				dmcsim.y[dbin-dnn]+=ddumcint;
				ddumposn=ddend;
				if(dnn==dbin) ddoneflag=0;
				dnn++;
if(ddebugflag<DEBUG){
fprintf(DBG,"Case 5. dwidth=%12.6e, dmnsim.y[%d]=%12.6e, ",dwidth,\
	(dbin-dnn+1),dmnsim.y[dbin-dnn+1]);
fprintf(DBG,"dmcsim.y[%d]=%12.6e\n",(dbin-dnn+1),dmcsim.y[dbin-dnn+1]);
fprintf(DBG,"ddumnint=%12.6e, ddumcint=%12.6e, ddumposn=%12.6e, dnn=%d\n",\
		ddumnint,ddumcint,ddumposn,dnn);
fflush(DBG);
}
				break; 
			case 6:
				ddoneflag=0;
if(ddebugflag<DEBUG){
fprintf(DBG,"Case 6\n");
fflush(DBG);
}
				break; 
			case 7:
				dmnsim.y[dbin-dnn]+=ddumnint;
				dmcsim.y[dbin-dnn]+=ddumcint;
				ddoneflag=0;
if(ddebugflag<DEBUG){
fprintf(DBG,"Case 7. dmnsim.y[%d]=%12.6e, ",(dbin-dnn),dmnsim.y[dbin-dnn]);
fprintf(DBG,"dmcsim.y[%d]=%12.6e\n",(dbin-dnn),dmcsim.y[dbin-dnn]);
fprintf(DBG,"ddumnint=%12.6e, ddumcint=%12.6e, ddumposn=%12.6e, dnn=%d\n",\
		ddumnint,ddumcint,ddumposn,dnn);
fflush(DBG);
}
				break; 
/* the following for emissions that provide only "spillover" */
			case 8:
				if(ddumposn==dposn){
					dwidth=ddumposn-ddend; 
					ddumnint-=dnslope*dwidth; 
					ddumcint-=dcslope*dwidth; 
					ddumposn=ddend;
					if(dnn==dbin) ddoneflag=0;
					dnn++;
					}
				else{
					if(ddend<(dposn-RESN)){
						printf(\
				"Fix switch case in do_sim at 1a.\n");
						exit(1);
						}
					else{
						ddumnint-=dnslope*PSTEP; 
						ddumcint-=dcslope*PSTEP; 
						ddumposn-=PSTEP;
						if(dnn==dbin) ddoneflag=0;
						dnn++;
						}
					}
if(ddebugflag<DEBUG){
fprintf(DBG,"Case 8. dwidth=%12.6e, dbin=%d",dwidth,dbin);
fprintf(DBG,"ddumnint=%12.6e, ddumcint=%12.6e, ddumposn=%12.6e, dnn=%d\n",\
		ddumnint,ddumcint,ddumposn,dnn);
fflush(DBG);
}
				break;
			default:
				printf("Go fix do_sim.\n");
				exit(1);
			}/* close switch-case */
		}/* close while for bins below the transition position */
	ddoneflag=1;
	ddumposn=dposn;
	ddumnint=MOL[da].pop*MOL[da].t[db].S[dc].ni[dd]; 
	ddumcint=MOL[da].pop*MOL[da].t[db].S[dc].ci[dd]; 
	dnn=0;
	if(dbin>dn-1) ddoneflag=0;
if(ddebugflag<DEBUG){
fprintf(DBG,"dbin=%d, dcslope=%12.6e, dnslope=%12.6e, ddumnint=%12.6e, ",\
		dbin,dcslope,dnslope,ddumnint);
fprintf(DBG,"ddumcint=%12.6e,ddoneflag=%d\n",ddumcint,ddoneflag);
fflush(DBG);
}
	while(ddoneflag==1){ /* fill bins at/above dbin */
		dsw=0;
		ddend=MINV+(dbin+dnn+1)*PSTEP;
		if(dposn==ddumposn) dsw+=1;
		if(ddend>(dposn+RESN)) dsw+=2;
		if(SCANTYPE==1) dsw+=4;
		if((dbin+dnn)<0) dsw=8;
if(ddebugflag<DEBUG){
fprintf(DBG,"AT/ABOVE ddumposn=%12.6e, ddend=%12.6e, ",ddumposn,ddend);
fprintf(DBG,"SCANTYPE=%d, dsw=%d, ",SCANTYPE,dsw);
fflush(DBG);
}
		switch (dsw){
/* the following for equipment that scans between "points" */
			case 0:
				dmnsim.y[dbin+dnn]+=\
					(ddumnint-0.5*dnslope*PSTEP)*PSTEP;
				dmcsim.y[dbin+dnn]+=\
					(ddumcint-0.5*dcslope*PSTEP)*PSTEP;
				ddumnint-=dnslope*PSTEP; 
				ddumcint-=dcslope*PSTEP; 
				ddumposn+=PSTEP;
				dnn++;
				if((dnn+dbin)==dn) ddoneflag=0;
if(ddebugflag<DEBUG){
fprintf(DBG,"Case 0. dmnsim.y[%d]=%12.6e, dmcsim.y[%d]=%12.6e\n",\
	(dbin+dnn-1),dmnsim.y[dbin+dnn-1],(dbin+dnn-1),dmcsim.y[dbin+dnn-1]);
fprintf(DBG,"ddumnint=%12.6e, ddumcint=%12.6e, ddumposn=%12.6e, dnn=%d\n",\
		ddumnint,ddumcint,ddumposn,dnn);
fflush(DBG);
}
				break;
			case 1:
				dwidth=ddend-ddumposn; 
				dmnsim.y[dbin+dnn]+=\
					(ddumnint-0.5*dnslope*dwidth)*dwidth;
				dmcsim.y[dbin+dnn]+=\
					(ddumcint-0.5*dcslope*dwidth)*dwidth;
				ddumnint-=dnslope*dwidth; 
				ddumcint-=dcslope*dwidth; 
				ddumposn=ddend;
				dnn++;
				if((dnn+dbin)==dn) ddoneflag=0;
if(ddebugflag<DEBUG){
fprintf(DBG,"Case 1. dwidth=%12.6e, dmnsim.y[%d]=%12.6e, ",dwidth,\
	(dbin+dnn-1),dmnsim.y[dbin+dnn-1]);
fprintf(DBG,"dmcsim.y[%d]=%12.6e\n",(dbin+dnn-1),dmcsim.y[dbin+dnn-1]);
fprintf(DBG,"ddumnint=%12.6e, ddumcint=%12.6e, ddumposn=%12.6e, dnn=%d\n",\
		ddumnint,ddumcint,ddumposn,dnn);
fflush(DBG);
}
				break;
			case 2:
				if(dnslope!=0) dwidth=ddumnint/dnslope; 
				else{
					if(dcslope!=0) dwidth=ddumcint/dcslope; 
					else dwidth=0;
					}
				dmnsim.y[dbin+dnn]+=ddumnint*0.5*dwidth;
				dmcsim.y[dbin+dnn]+=ddumcint*0.5*dwidth;
				ddoneflag=0;
if(ddebugflag<DEBUG){
fprintf(DBG,"Case 2. dwidth=%12.6e, dmnsim.y[%d]=%12.6e, ",dwidth,\
	(dbin+dnn),dmnsim.y[dbin+dnn]);
fprintf(DBG,"dmcsim.y[%d]=%12.6e\n",(dbin+dnn),dmcsim.y[dbin+dnn]);
fprintf(DBG,"ddumnint=%12.6e, ddumcint=%12.6e, ddumposn=%12.6e, dnn=%d\n",\
		ddumnint,ddumcint,ddumposn,dnn);
fflush(DBG);
}
				break;
			case 3:
				dmnsim.y[dbin+dnn]+=ddumnint*0.5*PSTEP;
				dmcsim.y[dbin+dnn]+=ddumcint*0.5*PSTEP;
				ddoneflag=0;
if(ddebugflag<DEBUG){
fprintf(DBG,"Case 3. dmnsim.y[%d]=%12.6e, ",(dbin+dnn),dmnsim.y[dbin+dnn]);
fprintf(DBG,"dmcsim.y[%d]=%12.6e\n",(dbin+dnn),dmcsim.y[dbin+dnn]);
fprintf(DBG,"ddumnint=%12.6e, ddumcint=%12.6e, ddumposn=%12.6e, dnn=%d\n",\
		ddumnint,ddumcint,ddumposn,dnn);
fflush(DBG);
}
				break; 
/* the following for equipment that jumps from point to point */
			case 4:
				ddumnint-=dnslope*PSTEP; 
				ddumcint-=dcslope*PSTEP; 
				dmnsim.y[dbin+dnn]+=ddumnint;
				dmcsim.y[dbin+dnn]+=ddumcint;
				ddumposn+=PSTEP;
				dnn++;
				if((dnn+dbin)==dn) ddoneflag=0;
if(ddebugflag<DEBUG){
fprintf(DBG,"Case 4. dmnsim.y[%d]=%12.6e, dmcsim.y[%d]=%12.6e\n",\
	(dbin+dnn-1),dmnsim.y[dbin+dnn-1],(dbin+dnn-1),dmcsim.y[dbin+dnn-1]);
fprintf(DBG,"ddumnint=%12.6e, ddumcint=%12.6e, ddumposn=%12.6e, dnn=%d\n",\
		ddumnint,ddumcint,ddumposn,dnn);
fflush(DBG);
}
				break; 
			case 5:
				dwidth=ddend-ddumposn; 
				ddumnint-=dnslope*dwidth; 
				ddumcint-=dcslope*dwidth; 
				dmnsim.y[dbin+dnn]+=ddumnint;
				dmcsim.y[dbin+dnn]+=ddumcint;
				ddumposn=ddend;
				dnn++;
				if((dnn+dbin)==dn) ddoneflag=0;
if(ddebugflag<DEBUG){
fprintf(DBG,"Case 5. dwidth=%12.6e, dmnsim.y[%d]=%12.6e, ",dwidth,\
	(dbin+dnn-1),dmnsim.y[dbin+dnn-1]);
fprintf(DBG,"dmcsim.y[%d]=%12.6e\n",(dbin+dnn-1),dmcsim.y[dbin+dnn-1]);
fprintf(DBG,"ddumnint=%12.6e, ddumcint=%12.6e, ddumposn=%12.6e, dnn=%d\n",\
		ddumnint,ddumcint,ddumposn,dnn);
fflush(DBG);
}
				break; 
			case 6:
				ddoneflag=0;
if(ddebugflag<DEBUG){
fprintf(DBG,"Case 6\n");
fflush(DBG);
}
				break; 
			case 7:
				dmnsim.y[dbin+dnn]+=ddumnint;
				dmcsim.y[dbin+dnn]+=ddumcint;
				ddoneflag=0;
if(ddebugflag<DEBUG){
fprintf(DBG,"Case 7. dmnsim.y[%d]=%12.6e, ",(dbin+dnn),dmnsim.y[dbin+dnn]);
fprintf(DBG,"dmcsim.y[%d]=%12.6e\n",(dbin+dnn),dmcsim.y[dbin+dnn]);
fprintf(DBG,"ddumnint=%12.6e, ddumcint=%12.6e, ddumposn=%12.6e, dnn=%d\n",\
		ddumnint,ddumcint,ddumposn,dnn);
fflush(DBG);
}
				break; 
/* the following for emissions that provide only "spillover" */
			case 8:
				if(ddumposn==dposn){
					dwidth=ddend-ddumposn; 
					ddumnint-=dnslope*dwidth; 
					ddumcint-=dcslope*dwidth; 
					ddumposn=ddend;
					dnn++;
					if((dnn+dbin)==dn) ddoneflag=0; 
					}
				else{
					if(ddend>(dposn+RESN)){
						printf(\
				"Fix switch case in do_sim at 1b.\n");
						exit(1);
						}
					else{
						ddumnint-=dnslope*PSTEP; 
						ddumcint-=dcslope*PSTEP; 
						ddumposn+=PSTEP;
						dnn++;
						if((dnn+dbin)==dn) ddoneflag=0; 
						}
					}
if(ddebugflag<DEBUG){
fprintf(DBG,"Case 8. dwidth=%12.6e, dbin=%d, ",dwidth,dbin);
fprintf(DBG,"ddumnint=%12.6e, ddumcint=%12.6e, ddumposn=%12.6e, dnn=%d\n",\
		ddumnint,ddumcint,ddumposn,dnn);
fflush(DBG);
}
				break;
			default:
				printf("Go fix do_sim.\n");
				exit(1);
			}/* close switch-case */
		}/* close fill bins above the transition position */
	} /* close if frequency in range condition */
				}/* close if non-zero frequency condition */
				}/* close simset-transitions loop */
			}/* close if include condition */
			}/* close loop over simulations */
		}/* close loop over transitions */
	}/* close loop over molecules */ 

/* start loop over previous simulation and write to file */ 
sprintf(dnfile,"%s_simulations/Mol_native_pop",PREF);
DNFILE=fopen(dnfile,"w");
if(DNFILE==NULL){
	printf("Error opening DNFILE for native pop. Exiting.\n");
	exit(1);
	}
sprintf(dcfile,"%s_simulations/Mol_cascade_pop",PREF);
DCFILE=fopen(dcfile,"w");
if(DCFILE==NULL){
	printf("Error opening DCFILE for cascade pop. Exiting.\n");
	exit(1);
	} 
fprintf(DNFILE,"#  File created by %s. \n",PROGRAM_NAME);
fprintf(DNFILE,"#  Contains simulated rovibronic emissions based on \n");
fprintf(DNFILE,"#      molecular populations specified by the user (only).\n");
fprintf(DNFILE,"#  See file Mol_cascade_pop for emissions due to cascade.\n");
fprintf(DNFILE,"#  See file %s_Sim_All for the whole simulation.\n#\n",PREF);
fprintf(DCFILE,"#  File created by %s. \n",PROGRAM_NAME);
fprintf(DCFILE,"#  Contains simulated rovibronic emissions based on \n");
fprintf(DCFILE,"#      cascade of molecular states (only).\n");
fprintf(DCFILE,"#  See file Mol_native_pop for emissions from molecular.\n");
fprintf(DCFILE,"#      populations specified by the user.\n");
fprintf(DCFILE,"#  See file %s_Sim_All for the whole simulation.\n#\n",PREF); 
for(da=0;da<dn;da++){
	dmnsim.x[da]=MINV+0.5*PSTEP+da*PSTEP;
	dmcsim.x[da]=MINV+0.5*PSTEP+da*PSTEP; 
	fprintf(DNFILE,"%18.12e\t%18.12e\n",dmnsim.x[da],dmnsim.y[da]);
	fprintf(DCFILE,"%18.12e\t%18.12e\n",dmcsim.x[da],dmcsim.y[da]); 
	}
fclose(DNFILE);
fclose(DCFILE);

/**************** Atomic emissions ******************/
for(da=0;da<NUMAT;da++){ /* loop over sets of atomic info */ 
for(db=0;db<AT[da].n;db++){ /* loop through transitions in this set */ 
if(AT[da].x[db]!=0){/* if the current frequency is not equal to zero */ 
/* see if line falls within requested simulation min/max */
if(UNITTYPE==0) dposn=AT[da].x[db]/10;
else dposn=1e8/AT[da].x[db];
ddoflag=0;
if(dposn<(MINV-RESN)) ddoflag=1;
if(dposn>(MAXV+RESN)) ddoflag=1; 
if(ddoflag==0){
/* find bin where line originates, also get native and cascade slopes */
	dbin=(int)(floor((dposn-MINV)/PSTEP));
	dnslope=AT[da].p*AT[da].A[db]*AT[da].pop[db]/RESN; 
/* increment through them and calculate the intensities for each bin -- 
   integrate if SCAN and point if JUMP */
	ddoneflag=1;
	ddumposn=dposn;
	ddumnint=AT[da].p*AT[da].A[db]*AT[da].pop[db]; 
	dnn=0;
if(ddebugflag<DEBUG){
fprintf(DBG,"dbin=%d, dposn=%12.6e, dnslope=%12.6e, ddumnint=%12.6e\n",\
		dbin,dposn,dnslope,ddumnint);
fprintf(DBG,"AT[%d].p=%12.6e, AT[%d].A[%d]=%12.6e, AT[%d].pop[%d]=%12.6e\n",\
	da,AT[da].p,da,db,AT[da].A[db],da,db,AT[da].pop[db]); 
fprintf(DBG,"RESN=%12.6e, MINV=%12.6e, PSTEP=%12.6e\n",RESN,MINV,PSTEP);
fflush(DBG);
} 
	if(dbin<0) ddoneflag=0;
	while(ddoneflag==1){ /* fill bins at/below dbin */ 
		dsw=0;
		ddend=MINV+(dbin-dnn)*PSTEP;
		if(dposn==ddumposn) dsw+=1;
		if(ddend<(dposn-RESN)) dsw+=2;
		if(SCANTYPE==1) dsw+=4;
		if((dbin-dnn)>(dn-1)) dsw=8;
		switch (dsw){ 
/* the following for equipment that scans between "points" */
			case 0:
				dasim.y[dbin-dnn]+=\
					(ddumnint-0.5*dnslope*PSTEP)*PSTEP;
				ddumnint-=dnslope*PSTEP; 
				ddumposn-=PSTEP;
				if(dnn==dbin) ddoneflag=0;
				dnn++;
if(ddebugflag<DEBUG){
fprintf(DBG,"AT/BELOW Case 0. dasim.y[%d]=%12.6e, ",\
		(dbin-dnn+1),dasim.y[dbin-dnn+1]);
fprintf(DBG,"ddumnint=%12.6e, ddumposn=%12.6e, dnn=%d\n",\
		ddumnint,ddumposn,dnn);
fflush(DBG);
}
				break;
			case 1:
				dwidth=ddumposn-ddend; 
				dasim.y[dbin-dnn]+=\
					(ddumnint-0.5*dnslope*dwidth)*dwidth;
				ddumnint-=dnslope*dwidth; 
				ddumposn=ddend;
				if(dnn==dbin) ddoneflag=0;
				dnn++;
if(ddebugflag<DEBUG){
fprintf(DBG,"AT/BELOW Case 1. dwidth=%12.6e, dasim.y[%d]=%12.6e, ",dwidth,\
	(dbin-dnn+1),dasim.y[dbin-dnn+1]);
fprintf(DBG,"ddumnint=%12.6e, ddumposn=%12.6e, dnn=%d\n",\
		ddumnint,ddumposn,dnn);
fflush(DBG);
}
				break;
			case 2:
				if(dnslope!=0) dwidth=ddumnint/dnslope; 
				else{
					printf("dammit (below)\n");
					exit(1);
					}
				dasim.y[dbin-dnn]+=ddumnint*0.5*dwidth;
				ddoneflag=0;
if(ddebugflag<DEBUG){
fprintf(DBG,"AT/BELOW Case 2. dwidth=%12.6e, dasim.y[%d]=%12.6e, ",dwidth,\
	(dbin-dnn),dasim.y[dbin-dnn]);
fprintf(DBG,"ddumnint=%12.6e, ddumposn=%12.6e, dnn=%d\n",\
		ddumnint,ddumposn,dnn);
fflush(DBG);
}
				break;
			case 3:
				dasim.y[dbin-dnn]+=ddumnint*0.5*PSTEP;
				ddoneflag=0;
if(ddebugflag<DEBUG){
fprintf(DBG,"AT/BELOW Case 3. dasim.y[%d]=%12.6e, ",\
		(dbin-dnn),dasim.y[dbin-dnn]);
fprintf(DBG,"ddumnint=%12.6e, ddumposn=%12.6e, dnn=%d\n",\
		ddumnint,ddumposn,dnn);
fflush(DBG);
}
				break; 
/* the following for equipment that jumps from point to point */
			case 4:
				ddumnint-=dnslope*PSTEP; 
				dasim.y[dbin-dnn]+=ddumnint;
				ddumposn-=PSTEP;
				if(dnn==dbin) ddoneflag=0;
				dnn++;
if(ddebugflag<DEBUG){
fprintf(DBG,"AT/BELOW Case 4. dasim.y[%d]=%12.6e, ",\
		(dbin-dnn+1),dmnsim.y[dbin-dnn+1]);
fprintf(DBG,"ddumnint=%12.6e, ddumposn=%12.6e, dnn=%d\n",\
		ddumnint,ddumposn,dnn);
fflush(DBG);
}
				break; 
			case 5:
				dwidth=ddumposn-ddend; 
				ddumnint-=dnslope*dwidth; 
				dasim.y[dbin-dnn]+=ddumnint;
				ddumposn=ddend;
				if(dnn==dbin) ddoneflag=0;
				dnn++;
if(ddebugflag<DEBUG){
fprintf(DBG,"AT/BELOW Case 5. dwidth=%12.6e, dasim.y[%d]=%12.6e, ",dwidth,\
	(dbin-dnn+1),dasim.y[dbin-dnn+1]);
fprintf(DBG,"ddumnint=%12.6e, ddumposn=%12.6e, dnn=%d\n",\
		ddumnint,ddumposn,dnn);
fflush(DBG);
}
				break; 
			case 6:
				ddoneflag=0;
if(ddebugflag<DEBUG){
fprintf(DBG,"AT/BELOW Case 6\n");
fflush(DBG);
}
				break; 
			case 7:
				dasim.y[dbin-dnn]+=ddumnint;
				ddoneflag=0;
if(ddebugflag<DEBUG){
fprintf(DBG,"AT/BELOW Case 7. dasim.y[%d]=%12.6e, ",\
		(dbin-dnn),dasim.y[dbin-dnn]);
fprintf(DBG,"ddumnint=%12.6e, ddumposn=%12.6e, dnn=%d\n",\
		ddumnint,ddumposn,dnn);
fflush(DBG);
}
				break; 
/* the following for emissions that provide only "spillover" */
			case 8:
				if(ddumposn==dposn){
					dwidth=ddumposn-ddend; 
					ddumnint-=dnslope*dwidth; 
					ddumcint-=dcslope*dwidth; 
					ddumposn=ddend;
					if(dnn==dbin) ddoneflag=0;
					dnn++;
					}
				else{
					if(ddend<(dposn-RESN)){
						printf(\
				"Fix switch case in do_sim at 1a.\n");
						exit(1);
						}
					else{
						ddumnint-=dnslope*PSTEP; 
						ddumcint-=dcslope*PSTEP; 
						ddumposn-=PSTEP;
						if(dnn==dbin) ddoneflag=0;
						dnn++;
						}
					}
if(ddebugflag<DEBUG){
fprintf(DBG,"AT/BELOW Case 8. dwidth=%12.6e, dbin=%d ",dwidth,dbin);
fprintf(DBG,"ddumnint=%12.6e, ddumposn=%12.6e, dnn=%d\n",\
		ddumnint,ddumposn,dnn);
fflush(DBG);
}
				break;
			default:
				printf("Go fix do_sim.\n");
				exit(1);
			}/* close switch-case */
		}/* close while for bins below the transition position */
	ddoneflag=1;
	ddumposn=dposn;
	ddumnint=AT[da].p*AT[da].A[db]*AT[da].pop[db]; 
	dnn=0;
	if(dbin>dn-1) ddoneflag=0;
	while(ddoneflag==1){ /* fill bins at/above dbin */
		dsw=0;
		ddend=MINV+(dbin+dnn+1)*PSTEP;
		if(dposn==ddumposn) dsw+=1;
		if(ddend>(dposn+RESN)) dsw+=2;
		if(SCANTYPE==1) dsw+=4;
		if((dbin+dnn)<0) dsw=8;
		switch (dsw){
/* the following for equipment that scans between "points" */
			case 0:
				dasim.y[dbin+dnn]+=\
					(ddumnint-0.5*dnslope*PSTEP)*PSTEP;
				ddumnint-=dnslope*PSTEP; 
				ddumposn+=PSTEP;
				dnn++; 
				if((dnn+dbin)==dn) ddoneflag=0; 
if(ddebugflag<DEBUG){
fprintf(DBG,"AT/ABOVE Case 0. dasim.y[%d]=%12.6e ",\
		(dbin+dnn-1),dasim.y[dbin+dnn-1]);
fprintf(DBG,"ddumnint=%12.6e, ddumposn=%12.6e, dnn=%d\n",\
		ddumnint,ddumposn,dnn);
fflush(DBG);
}
				break;
			case 1:
				dwidth=ddend-ddumposn; 
				dasim.y[dbin+dnn]+=\
					(ddumnint-0.5*dnslope*dwidth)*dwidth;
				ddumnint-=dnslope*dwidth; 
				ddumposn=ddend; 
				dnn++; 
				if((dnn+dbin)==dn) ddoneflag=0;
if(ddebugflag<DEBUG){
fprintf(DBG,"AT/ABOVE Case 1. dwidth=%12.6e, dasim.y[%d]=%12.6e, ",dwidth,\
	(dbin+dnn-1),dasim.y[dbin+dnn-1]);
fprintf(DBG,"ddumnint=%12.6e, ddumposn=%12.6e, dnn=%d\n",\
		ddumnint,ddumposn,dnn);
fflush(DBG);
}
				break;
			case 2:
				if(dnslope!=0) dwidth=ddumnint/dnslope; 
				else{
					printf("dammit...\n");
					exit(1);
					}
				dasim.y[dbin+dnn]+=ddumnint*0.5*dwidth;
				ddoneflag=0;
if(ddebugflag<DEBUG){
fprintf(DBG,"AT/ABOVE Case 2. dwidth=%12.6e, dasim.y[%d]=%12.6e, ",dwidth,\
	(dbin+dnn),dasim.y[dbin+dnn]);
fprintf(DBG,"ddumnint=%12.6e, dnslope=%12.6e, ddumposn=%12.6e, dnn=%d\n",\
		ddumnint,dnslope,ddumposn,dnn);
fflush(DBG);
}
				break;
			case 3:
				dasim.y[dbin+dnn]+=ddumnint*0.5*PSTEP;
				ddoneflag=0;
if(ddebugflag<DEBUG){
fprintf(DBG,"AT/ABOVE Case 3. dasim.y[%d]=%12.6e, ",\
		(dbin+dnn),dasim.y[dbin+dnn]);
fprintf(DBG,"ddumnint=%12.6e, ddumposn=%12.6e, dnn=%d\n",\
		ddumnint,ddumposn,dnn);
fflush(DBG);
}
				break; 
/* the following for equipment that jumps from point to point */
			case 4:
				ddumnint-=dnslope*PSTEP; 
				dasim.y[dbin+dnn]+=ddumnint;
				ddumposn+=PSTEP; 
				dnn++; 
				if((dnn+dbin)==dn) ddoneflag=0;
if(ddebugflag<DEBUG){
fprintf(DBG,"AT/ABOVE Case 4. dasim.y[%d]=%12.6e ",\
		(dbin+dnn-1),dasim.y[dbin+dnn-1]);
fprintf(DBG,"ddumnint=%12.6e, ddumposn=%12.6e, dnn=%d\n",\
		ddumnint,ddumposn,dnn);
fflush(DBG);
}
				break; 
			case 5:
				dwidth=ddend-ddumposn; 
				ddumnint-=dnslope*dwidth; 
				dasim.y[dbin+dnn]+=ddumnint;
				ddumposn=ddend; 
				dnn++; 
				if((dnn+dbin)==dn) ddoneflag=0;
if(ddebugflag<DEBUG){
fprintf(DBG,"AT/ABOVE Case 5. dwidth=%12.6e, dasim.y[%d]=%12.6e, ",dwidth,\
	(dbin+dnn-1),dasim.y[dbin+dnn-1]);
fprintf(DBG,"ddumnint=%12.6e, ddumposn=%12.6e, dnn=%d\n",\
		ddumnint,ddumposn,dnn);
fflush(DBG);
}
				break; 
			case 6:
				ddoneflag=0;
if(ddebugflag<DEBUG){
fprintf(DBG,"AT/ABOVE Case 6\n");
fflush(DBG);
}
				break; 
			case 7:
				dasim.y[dbin+dnn]+=ddumnint;
				ddoneflag=0;
if(ddebugflag<DEBUG){
fprintf(DBG,"AT/ABOVE Case 7. dasim.y[%d]=%12.6e, ",\
		(dbin+dnn),dasim.y[dbin+dnn]);
fprintf(DBG,"ddumnint=%12.6e, ddumposn=%12.6e, dnn=%d\n",\
		ddumnint,ddumposn,dnn);
fflush(DBG);
}
				break; 
/* the following for emissions that provide only "spillover" */
			case 8:
				if(ddumposn==dposn){
					dwidth=ddend-ddumposn; 
					ddumnint-=dnslope*dwidth; 
					ddumposn=ddend; 
					dnn++; 
					if((dnn+dbin)==dn) ddoneflag=0; 
					}
				else{
					if(ddend>(dposn+RESN)){
						printf(\
				"Fix switch case in do_sim at 1b.\n");
						exit(1);
						}
					else{
						ddumnint-=dnslope*PSTEP; 
						ddumposn+=PSTEP;
						dnn++;
						if((dnn+dbin)==dn) ddoneflag=0; 
						}
					}
if(ddebugflag<DEBUG){
fprintf(DBG,"AT/ABOVE Case 8. dwidth=%12.6e dbin=%d",dwidth,dbin);
fprintf(DBG,"ddumnint=%12.6e, ddumposn=%12.6e, dnn=%d\n",\
		ddumnint,ddumposn,dnn);
fflush(DBG);
}
				break;
			default:
				printf("Go fix do_sim.\n");
				exit(1);
			}/* close switch-case */
		}/* close fill bins above the transition position */
		} /* close if frequency in range condition */
	}/* close if non-zero frequency condition */
	} /* close atomic transitions loop */
	}/* close atomic sets loop */
/* start loop over previous simulation and write to file */ 
sprintf(dnfile,"%s_simulations/Atomic_lines",PREF);
DNFILE=fopen(dnfile,"w");
if(DNFILE==NULL){
	printf("Error opening DNFILE for native pop. Exiting.\n");
	exit(1);
	}
fprintf(DNFILE,"#  File created by %s. \n",PROGRAM_NAME);
fprintf(DNFILE,"#  Contains simulated spectrum based on \n");
fprintf(DNFILE,"#      atomic information specified by the user (only).\n");
fprintf(DNFILE,"#  See file Mol_cascade_pop for emissions due to cascade.\n");
fprintf(DNFILE,"#  See file Mol_native_pop for emissions from molecular.\n");
fprintf(DNFILE,"#      populations specified by the user.\n");
fprintf(DNFILE,"#  See file %s_Sim_All for the whole simulation.\n#\n",PREF);
for(da=0;da<dn;da++){
	dasim.x[da]=MINV+0.5*PSTEP+da*PSTEP;
	fprintf(DNFILE,"%18.12e\t%18.12e\n",dasim.x[da],dasim.y[da]);
	}
fclose(DNFILE);
/* add all simulations together -- find maximum intensity and adjust */
for(da=0;da<dn;da++){
	dtotsim.y[da]=dmnsim.y[da]+dmcsim.y[da]+dasim.y[da]+dosim.y[da];
	if(dintmax<dtotsim.y[da]) dintmax=dtotsim.y[da];
	}
if(MAXINT>0){
for(da=0;da<dn;da++){
	dtotsim.y[da]*=(MAXINT/dintmax);
	}}
/* write total simulation to file */
/* start loop over previous simulation and write to file */ 
sprintf(dnfile,"%s_simulations/%s_Sim_All",PREF,PREF);
DNFILE=fopen(dnfile,"w");
if(DNFILE==NULL){
	printf("Error opening DNFILE for native pop. Exiting.\n");
	exit(1);
	}
fprintf(DNFILE,"#  File created by %s. \n",PROGRAM_NAME);
fprintf(DNFILE,"#  Contains total simulated spectrum based on \n");
fprintf(DNFILE,"#      all information specified by the user as well as\n");
fprintf(DNFILE,"#      diatomic emissions due to cascade (as possible).\n");
fprintf(DNFILE,"#  See documentation for information about other files.\n");
fprintf(DNFILE,"#\n");
for(da=0;da<dn;da++){
	dtotsim.x[da]=MINV+0.5*PSTEP+da*PSTEP;
	fprintf(DNFILE,"%18.12e\t%18.12e\n",dtotsim.x[da],dtotsim.y[da]);
	}
fclose(DNFILE); 
return;
}

