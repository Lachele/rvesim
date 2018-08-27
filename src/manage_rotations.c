#include "rvesim.h"

/**************** manage_rotations *****************/

/* This function determines the characteristics of the molecule as they
   relate to rotational levels and then calls another function to calculate
   energy levels and populations.  */
void manage_rotations(){

int mra=0,mrb=0,mrstatetype=0;

if(DEBUG>mrdebugflag){
fprintf(DBG,"\nTop of manage rotations\n\n");
fflush(DBG);
}
/* loop over molecules */
for(mra=0;mra<NUMMOL;mra++){
/* loop over states */
	for(mrb=0;mrb<MOL[mra].states;mrb++){
/* determine type of state */
		if(MOL[mra].s[mrb].Ca>-1){
			if(CA[MOL[mra].s[mrb].Ca].L==0) mrstatetype=0;
			if(CA[MOL[mra].s[mrb].Ca].L>0) mrstatetype=2;
			if(CA[MOL[mra].s[mrb].Ca].S==0) mrstatetype+=0;
			if(CA[MOL[mra].s[mrb].Ca].S>0) mrstatetype+=1;
			}
		if(MOL[mra].s[mrb].Cb>-1){
			if(CB[MOL[mra].s[mrb].Cb].L==0) mrstatetype=0;
			if(CB[MOL[mra].s[mrb].Cb].L>0) mrstatetype=2;
			if(CB[MOL[mra].s[mrb].Cb].S==0) mrstatetype+=0;
			if(CB[MOL[mra].s[mrb].Cb].S>0) mrstatetype+=1;
			}
if(DEBUG>mrdebugflag){
fprintf(DBG,"mrstatetype is %d\n",mrstatetype);
fflush(DBG);
}
/* call appropriate function */
		switch(mrstatetype){
			case 0:
if(DEBUG>mrdebugflag){
fprintf(DBG,"Singlet rotation function called (case 0).\n");
fprintf(DBG,"Mol: %s; State: %s.\n",MOL[mra].Mol,MOL[mra].s[mrb].Name);
fflush(DBG);
}
				singlet_JK(mra,mrb);
if(DEBUG>mrdebugflag){
fprintf(DBG,"Returned from Singlet rotation function (case 0).\n");
fflush(DBG);
} 
				break;
			case 1:
if(DEBUG>mrdebugflag){
fprintf(DBG,"Multiplet L=0 rotation function called.\n");
fprintf(DBG,"Mol: %s; State: %s.\n",MOL[mra].Mol,MOL[mra].s[mrb].Name);
fflush(DBG);
} 
				multiplet_zeroL_JK(mra,mrb); 
				break;
			case 2:
if(DEBUG>mrdebugflag){
fprintf(DBG,"Singlet rotation function called. (case 2)\n");
fprintf(DBG,"Mol: MOL[%d].Mol %s, ",mra,MOL[mra].Mol);
fprintf(DBG,"State: MOL[%d].s[%d].Name %s\n",mra,mrb,MOL[mra].s[mrb].Name);
fflush(DBG);
}
				singlet_JK(mra,mrb);
if(DEBUG>mrdebugflag){
fprintf(DBG,"Returned from Singlet rotation function (case 0).\n");
fflush(DBG);
}
				break;
			case 3:
if(DEBUG>mrdebugflag){
fprintf(DBG,"Multiplet L>0 rotation function called.\n");
fprintf(DBG,"Mol: %s; State: %s.\n",MOL[mra].Mol,MOL[mra].s[mrb].Name);
fflush(DBG);
} 
				if(MOL[mra].s[mrb].Ca>-1){
if(DEBUG>mrdebugflag){
fprintf(DBG,"\tCase a \n");
fflush(DBG);
} 
					multiplet_nonzeroL_JK_a(mra,mrb); 
					}
				if(MOL[mra].s[mrb].Cb>-1){
if(DEBUG>mrdebugflag){
fprintf(DBG,"\tCase b \n");
fflush(DBG);
} 
					multiplet_nonzeroL_JK_b(mra,mrb); 
					} 
				break;
			default:
printf("Error in manage_rot switch case for mra=%d, mrb=%d. Exiting.\n",\
		mra,mrb);
				exit(1);
			} /* close switch case */
if(DEBUG>mrdebugflag){
fprintf(DBG,"\n switch-case just closed mra=%d, mrb=%d\n",mra,mrb);
fflush(DBG);
} 
		} /* close loop over states */
if(DEBUG>mrdebugflag){
fprintf(DBG,"\nstates loop closed mra=%d, mrb=%d\n",mra,mrb);
fflush(DBG);
} 
	} /* close loop over molecules */
if(DEBUG>mrdebugflag){
fprintf(DBG,"\n about to return in manage_rotations\n");
fflush(DBG);
} 
return;
} 

/**************** multiplet_nonzero_JK_a *****************/

/* This function calculates energy levels and rotational distributions
   (if needed) for all multiplet non-Sigma states (Hund's Case a) */
void multiplet_nonzeroL_JK_a(int mnaM, int mnaS){ 
int mnaa=0,mnab=0,mnaJck=1,mnaL=0,mnaV=0,mnaVnum=0,mnaOnum=0;
int mnaCa=0,mnaJDnum=0,mnaD=1,mnaOOD=0,mnaOODend=0,mnaf=0;
double mnaT=0,mnaDv=0,mnaBe=0,mnaAe=0,mnaJMAX=0,mnaJMAXc=0,mnaJMAXf=0;
double mnaaA=0,mnaJcut=0,mnawe=0,mnawexe=0,mnaInt=0,mnaImax=0;
double mnabeta=0,mnaFv=0,mnaTe=0,mnaDe=0,mnaBv=0,mnaTeVib=0;
double mnaSpin=0,mnaJcurr=0,mnaOm=0,mnaBvDv=0;
int mnaJnum=0,mnaJcountnum=0,mnaAlloc=0,mnaO=0,mnaOO=0;
char mnafnstr[1000];
FILE *MNAROT;

if(DEBUG>mnadebugflag){
fprintf(DBG,"\n\nTop of multiplet_nonzeroL_JK_a. Molecule %s; State %s.\n\n",\
		MOL[mnaM].Mol,MOL[mnaM].s[mnaS].Name);
fflush(DBG);
} 
/* Find lowest-level temperature designation */ 
if(MOL[mnaM].s[mnaS].T[0]!=0){
	mnaT=MOL[mnaM].s[mnaS].T[0];
if(DEBUG>mnadebugflag){
fprintf(DBG,"Temperature defined at state level: %f",mnaT);
fflush(DBG);
}
	}
else{
	if(MOL[mnaM].T!=0){
		mnaT=MOL[mnaM].T;
if(DEBUG>mnadebugflag){
fprintf(DBG,"Temperature defined at molecule level: %f",mnaT);
fflush(DBG);
}
		}
	else{
		mnaT=TEMP; 
if(DEBUG>mnadebugflag){
fprintf(DBG,"Temperature defined at global level: %f",mnaT);
fflush(DBG);
}
		}
	}
if(mnaT==0){
	fprintf(PAR,"WARNING!! temperature defined as zero for state ");
	fprintf(PAR,"%s of molecule %s.\n\n",\
		MOL[mnaM].s[mnaS].Name,MOL[mnaM].Mol);
	} 
/* check for state identity sanity and complain if not sane */
if(MOL[mnaM].s[mnaS].Ca==-1){ /* if somehow this isn't Hund's Case a */
	printf("Error in mutiplet_nonzeroL_JK_a.  Any state this function\n");
	printf("calculates should be Hund's Case A.\nState ");
	printf("%s of molecule %s doesn't think it's a Hund's Case a.\n",\
		MOL[mnaM].s[mnaS].Name,MOL[mnaM].Mol);
	printf("Fatal error.  Exiting.\n");
	exit(1);
	}
	mnaCa=MOL[mnaM].s[mnaS].Ca;
	mnaBe=CA[mnaCa].Be;
	mnaAe=CA[mnaCa].ae; 
	mnaTe=CA[mnaCa].Te; 
	mnawe=CA[mnaCa].we; 
	mnawexe=CA[mnaCa].wexe; 
	mnaL=CA[mnaCa].L; 
	mnaSpin=CA[mnaCa].S; 
	mnaOnum=CA[mnaCa].O; 
	if(mnaOnum==0){mnaOnum=(int)(2*mnaSpin+1);}
if(DEBUG>mnadebugflag){
fprintf(DBG,"State is Case (a).\nBe=%f ae=%f ",mnaBe,mnaAe);
fprintf(DBG,"Te=%f; we=%f; wexe=%f;\n",mnaTe,mnawe,mnawexe);
fprintf(DBG,"mnaL=%d; mnaSpin=%f; mnaOnum=%d ",mnaL,mnaSpin,mnaOnum);
fflush(DBG);
}
	
/* Calculate some stuff that will be useful later */
if(mnawe!=0) mnaDe=4*pow(mnaBe,3)/(mnawe*mnawe); /* Centrifugal distortion */
if(DEBUG>mnadebugflag){
fprintf(DBG,"mnaDe is %f;\n",mnaDe);
fflush(DBG);
} 
if(mnaL==0){ /* if this is a Sigma state, there's a problem */
printf("How did multiplet_nonzeroL_JK_a get called if L==0? Exiting.\n");
	exit(1);
	} 
/* Print info to the parameter file */
fprintf(PAR,"STATE %s is a  ",MOL[mnaM].s[mnaS].Name);
if(DEBUG>mnadebugflag){
fprintf(DBG,"STATE %s is a  ",MOL[mnaM].s[mnaS].Name);
fflush(DBG);
}
	switch(mnaL){
		case 1:
if(DEBUG>mnadebugflag){
fprintf(DBG,"Pi state -- Hund's Case (a)\n ");
fflush(DBG);
}
			fprintf(PAR,"Pi state\n");
			break;
		case 2:
if(DEBUG>mnadebugflag){
fprintf(DBG,"Delta state -- Hund's Case (a)\n ");
fflush(DBG);
}
			fprintf(PAR,"Delta state\n");
			break;
		case 3:
if(DEBUG>mnadebugflag){
fprintf(DBG,"Phi state -- Hund's Case (a)\n ");
fflush(DBG);
}
			fprintf(PAR,"Phi state\n");
			break;
		default:
fprintf(PAR,"state with Lambda>3 -- an exotic state\n");
fprintf(PAR,"Interpret the results from this program with care.\n");
			break; 
		}
 
/* There may be multiple Omegas here */ 
if(MOL[mnaM].s[mnaS].Dist!=-1){ 
/* Estimate J max.  This comes from taking the derivative of the 
   population distribution and setting it equal to zero. */ 
	mnaBv=mnaBe-mnaAe/2; 
	mnaaA=mnaBv/(kB*mnaT);
	mnaJMAX=1/sqrt(2*mnaaA) - 0.5;
	mnaJMAXc=ceil(mnaJMAX);
	mnaJMAXf=floor(mnaJMAX);
	mnaaA=mnaBv*mnaJMAXc*(mnaJMAXc+1)/(kB*mnaT);
	mnaImax=(2*mnaJMAXc+1)*exp(-mnaaA); 
	mnaaA=mnaBv*mnaJMAXf*(mnaJMAXf+1)/(kB*mnaT);
	mnaInt=(2*mnaJMAXf+1)*exp(-mnaaA); 
if(DEBUG>mnadebugflag){
fprintf(DBG,"mnaBv=%f; mnaaA=%f; mnaJMAX=%f; mnaJMAXc=%f; mnaJMAXf=%f\n",\
		mnaBv,mnaaA,mnaJMAX,mnaJMAXc,mnaJMAXf);
fprintf(DBG,"mnaImax is %f; mnaInt is %f ---  ",mnaImax,mnaInt);
fflush(DBG);
}
	if((mnaImax)<(mnaInt)){mnaJMAX=mnaJMAXf;}
	else{mnaJMAX=mnaJMAXc;}
if(DEBUG>mnadebugflag){
fprintf(DBG,"mnaJMAX is now %f\n ",mnaJMAX);
fflush(DBG);
} 
/* Find an upper limit for J/K corresponding to the user specification.  To
   find out where this equation comes from, see the documentation, 
   particularly the documentation found in files or directories with the
   word "trick" in the name.  Note that this trick is only good down to
   a cutoff intensity of about 1/10,000 of the maximum value.  */ 
	mnaaA=mnaBv/(kB*mnaT);
	mnaJcut=sqrt((-2/mnaaA)/(JKm + JKb/(log(JKCUT)))) + 1/(2*mnaaA) - 0.5;
	mnaJcut=ceil(mnaJcut);
if(DEBUG>mnadebugflag){
fprintf(DBG,"JKCUT is %f; mnaJcut is %f\n",JKCUT,mnaJcut);
fflush(DBG);
}
/* check this number and chastise programmer if not good...  Also, if not 
   good, find a better number using a less elegant method. */ 
	mnaaA=mnaBv*mnaJcut*(mnaJcut+1)/(kB*mnaT);
	mnaInt=(2*mnaJcut+1)*exp(-mnaaA);
	mnaaA=mnaBv*mnaJMAX*(mnaJMAX+1)/(kB*mnaT);
	mnaImax=(2*mnaJMAX+1)*exp(-mnaaA);
if(DEBUG>mnadebugflag){
fprintf(DBG,"mnaInt=%f;  mnaImax=%f \n ",mnaInt,mnaImax);
fflush(DBG);
}
	if((mnaInt/mnaImax)>JKCUT){
		if(JKCUT>0){ /* put back to 0.0001 */
			printf("Dammit... Go fix the J/K business... \n");
			printf("The calculated mnaJcut is %5.4f ",mnaJcut);
			}
		mnaJck=1;
		while(mnaJck==1){
			mnaJcut++; 
			mnaaA=mnaBv*mnaJcut*(mnaJcut+1)/(kB*mnaT);
			mnaInt=(2*mnaJcut+1)*exp(-mnaaA); 
			if((mnaInt/mnaImax)<JKCUT){
				mnaJck=0;
				}
			} /*vvvvvvvvv put back to 0.0001 vvvvvvvv*/
		if(JKCUT>0) printf("The refined one is %5.4f \n",mnaJcut);
		} 
	}
/* Check for the distribution to be defined in a file. */ 
if(MOL[mnaM].s[mnaS].Dist==-1){ 
	mnaD=0;
	mnaJDnum=(int)(CA[mnaCa].Jmax-CA[mnaCa].Jmin+1); 
	mnaJcut=(int)CA[mnaCa].Jmax;
if(DEBUG>mnadebugflag){
fprintf(DBG,"Distribution is by file.\n ");
fflush(DBG);
}
	} 
/* Allocate memory for J/K array.  Assign s/a and +/- flags. */ 
mnaVnum=MOL[mnaM].s[mnaS].v[0].vnum; /* number of vib levels */
mnaAlloc=mnaVnum*mnaOnum; /* allocate for vib levels and Omegas */
MOL[mnaM].s[mnaS].nV=mnaVnum; 
MOL[mnaM].s[mnaS].nr=mnaAlloc;
MOL[mnaM].s[mnaS].r=(rotset*)calloc(mnaAlloc,sizeof(rotset)); 
MOL[mnaM].s[mnaS].rc=(rotset*)calloc((30*mnaOnum),sizeof(rotset)); 
mnaJnum=mnaJcut+1; /* number of J's to consider */
for(mnaO=0;mnaO<mnaOnum;mnaO++){
	for(mnaa=0;mnaa<mnaVnum;mnaa++){
		mnaOO=mnaO*mnaVnum+mnaa;
		MOL[mnaM].s[mnaS].r[mnaOO].j=mnaJnum;
		MOL[mnaM].s[mnaS].r[mnaOO].J=\
			(double*)calloc((mnaJnum),sizeof(double));
/* There aren't any extra E's here.  If someone ever wants to add in
   a splitting factor for the +/- wavefunctions, then double the space 
   available (mnaJnum) and change the indexing and energy calculation
   in the loops below.*/
		MOL[mnaM].s[mnaS].r[mnaOO].EJ=\
			(double*)calloc((mnaJnum),sizeof(double)); 
		MOL[mnaM].s[mnaS].r[mnaOO].PJ=\
			(double*)calloc((mnaJnum),sizeof(double));
		MOL[mnaM].s[mnaS].r[mnaOO].CJ=\
			(double*)calloc((mnaJnum),sizeof(double)); 
		} 
	}
/* Loop through J/K values, calculating energies and populations 
   (include nuclear effects if applicable). */
for(mnaO=0;mnaO<mnaOnum;mnaO++){
	if(CA[MOL[mnaM].s[mnaS].Ca].O==0){
		mnaOm=mnaL-mnaSpin+mnaO;
		}
	else{
		mnaOm=CA[MOL[mnaM].s[mnaS].Ca].lO[mnaO];
		}
	for(mnaa=0;mnaa<mnaVnum;mnaa++){
		mnaOO=mnaO*mnaVnum+mnaa;
		mnaV=MOL[mnaM].s[mnaS].v[0].vlo+mnaa;
		mnaImax=0;
		mnaTeVib=mnaTe+(mnaV+0.5)*mnawe-(mnaV+0.5)*(mnaV+0.5)*mnawexe; 
if(DEBUG>mnadebugflag){
fprintf(DBG,"mnaVnum is %d; mnaTeVib=%f \n",mnaVnum,mnaTeVib);
fflush(DBG);
}
	mnaJcurr=mnaOm;
	if(mnawe!=0){
		mnaBv=mnaBe-mnaAe*(mnaV+0.5);
		mnaDv=mnaDe+mnabeta*(mnaV+0.5);
		mnaBvDv=mnaBv/(2*mnaDv);
		MOL[mnaM].s[mnaS].r[mnaOO].Jdissoc=\
			floor((sqrt(1+4*mnaBvDv)-1)/2);
		} 
	else MOL[mnaM].s[mnaS].r[mnaOO].Jdissoc=RAND_MAX;
	if(MOL[mnaM].s[mnaS].r[mnaOO].Jdissoc>(mnaOm+mnaJnum)){
		mnaJcountnum=mnaJnum; 
		}
	else{
		mnaJcountnum=(int)(MOL[mnaM].s[mnaS].r[mnaOO].Jdissoc-mnaOm);
fprintf(PAR,"\nWARNING!!  MOLECULE DISSOCIATION (see below)\nAccording to");
fprintf(PAR,"the spectroscopic constants in the state file,\n state");
fprintf(PAR,"%s(%.1f) of Molecule %s in vibration level %d dissociates\n",\
		MOL[mnaM].s[mnaS].Name,mnaOm,MOL[mnaM].Mol,mnaa);
fprintf(PAR,"at rotational level J=%.1f,",MOL[mnaM].s[mnaS].r[mnaOO].Jdissoc);
fprintf(PAR," which is significantly populated at temperature %.1f.\n",mnaT);
fprintf(PAR,"Interpret results from this simulation carefully\n");
		}
	for(mnab=0;mnab<mnaJcountnum;mnab++){
		MOL[mnaM].s[mnaS].r[mnaOO].J[mnab]=mnaJcurr;
		mnaFv=mnaBv*(mnaJcurr*(mnaJcurr+1) - mnaOm*mnaOm)-\
			mnaDv*mnaJcurr*mnaJcurr*(mnaJcurr+1)*(mnaJcurr+1); 
		MOL[mnaM].s[mnaS].r[mnaOO].EJ[mnab]=mnaFv;
		if(mnaD!=0){
			mnaInt=(2*mnaJcurr+1)*exp(-mnaFv/(kB*mnaT)); 
			MOL[mnaM].s[mnaS].r[mnaOO].PJ[mnab]=mnaInt;
			}
if(DEBUG>mnadebugflag){
fprintf(DBG,"mnab is %d; mnaJcurr is %f; mnaFv=%f; mnaInt=%f\n",\
		mnab,mnaJcurr,mnaFv,mnaInt);
fprintf(DBG,"MOL[%d].s[%d].r[%d].EJ[%d]=%f; MOL[%d].s[%d].r[%d].PJ[%d]=%f\n",\
	mnaM,mnaS,mnaa,mnab,MOL[mnaM].s[mnaS].r[mnaOO].EJ[mnab],\
	mnaM,mnaS,mnaa,mnab,MOL[mnaM].s[mnaS].r[mnaOO].PJ[mnab]);
fflush(DBG);
}
		mnaImax+=mnaInt;
if(DEBUG>mnadebugflag){
fprintf(DBG,"\n mnaImax is %f\n\n",mnaImax);
fflush(DBG);
}
		mnaJcurr++;
		}/* close J number loop */
	if(mnaD==0){ /* this means the dist is read from file */
		mnaOOD=(mnaO*mnaVnum+mnaa)*mnaJDnum;
		mnaOODend=(mnaO*mnaVnum+mnaa+1)*mnaJDnum;
		mnaf=(int)(CA[mnaCa].Jmin-mnaOm);
	for(mnab=mnaOOD;mnab<mnaOODend;mnab++){
		MOL[mnaM].s[mnaS].r[mnaa].PJ[mnaf]=\
			CA[mnaCa].Jpop[mnab];
		mnaInt=MOL[mnaM].s[mnaS].r[mnaa].PJ[mnaf];
		mnaImax+=mnaInt;
		mnaf++;
		}
		}
	mnaJcurr=mnaOm;
	for(mnab=0;mnab<mnaJcountnum;mnab++){ 
		MOL[mnaM].s[mnaS].r[mnaOO].EJ[mnab]+=mnaTeVib;
		MOL[mnaM].s[mnaS].r[mnaOO].PJ[mnab]/=mnaImax;
if(DEBUG>mnadebugflag){
fprintf(DBG,"%d\t%f\t%f\t%f\n",mnab,mnaJcurr,\
		MOL[mnaM].s[mnaS].r[mnaOO].EJ[mnab],\
		MOL[mnaM].s[mnaS].r[mnaOO].PJ[mnab]);
fflush(DBG);
}
		mnaJcurr++;
		} /* close J number loop */
if(DEBUG>mnadebugflag){
fprintf(DBG,"\n out of for loop \n");
fflush(DBG);
}
		} /* close vibration number loop */
	} /* close Omega number loop */
sprintf(mnafnstr,"mkdir %s_molecules/%s",PREF,MOL[mnaM].Mol);
system(mnafnstr); 
sprintf(mnafnstr,"%s_molecules/%s/%s_%s_rot.dat",PREF,MOL[mnaM].Mol,\
		MOL[mnaM].Mol,MOL[mnaM].s[mnaS].Name);
MNAROT=fopen(mnafnstr,"w");
if(MNAROT==NULL){
	printf("Error opening output file (singlet_JK) %s.  Exiting.\n",\
			mnafnstr);
	exit(1);
	}
fprintf(MNAROT,"## File generated by program %s.\n",PROGRAM_NAME);
fprintf(MNAROT,"## This file contains ");
fprintf(MNAROT,"rotational state information about\n");
fprintf(MNAROT,"## state %s of molecule %s.\n",\
		MOL[mnaM].s[mnaS].Name,MOL[mnaM].Mol);
fprintf(MNAROT,"##\tTo get the value of J add the Counter to the value\n");
fprintf(MNAROT,"##\tof Omega for that set of columns,\n");
fprintf(MNAROT,"## The columns are, in order:\n## Counter ");
for(mnaO=0;mnaO<mnaOnum;mnaO++){ 
	if(CA[MOL[mnaM].s[mnaS].Ca].O==0){
		fprintf(MNAROT," OMEGA=%.1f<<",(mnaL-mnaSpin+mnaO));
		}
	else{
		fprintf(MNAROT," OMEGA=%.1f<<",\
				CA[MOL[mnaM].s[mnaS].Ca].lO[mnaO]);
		}
	for(mnaa=0;mnaa<mnaVnum;mnaa++){
		mnaV=MOL[mnaM].s[mnaS].v[0].vlo+mnaa;
		fprintf(MNAROT," F(v=%d) Pop(v=%d) ",mnaV,mnaV); 
		}
	fprintf(MNAROT,">> ");
	}
fprintf(MNAROT,"\n");
mnaJcurr=mnaOm;
for(mnab=0;mnab<mnaJcountnum;mnab++){ 
	fprintf(MNAROT,"%7d",mnab);
	for(mnaO=0;mnaO<mnaOnum;mnaO++){
		for(mnaa=0;mnaa<mnaVnum;mnaa++){
			mnaOO=mnaO*mnaVnum+mnaa;
			fprintf(MNAROT,"\t%18.10e\t%18.10e",\
				MOL[mnaM].s[mnaS].r[mnaOO].EJ[mnab],\
				MOL[mnaM].s[mnaS].r[mnaOO].PJ[mnab]);
			}
		}
	fprintf(MNAROT,"\n");
	mnaJcurr++;
	}
fflush(MNAROT);
fclose(MNAROT);
if(DEBUG>mnadebugflag){
fprintf(DBG,"\n out of mnaa(Vnum) loop.  mnaa=%d\n",mnaa);
fflush(DBG);
}
return;
}

/**************** multiplet_nonzero_JK_b *****************/

/* This function calculates energy levels and rotational distributions
   (if needed) for all multiplet non-Sigma states (Hund's Case b) */
void multiplet_nonzeroL_JK_b(int mnbM, int mnbS){ 
int mnba=0,mnbb=0,mnbc=0,mnbKck=1,mnbL=0,mnbV=0,mnbVnum=0;
int mnbCb=0,mnbJDnum=0,mnbJDmin=0,mnbD=1;
double mnbT=0,mnbDv=0,mnbBe=0,mnbAe=0,mnbKMAX=0,mnbKMAXc=0,mnbKMAXf=0;
double mnbaA=0,mnbKcut=0,mnbwe=0,mnbwexe=0,mnbInt=0,mnbImax=0;
double mnbbeta=0,mnbFv=0,mnbTe=0,mnbDe=0,mnbBv=0,mnbTeVib=0;
double mnbSpin=0,mnbJcurr=0,mnbBvDv=0,mnbKcountnum=0;
int mnbKnum=0,mnbJnum=0,mnbcc=0;
char mnbfnstr[1000];
FILE *MNBROT;

if(DEBUG>mnbdebugflag){
fprintf(DBG,"\n\nTop of multiplet_nonzeroL_JK_b. Molecule %s; State %s.\n\n",\
		MOL[mnbM].Mol,MOL[mnbM].s[mnbS].Name);
fflush(DBG);
} 
/* Find lowest-level temperature designation */ 
if(MOL[mnbM].s[mnbS].T[0]!=0){
	mnbT=MOL[mnbM].s[mnbS].T[0];
if(DEBUG>mnbdebugflag){
fprintf(DBG,"Temperature defined at state level: %f",mnbT);
fflush(DBG);
}
	}
else{
	if(MOL[mnbM].T!=0){
		mnbT=MOL[mnbM].T;
if(DEBUG>mnbdebugflag){
fprintf(DBG,"Temperature defined at molecule level: %f",mnbT);
fflush(DBG);
}
		}
	else{
		mnbT=TEMP; 
if(DEBUG>mnbdebugflag){
fprintf(DBG,"Temperature defined at global level: %f",mnbT);
fflush(DBG);
}
		}
	}
if(mnbT==0){
	fprintf(PAR,"WARNING!! temperature defined as zero for state ");
	fprintf(PAR,"%s of molecule %s.\n\n",\
		MOL[mnbM].s[mnbS].Name,MOL[mnbM].Mol);
	} 

/* check for state identity sanity and complain if not sane */
if(MOL[mnbM].s[mnbS].Cb==-1){ /* if somehow this isn't Hund's Case b */
	printf("Error in mutiplet_nonzeroL_JK_b.  Any state this function\n");
	printf("calculates should be Hund's Case B.\nState ");
	printf("%s of molecule %s doesn't think it's a Hund's Case b.\n",\
		MOL[mnbM].s[mnbS].Name,MOL[mnbM].Mol);
	printf("Fatal error.  Exiting.\n");
	exit(1);
	}
	mnbCb=MOL[mnbM].s[mnbS].Cb;
	mnbBe=CB[mnbCb].Be;
	mnbAe=CB[mnbCb].ae; 
	mnbTe=CB[mnbCb].Te; 
	mnbwe=CB[mnbCb].we; 
	mnbwexe=CB[mnbCb].wexe; 
	mnbL=CB[mnbCb].L; 
	mnbSpin=CB[mnbCb].S; 
if(DEBUG>mnbdebugflag){
fprintf(DBG,"State is Case (b).\nBe=%f ae=%f ",mnbBe,mnbAe);
fprintf(DBG,"Te=%f; we=%f; wexe=%f; ",mnbTe,mnbwe,mnbwexe);
fprintf(DBG," mnbL=%d; mnbSpin=%f ",mnbL,mnbSpin);
fflush(DBG);
}
	
/* Calculate some stuff that will be useful later */
if(mnbwe!=0) mnbDe=4*pow(mnbBe,3)/(mnbwe*mnbwe); /* Centrifugal distortion */
if(DEBUG>mnbdebugflag){
fprintf(DBG,"mnbDe is %f;\n",mnbDe);
fflush(DBG);
} 
if(mnbL==0){ /* if this isn't a Sigma state, there's a problem */
printf("How did multiplet_nonzeroL_JK_b get called if L==0? Exiting.\n");
	exit(1);
	} 
/* Print info to the parameter file */
fprintf(PAR,"STATE %s is a  ",MOL[mnbM].s[mnbS].Name);
if(DEBUG>mnbdebugflag){
fprintf(DBG,"STATE %s is a  ",MOL[mnbM].s[mnbS].Name);
fflush(DBG);
}
	switch(mnbL){
		case 1:
if(DEBUG>mnbdebugflag){
fprintf(DBG,"Pi state -- Hund's Case (b)\n ");
fflush(DBG);
}
			fprintf(PAR,"Pi state\n");
			break;
		case 2:
if(DEBUG>mnbdebugflag){
fprintf(DBG,"Delta state -- Hund's Case (b)\n ");
fflush(DBG);
}
			fprintf(PAR,"Delta state\n");
			break;
		case 3:
if(DEBUG>mnbdebugflag){
fprintf(DBG,"Phi state -- Hund's Case (b)\n ");
fflush(DBG);
}
			fprintf(PAR,"Phi state\n");
			break;
		default:
fprintf(PAR,"state with Lambda>3 -- an exotic state\n");
fprintf(PAR,"Interpret the results from this program with care.\n");
			break; 
		}
 
if(MOL[mnbM].s[mnbS].Dist!=-1){ 
/* Estimate J/K max.  This comes from taking the derivative of the 
   population distribution and setting it equal to zero. */ 
	mnbBv=mnbBe-mnbAe/2; 
	mnbaA=mnbBv/(kB*mnbT);
	mnbKMAX=1/sqrt(2*mnbaA) - 0.5;
	mnbKMAXc=ceil(mnbKMAX);
	mnbKMAXf=floor(mnbKMAX);
	mnbaA=mnbBv*mnbKMAXc*(mnbKMAXc+1)/(kB*mnbT);
	mnbImax=(2*mnbKMAXc+1)*exp(-mnbaA); 
	mnbaA=mnbBv*mnbKMAXf*(mnbKMAXf+1)/(kB*mnbT);
	mnbInt=(2*mnbKMAXf+1)*exp(-mnbaA); 
if(DEBUG>mnbdebugflag){
fprintf(DBG,"mnbBv=%f; mnbaA=%f; mnbKMAX=%f; mnbKMAXc=%f; mnbJMAXf=%f\n",\
		mnbBv,mnbaA,mnbKMAX,mnbKMAXc,mnbKMAXf);
fprintf(DBG,"mnbImax is %f; mnbInt is %f ---  ",mnbImax,mnbInt);
fflush(DBG);
}
	if((mnbImax)<(mnbInt)){mnbKMAX=mnbKMAXf;}
	else{mnbKMAX=mnbKMAXc;}
if(DEBUG>mnbdebugflag){
fprintf(DBG,"mnbKMAX is now %f\n ",mnbKMAX);
fflush(DBG);
} 
/* Find an upper limit for J/K corresponding to the user specification.  To
   find out where this equation comes from, see the documentation, 
   particularly the documentation found in files or directories with the
   word "trick" in the name.  Note that this trick is only good down to
   a cutoff intensity of about 1/10,000 of the maximum value.  */ 
	mnbaA=mnbBv/(kB*mnbT);
	mnbKcut=sqrt((-2/mnbaA)/(JKm + JKb/(log(JKCUT)))) + 1/(2*mnbaA) - 0.5;
	mnbKcut=ceil(mnbKcut);
if(DEBUG>mnbdebugflag){
fprintf(DBG,"JKCUT is %f; mnbKcut is %f\n",JKCUT,mnbKcut);
fflush(DBG);
}
/* check this number and chastise programmer if not good...  Also, if not 
   good, find a better number using a less elegant method. */ 
	mnbaA=mnbBv*mnbKcut*(mnbKcut+1)/(kB*mnbT);
	mnbInt=(2*mnbKcut+1)*exp(-mnbaA);
	mnbaA=mnbBv*mnbKMAX*(mnbKMAX+1)/(kB*mnbT);
	mnbImax=(2*mnbKMAX+1)*exp(-mnbaA);
if(DEBUG>mnbdebugflag){
fprintf(DBG,"mnbInt=%f;  mnbImax=%f \n ",mnbInt,mnbImax);
fflush(DBG);
}
	if((mnbInt/mnbImax)>JKCUT){
		if(JKCUT>0){ /* put back to 0.0001 */
			printf("Dammit... Go fix the J/K business... \n");
			printf("The calculated mnbKcut is %5.4f ",mnbKcut);
			}
		mnbKck=1;
		while(mnbKck==1){
			mnbKcut++; 
			mnbaA=mnbBv*mnbKcut*(mnbKcut+1)/(kB*mnbT);
			mnbInt=(2*mnbKcut+1)*exp(-mnbaA); 
			if((mnbInt/mnbImax)<JKCUT){
				mnbKck=0;
				}
			} /*vvvvvvvvv put back to 0.0001 vvvvvvvv*/
		if(JKCUT>0) printf("The refined one is %5.4f \n",mnbKcut);
		} 
	}
/* There are no Omegas here -- check for the distribution
   to be defined in a file. */ 
if(MOL[mnbM].s[mnbS].Dist==-1){ 
	mnbD=0;
	mnbKcut=(int)CB[mnbCb].Kmax;
	mnbJDmin=(int)CB[mnbCb].Kmin;
	mnbJDnum=(int)CB[mnbCb].Jnum; 
if(DEBUG>mnbdebugflag){
fprintf(DBG,"Distribution is by file.\n ");
fflush(DBG);
}
	} 
/* Allocate memory for J/K array.  Assign s/a and +/- flags. */ 
mnbVnum=MOL[mnbM].s[mnbS].v[0].vnum; /* in other functions, where there might be
multiple Omegas, this is vnum*numomega (so far as memory allocation goes) */
MOL[mnbM].s[mnbS].nV=mnbVnum;
MOL[mnbM].s[mnbS].nr=mnbVnum; 
MOL[mnbM].s[mnbS].r=(rotset*)calloc(mnbVnum,sizeof(rotset)); 
MOL[mnbM].s[mnbS].rc=(rotset*)calloc(30,sizeof(rotset)); 
mnbKnum=mnbKcut+1; /* number of K's to consider */
mnbJnum=mnbKnum*((int)(2*(float)mnbSpin+1)); /* number of J's to consider -- 
	this is the multiplcity times the number of K's considered */
for(mnba=0;mnba<mnbVnum;mnba++){
	MOL[mnbM].s[mnbS].r[mnba].j=mnbJnum;
	MOL[mnbM].s[mnbS].r[mnba].k=mnbKnum;
	MOL[mnbM].s[mnbS].r[mnba].J=(double*)calloc((mnbJnum),sizeof(double));
/* The extra E's are redundant here, but if someone ever wants to add in
   a splitting factor for the different J's relative to the nearest K, then
   the space is available */
	MOL[mnbM].s[mnbS].r[mnba].EJ=(double*)calloc((mnbJnum),sizeof(double)); 
	MOL[mnbM].s[mnbS].r[mnba].PJ=(double*)calloc((mnbJnum),sizeof(double));
	MOL[mnbM].s[mnbS].r[mnba].CJ=(double*)calloc((mnbJnum),sizeof(double)); 
	} 
/* Loop through J/K values, calculating energies and populations 
   (include nuclear effects if applicable). */
for(mnba=0;mnba<mnbVnum;mnba++){
	mnbV=MOL[mnbM].s[mnbS].v[0].vlo+mnba;
	mnbImax=0;
	mnbTeVib=mnbTe+(mnbV+0.5)*mnbwe-(mnbV+0.5)*(mnbV+0.5)*mnbwexe; 
	if(mnbwe!=0){
		mnbBv=mnbBe-mnbAe*(mnbV+0.5);
		mnbDv=mnbDe+mnbbeta*(mnbV+0.5);
		mnbBvDv=mnbBv/(2*mnbDv); 
		MOL[mnbM].s[mnbS].r[mnba].Kdissoc=\
			floor((sqrt(1+4*mnbBvDv)-1)/2);
		}
	else MOL[mnbM].s[mnbS].r[mnba].Kdissoc=RAND_MAX;
	if(MOL[mnbM].s[mnbS].r[mnba].Kdissoc>((double)mnbKnum)){
		mnbKcountnum=mnbKnum; 
		} 
	else{
		mnbKcountnum=(int)(MOL[mnbM].s[mnbS].r[mnba].Kdissoc);
fprintf(PAR,"\nWARNING!!  MOLECULE DISSOCIATION (see below)\nAccording to");
fprintf(PAR,"the spectroscopic constants in the state file,\n state");
fprintf(PAR,"%s of Molecule %s in vibration level %d dissociates\n",\
		MOL[mnbM].s[mnbS].Name,MOL[mnbM].Mol,mnba);
fprintf(PAR,"at rotational level K=%.1f,",MOL[mnbM].s[mnbS].r[mnba].Kdissoc);
fprintf(PAR," which is significantly populated at temperature %.1f.\n",mnbT);
fprintf(PAR,"Interpret results from this simulation carefully\n");
		} 
if(DEBUG>mnbdebugflag){
fprintf(DBG,"mnbVnum=%d; mnbTeVib=%f mnbBv=%12.6e, mnbDv=%12.6e\n",\
		mnbVnum,mnbTeVib,mnbBv,mnbDv);
fflush(DBG);
}
	for(mnbb=0;mnbb<mnbKcountnum;mnbb++){
		mnbFv=mnbBv*mnbb*(mnbb+1)-mnbDv*mnbb*mnbb*(mnbb+1)*(mnbb+1); 
		for(mnbc=0;mnbc<((int)(2*mnbSpin+1));mnbc++){
			mnbJcurr=(double)mnbb-mnbSpin+(double)mnbc; 
			mnbcc=mnbb*((int)(2*mnbSpin+1))+mnbc;
			if(mnbD!=0){
				mnbInt=(2*mnbJcurr+1)*exp(-mnbFv/(kB*mnbT)); 
				MOL[mnbM].s[mnbS].r[mnba].PJ[mnbcc]=mnbInt;
				}
			MOL[mnbM].s[mnbS].r[mnba].EJ[mnbcc]=mnbFv; 
			if((mnbb<mnbL)||(mnbJcurr<0)){
				MOL[mnbM].s[mnbS].r[mnba].EJ[mnbcc]=0;
				MOL[mnbM].s[mnbS].r[mnba].PJ[mnbcc]=0;
				}
if(DEBUG>mnbdebugflag){
fprintf(DBG,\
	"mnbb is %d; mnbbeta=%12.6e, mnbJcurr is %f; mnbFv=%f; mnbInt=%f\n",\
		mnbb,mnbbeta,mnbJcurr,mnbFv,mnbInt);
fprintf(DBG,"MOL[%d].s[%d].r[%d].EJ[%d]=%f; MOL[%d].s[%d].r[%d].PJ[%d]=%f\n",\
	mnbM,mnbS,mnba,mnbcc,MOL[mnbM].s[mnbS].r[mnba].EJ[mnbcc],\
	mnbM,mnbS,mnba,mnbcc,MOL[mnbM].s[mnbS].r[mnba].PJ[mnbcc]);
fflush(DBG);
}
			mnbImax+=mnbInt;
if(DEBUG>mnbdebugflag){
fprintf(DBG,"\n mnbImax is %f\n\n",mnbImax);
fflush(DBG);
}
			}
		}
	if(mnbD==0){
	for(mnbb=(mnba*mnbJDnum);mnbb<((mnba+1)*mnbJDnum);mnbb++){
/* the following is algebra for:  K(2S+1) + J - [K-S] */
		mnbcc=(int)(mnbSpin*(2*CB[mnbCb].K[mnbb]+1)+CB[mnbCb].J[mnbb]);
		MOL[mnbM].s[mnbS].r[mnba].PJ[mnbcc]=CB[mnbCb].Jpop[mnbb];
		mnbInt=MOL[mnbM].s[mnbS].r[mnba].PJ[mnbcc];
		mnbImax+=mnbInt;
		}
		}
	for(mnbb=0;mnbb<mnbKcountnum;mnbb++){
		for(mnbc=0;mnbc<((int)(2*mnbSpin+1));mnbc++){
			mnbJcurr=(double)mnbb-mnbSpin+(double)mnbc; 
			mnbcc=mnbb*((int)(2*mnbSpin+1))+mnbc;
			MOL[mnbM].s[mnbS].r[mnba].EJ[mnbcc]+=mnbTeVib;
			MOL[mnbM].s[mnbS].r[mnba].PJ[mnbcc]/=mnbImax;
			if((mnbb<mnbL)||(mnbJcurr<0)){
				MOL[mnbM].s[mnbS].r[mnba].EJ[mnbcc]=0;
				MOL[mnbM].s[mnbS].r[mnba].PJ[mnbcc]=0;
				}
if(DEBUG>mnbdebugflag){
fprintf(DBG,"%d\t%f\t%d\t%f\t%f\n",mnbb,mnbJcurr,mnbcc,\
		MOL[mnbM].s[mnbS].r[mnba].EJ[mnbcc],\
		MOL[mnbM].s[mnbS].r[mnba].PJ[mnbcc]);
fflush(DBG);
}
			}
		} 
if(DEBUG>mnbdebugflag){
fprintf(DBG,"\n out of for loop \n");
fflush(DBG);
}
	} 
sprintf(mnbfnstr,"mkdir %s_molecules/%s",PREF,MOL[mnbM].Mol);
system(mnbfnstr); 
sprintf(mnbfnstr,"%s_molecules/%s/%s_%s_rot.dat",PREF,MOL[mnbM].Mol,\
		MOL[mnbM].Mol,MOL[mnbM].s[mnbS].Name);
MNBROT=fopen(mnbfnstr,"w");
if(MNBROT==NULL){
	printf("Error opening output file (singlet_JK) %s.  Exiting.\n",\
			mnbfnstr);
	exit(1);
	}
fprintf(MNBROT,"## File generated by program %s.\n",PROGRAM_NAME);
fprintf(MNBROT,"## This file contains ");
fprintf(MNBROT,"rotational state information about\n");
fprintf(MNBROT,"## state %s of molecule %s.\n",\
		MOL[mnbM].s[mnbS].Name,MOL[mnbM].Mol);
fprintf(MNBROT,"## The columns are, in order:\n## K  J");
for(mnba=0;mnba<mnbVnum;mnba++){
	mnbV=MOL[mnbM].s[mnbS].v[0].vlo+mnba;
	fprintf(MNBROT,"\tF(v=%d)\tPop(v=%d)",mnbV,mnbV); 
	}
fprintf(MNBROT,"\n");
for(mnbb=0;mnbb<mnbKcountnum;mnbb++){
	for(mnbc=0;mnbc<((int)(2*mnbSpin+1));mnbc++){
		mnbJcurr=(double)mnbb-mnbSpin+(double)mnbc; 
		mnbcc=mnbb*((int)(2*mnbSpin+1))+mnbc;
		fprintf(MNBROT,"%d\t%f",mnbb,mnbJcurr);
		for(mnba=0;mnba<mnbVnum;mnba++){
			fprintf(MNBROT,"\t%18.10e\t%18.10e",\
				MOL[mnbM].s[mnbS].r[mnba].EJ[mnbcc],\
				MOL[mnbM].s[mnbS].r[mnba].PJ[mnbcc]);
			}
		fprintf(MNBROT,"\n");
		}
	}
fflush(MNBROT);
fclose(MNBROT);
if(DEBUG>mnbdebugflag){
fprintf(DBG,"\n out of mnba(Vnum) loop.  mnba=%d\n",mnba);
fflush(DBG);
}
return;
} 

/**************** multiplet_zeroL_JK *****************/

/* This function calculates energy levels and rotational distributions
   (if needed) for all multiplet sigma states (always Hund's Case b) */
void multiplet_zeroL_JK(int msM, int msS){ 
int msa=0,msb=0,msc=0,msKck=1,msg=0,msL=0,msp=0,msV=0,msVnum=0,msItst=0;
int msCb=0,msD=1,msKnum=0,msJnum=0,mscc=0,msJDmin=0,msJDnum=0;
double msT=0,msDv=0,msBe=0,msAe=0,msKMAX=0,msKMAXc=0,msKMAXf=0,msaA=0,msKcut=0;
double mswe=0,mswexe=0,msI=0,msInt=0,msImax=0,msbeta=0,msFv=0,msIs=0;
double msIa=0,msTe=0,msDe=0,msBv=0,msTeVib=0,msSpin=0,msIfac=0,msJcurr=0;
double msBvDv=0;
int msKcountnum=0;
char msfnstr[1000];
FILE *MSROT;

if(DEBUG>msdebugflag){
fprintf(DBG,"\n\nTop of multiplet_zeroL_JK.  Molecule %s;  State %s.\n\n",\
		MOL[msM].Mol,MOL[msM].s[msS].Name);
fflush(DBG);
} 
/* Find lowest-level temperature designation */ 
if(MOL[msM].s[msS].T[0]!=0){
	msT=MOL[msM].s[msS].T[0];
if(DEBUG>msdebugflag){
fprintf(DBG,"Temperature defined at state level: %f",msT);
fflush(DBG);
}
	}
else{
	if(MOL[msM].T!=0){
		msT=MOL[msM].T;
if(DEBUG>msdebugflag){
fprintf(DBG,"Temperature defined at molecule level: %f",msT);
fflush(DBG);
}
		}
	else{
		msT=TEMP; 
if(DEBUG>msdebugflag){
fprintf(DBG,"Temperature defined at global level: %f",msT);
fflush(DBG);
}
		}
	}
if(msT==0){
	fprintf(PAR,"WARNING!! temperature defined as zero for state ");
	fprintf(PAR,"%s of molecule %s.\n\n",\
		MOL[msM].s[msS].Name,MOL[msM].Mol);
	} 

/* check for state identity sanity and complain if not sane */
if(MOL[msM].s[msS].Cb==-1){ /* if somehow this isn't Hund's Case b */
	printf("Error in mutiplet_zeroL_JK.  Any state this function\n");
	printf("calculates should be Hund's Case B.\nState ");
	printf("%s of molecule %s doesn't think it's a Hund's Case b.\n",\
		MOL[msM].s[msS].Name,MOL[msM].Mol);
	printf("Fatal error.  Exiting.\n");
	exit(1);
	}
	msCb=MOL[msM].s[msS].Cb;
	msBe=CB[msCb].Be;
	msAe=CB[msCb].ae; 
	msTe=CB[msCb].Te; 
	mswe=CB[msCb].we; 
	mswexe=CB[msCb].wexe; 
	msg=CB[msCb].g; 
	msI=CB[msCb].I; 
	msL=CB[msCb].L; 
	msp=CB[msCb].p; 
	msSpin=CB[msCb].S; 
if(DEBUG>msdebugflag){
fprintf(DBG,"State is Case (b).\nBe=%f ae=%f ",msBe,msAe);
fprintf(DBG,"Te=%f; we=%f; wexe=%f; \nmsg=%d",msTe,mswe,mswexe,msg);
fprintf(DBG," msI=%f; msL=%d; msp=%d; msSpin=%f ",msI,msL,msp,msSpin);
fflush(DBG);
}
	
/* Calculate some stuff that will be useful later */
if(mswe!=0) msDe=4*pow(msBe,3)/(mswe*mswe); /* Centrifugal distortion */
if(DEBUG>msdebugflag){
fprintf(DBG,"msDe is %f;\n",msDe);
fflush(DBG);
}
if(msg!=0){
	msIs=1;
	msIa=msI/(msI+1);
if(DEBUG>msdebugflag){
fprintf(DBG,"msIs=%f;  msIa=%f; \n",msIs,msIa);
fflush(DBG);
}
	} 
if(msL!=0){ /* if this isn't a Sigma state, there's a problem */
	printf("How did multiplet_zeroL_JK get called if L!=0? Exiting.\n");
	exit(1);
	}
fprintf(PAR,"STATE %s is a  Sigma",MOL[msM].s[msS].Name);
if(DEBUG>msdebugflag){
fprintf(DBG,"STATE %s is a Sigma",MOL[msM].s[msS].Name);
fflush(DBG);
}
if(msg!=0){ /* if this is a homonuclear molecule */
	if(msp==+1){ /* and this is a Sigma+ state */
		fprintf(PAR,"+ ");
if(DEBUG>msdebugflag){
fprintf(DBG,"+ ");
fflush(DBG);
}
		if(msg==+1){ /* and the state is gerade */
		fprintf(PAR,"gerade state and even K's are (+), ");
if(DEBUG>msdebugflag){
fprintf(DBG,"gerade state and even K's are (+), ");
fflush(DBG);
}
			MOL[msM].s[msS].pmsymm=+1; 
			if(fmod(msI,1.0)==0){
				MOL[msM].s[msS].Isymm=+1;
				fprintf(PAR,"s (boson) \n");
if(DEBUG>msdebugflag){
fprintf(DBG,"s (boson)\n");
fflush(DBG);
}
					}
			if(fmod(msI,1.0)==0.5){
				MOL[msM].s[msS].Isymm=-1;
				fprintf(PAR,"a (fermion) \n");
if(DEBUG>msdebugflag){
fprintf(DBG,"a (fermion)\n");
fflush(DBG);
}
				}
			}
		else{ /* and the state is ungerade */
		fprintf(PAR,"ungerade state and even K's are (+), ");
if(DEBUG>msdebugflag){
fprintf(DBG,"ungerade state and even K's are (+), ");
fflush(DBG);
}
			MOL[msM].s[msS].pmsymm=+1; 
			if(fmod(msI,1.0)==0){
				MOL[msM].s[msS].Isymm=-1;
				fprintf(PAR,"a (boson) \n");
if(DEBUG>msdebugflag){
fprintf(DBG,"a (boson) \n");
fflush(DBG);
}
				}
			if(fmod(msI,1.0)==0.5){
				MOL[msM].s[msS].Isymm=+1;
				fprintf(PAR,"s (fermion) \n");
if(DEBUG>msdebugflag){
fprintf(DBG,"s (fermion)\n");
fflush(DBG);
}
				}
			}
		}
	else{ /* and this is a Sigma- state */
		fprintf(PAR,"- ");
if(DEBUG>msdebugflag){
fprintf(DBG,"- ");
fflush(DBG);
}
		if(msg==+1){ /* and the state is gerade */
		fprintf(PAR,"gerade state and even K's are (-), ");
if(DEBUG>msdebugflag){
fprintf(DBG,"gerade state and even K's are (-), ");
fflush(DBG);
}
			MOL[msM].s[msS].pmsymm=-1; 
			if(fmod(msI,1.0)==0){
				MOL[msM].s[msS].Isymm=-1;
				fprintf(PAR,"a (boson) \n");
if(DEBUG>msdebugflag){
fprintf(DBG,"a (boson) \n");
fflush(DBG);
}
				}
			if(fmod(msI,1.0)==0.5){
				MOL[msM].s[msS].Isymm=+1;
				fprintf(PAR,"s (fermion) \n");
if(DEBUG>msdebugflag){
fprintf(DBG,"s (fermion)\n");
fflush(DBG);
}
				}
			}
		else{ /* and the state is ungerade */
		fprintf(PAR,"ungerade state and even K's are (-), ");
if(DEBUG>msdebugflag){
fprintf(DBG,"ungerade state and even K's are (-), ");
fflush(DBG);
}
			MOL[msM].s[msS].pmsymm=-1;
			if(fmod(msI,1.0)==0){
				MOL[msM].s[msS].Isymm=+1;
				fprintf(PAR,"s (boson) \n");
if(DEBUG>msdebugflag){
fprintf(DBG,"s (boson)\n");
fflush(DBG);
}
				}
			if(fmod(msI,1.0)==0.5){
				MOL[msM].s[msS].Isymm=-1;
				fprintf(PAR,"a (fermion) \n");
if(DEBUG>msdebugflag){
fprintf(DBG,"a (fermion)\n");
fflush(DBG);
}
				}
			} 
		}
	}
else{  /* and this a heteronuclear molecule */
	if(msp==+1){
		fprintf(PAR,"+ state and even K's are (+)\n");
		MOL[msM].s[msS].pmsymm=+1; /* and a Sigma+ state */
if(DEBUG>msdebugflag){
fprintf(DBG,"+ state and even K's are (+)\n");
fflush(DBG);
}
		}
	else{
		fprintf(PAR,"- state and even K's are (-)\n");
		MOL[msM].s[msS].pmsymm=-1; /* and a Sigma- state */
if(DEBUG>msdebugflag){
fprintf(DBG,"- state and even K's are (-)\n");
fflush(DBG);
}
		}
	}
if(MOL[msM].s[msS].Dist!=-1){ 
/* Estimate J/K max.  This comes from taking the derivative of the 
   population distribution and setting it equal to zero. */ 
	msBv=msBe-msAe/2; 
	msaA=msBv/(kB*msT);
	msKMAX=1/sqrt(2*msaA) - 0.5;
	msKMAXc=ceil(msKMAX);
	msKMAXf=floor(msKMAX);
	msaA=msBv*msKMAXc*(msKMAXc+1)/(kB*msT);
	msImax=(2*msKMAXc+1)*exp(-msaA); 
	msaA=msBv*msKMAXf*(msKMAXf+1)/(kB*msT);
	msInt=(2*msKMAXf+1)*exp(-msaA); 
if(DEBUG>msdebugflag){
fprintf(DBG,"msBv=%f; msaA=%f; msKMAX=%f; msKMAXc=%f; msJMAXf=%f\n",\
		msBv,msaA,msKMAX,msKMAXc,msKMAXf);
fprintf(DBG,"msImax is %f; msInt is %f ---  ",msImax,msInt);
fflush(DBG);
}
	if((msImax)<(msInt)){msKMAX=msKMAXf;}
	else{msKMAX=msKMAXc;}
if(DEBUG>msdebugflag){
fprintf(DBG,"msKMAX is now %f\n ",msKMAX);
fflush(DBG);
} 
/* Find an upper limit for J/K corresponding to the user specification.  To
   find out where this equation comes from, see the documentation, 
   particularly the documentation found in files or directories with the
   word "trick" in the name.  Note that this trick is only good down to
   a cutoff intensity of about 1/10,000 of the maximum value.  */ 
	msaA=msBv/(kB*msT);
	msKcut=sqrt((-2/msaA)/(JKm + JKb/(log(JKCUT)))) + 1/(2*msaA) - 0.5;
	msKcut=ceil(msKcut);
if(DEBUG>msdebugflag){
fprintf(DBG,"JKCUT is %f; msKcut is %f\n",JKCUT,msKcut);
fflush(DBG);
}
/* check this number and chastise programmer if not good...  Also, if not 
   good, find a better number using a less elegant method. */ 
	msaA=msBv*msKcut*(msKcut+1)/(kB*msT);
	msInt=(2*msKcut+1)*exp(-msaA);
	msaA=msBv*msKMAX*(msKMAX+1)/(kB*msT);
	msImax=(2*msKMAX+1)*exp(-msaA);
if(DEBUG>msdebugflag){
fprintf(DBG,"msInt=%f;  msImax=%f \n ",msInt,msImax);
fflush(DBG);
}
	if((msInt/msImax)>JKCUT){
		if(JKCUT>0){ /* put back to 0.0001 */
			printf("Dammit... Go fix the J/K business... \n");
			printf("The calculated msKcut is %5.4f ",msKcut);
			}
		msKck=1;
		while(msKck==1){
			msKcut++; 
			msaA=msBv*msKcut*(msKcut+1)/(kB*msT);
			msInt=(2*msKcut+1)*exp(-msaA); 
			if((msInt/msImax)<JKCUT){
				msKck=0;
				}
			} /*vvvvvvvvv put back to 0.0001 vvvvvvvv*/
		if(JKCUT>0) printf("The refined one is %5.4f \n",msKcut);
		} 
	}/* close if distribution not in a file condition */
/* There can't be more than one omega here -- check for the distribution
   to be defined in a file.  If it is, go back. */ 
if(MOL[msM].s[msS].Dist==-1){ 
	msD=0;
	msKcut=(int)CB[msCb].Kmax; 
	msJDmin=(int)CB[msCb].Kmin;
	msJDnum=(int)CB[msCb].Jnum;
if(DEBUG>msdebugflag){
fprintf(DBG,"Distribution is by file.\n ");
fflush(DBG);
}
	} 
/* Allocate memory for J/K array.  Assign s/a and +/- flags. */ 
msVnum=MOL[msM].s[msS].v[0].vnum; /* in other functions, where there might be
multiple Omegas, this is vnum*numomega (so far as memory allocation goes) */
MOL[msM].s[msS].nV=msVnum;
MOL[msM].s[msS].nr=msVnum; 
MOL[msM].s[msS].r=(rotset*)calloc(msVnum,sizeof(rotset)); 
MOL[msM].s[msS].rc=(rotset*)calloc(30,sizeof(rotset)); 
msKnum=msKcut+1; /* number of K's to consider */
msJnum=msKnum*((int)(2*(float)msSpin+1)); /* number of J's to consider -- this
     is the multiplcity times the number of K's considered */
for(msa=0;msa<msVnum;msa++){
	MOL[msM].s[msS].r[msa].j=msJnum;
	MOL[msM].s[msS].r[msa].k=msKnum;
	MOL[msM].s[msS].r[msa].J=(double*)calloc((msJnum),sizeof(double));
/* The extra E's are redundant here, but if someone ever wants to add in
   a splitting factor for the different J's relative to the nearest K, then
   the space is available */
	MOL[msM].s[msS].r[msa].EJ=(double*)calloc((msJnum),sizeof(double)); 
	MOL[msM].s[msS].r[msa].PJ=(double*)calloc((msJnum),sizeof(double));
	MOL[msM].s[msS].r[msa].CJ=(double*)calloc((msJnum),sizeof(double)); 
	} 
/* Loop through J/K values, calculating energies and populations 
   (include nuclear effects if applicable). */
for(msa=0;msa<msVnum;msa++){
	msV=MOL[msM].s[msS].v[0].vlo+msa;
	msImax=0;
	msTeVib=msTe+(msV+0.5)*mswe-(msV+0.5)*(msV+0.5)*mswexe; 
if(DEBUG>msdebugflag){
fprintf(DBG,"msVnum is %d; msTeVib=%f \n",msVnum,msTeVib);
fflush(DBG);
}
	if(mswe!=0){
		msBv=msBe-msAe*(msV+0.5);
		msDv=msDe+msbeta*(msV+0.5); 
		msBvDv=msBv/(2*msDv); 
		MOL[msM].s[msS].r[msa].Kdissoc=floor((sqrt(1+4*msBvDv)-1)/2);
		}
	else MOL[msM].s[msS].r[msa].Kdissoc=RAND_MAX;
if(DEBUG>msdebugflag){
fprintf(DBG,"msBv=%12.6e, msDv=%12.6e, msBvDv=%12.6e\n",msBv,msDv,msBvDv);
fprintf(DBG,"msKcountnum=%d, msKnum=%d\n",msKcountnum,msKnum);
fprintf(DBG,"MOL[%d].s[%d].r[%d].Kdissoc=%.1f\n",msM,msS,msa,\
		MOL[msM].s[msS].r[msa].Kdissoc);
fflush(DBG);
}
	if(MOL[msM].s[msS].r[msa].Kdissoc>((double)msKnum)){
		msKcountnum=msKnum; 
		} 
	else{
		msKcountnum=(int)(MOL[msM].s[msS].r[msa].Kdissoc);
fprintf(PAR,"\nWARNING!!  MOLECULE DISSOCIATION (see below)\nAccording to");
fprintf(PAR,"the spectroscopic constants in the state file,\n state");
fprintf(PAR,"%s of Molecule %s in vibration level %d dissociates\n",\
		MOL[msM].s[msS].Name,MOL[msM].Mol,msa);
fprintf(PAR,"at rotational level K=%.1f,",MOL[msM].s[msS].r[msa].Kdissoc);
fprintf(PAR," which is significantly populated at temperature %.1f.\n",msT);
fprintf(PAR,"Interpret results from this simulation carefully\n");
		} 
if(DEBUG>msdebugflag){
fprintf(DBG,"msVnum=%d; msTeVib=%f msBv=%12.6e, msDv=%12.6e\n",\
		msVnum,msTeVib,msBv,msDv);
fprintf(DBG,"msKcountnum=%d, msKnum=%d\n",msKcountnum,msKnum);
fprintf(DBG,"MOL[%d].s[%d].r[%d].Kdissoc=%.1f\n",msM,msS,msa,\
		MOL[msM].s[msS].r[msa].Kdissoc);
fflush(DBG);
}
	for(msb=0;msb<msKcountnum;msb++){ 
		msFv=msBv*msb*(msb+1)-msDv*msb*msb*(msb+1)*(msb+1);
		msIfac=1;
		if((msg!=0)&&(msD!=0)){ 
			msItst=msb%2 + MOL[msM].s[msS].Isymm; 
if(DEBUG>msdebugflag){
fprintf(DBG,"msg is %d\n",msg);
fprintf(DBG,"msItst is %d\n",msItst);
fflush(DBG);
}
			switch(msItst){ /* see "switch_case.txt" in Tricks */
				case 0:
					msIfac=msIs;
if(DEBUG>msdebugflag){
fprintf(DBG,"msIfac=%f",msIfac);
fflush(DBG);
}
					break;
				case 1:
					msIfac=msIs;
if(DEBUG>msdebugflag){
fprintf(DBG,"msIfac=%f",msIfac);
fflush(DBG);
}
					break;
				case -1:
					msIfac=msIa;
if(DEBUG>msdebugflag){
fprintf(DBG,"msIfac=%f",msIfac);
fflush(DBG);
}
					break;
				case 2:
					msIfac=msIa;
if(DEBUG>msdebugflag){
fprintf(DBG,"msIfac=%f",msIfac);
fflush(DBG);
}
					break;
				default:
printf("problem with nuclear stats in singlet_JK. Exiting. \n");
					exit(1); 
				}
			} 
		for(msc=0;msc<((int)(2*msSpin+1));msc++){
			msJcurr=(double)msb-msSpin+(double)msc; 
			mscc=msb*((int)(2*msSpin+1))+msc;
			if(msD!=0){
				msInt=msIfac*(2*msJcurr+1)*exp(-msFv/(kB*msT)); 
if(DEBUG>msdebugflag){
fprintf(DBG,"msIfac=%.1f, msJcurr=%.1f, msFv=%.2f, (kB*msT)=%.1f\n",\
		msIfac,msJcurr,msFv,kB*msT); 
fprintf(DBG,"msD!=0, msInt=%f\n",msInt);
fflush(DBG);
}
				}
			MOL[msM].s[msS].r[msa].EJ[mscc]=msFv;
			if(msD!=0) MOL[msM].s[msS].r[msa].PJ[mscc]=msInt; 
			if(msJcurr<0){
				MOL[msM].s[msS].r[msa].EJ[mscc]=0;
				MOL[msM].s[msS].r[msa].PJ[mscc]=0;
				}
if(DEBUG>msdebugflag){
fprintf(DBG,"msb is %d; msJcurr is %f; msFv=%f; msInt=%f\n",\
		msb,msJcurr,msFv,msInt);
fprintf(DBG,"MOL[%d].s[%d].r[%d].EJ[%d]=%f; MOL[%d].s[%d].r[%d].PJ[%d]=%f\n",\
	msM,msS,msa,mscc,MOL[msM].s[msS].r[msa].EJ[mscc],msM,msS,msa,mscc,\
	MOL[msM].s[msS].r[msa].PJ[mscc]);
fflush(DBG);
}
			msImax+=msInt;
if(DEBUG>msdebugflag){
fprintf(DBG,"\n msImax is %f\n\n",msImax);
fflush(DBG);
}
			}
		}
	if(msD==0){
	for(msb=(msa*msJDnum);msb<((msa+1)*msJDnum);msb++){
/* the following is algebra for:  K(2S+1) + J - [K-S] */
		mscc=(int)(msSpin*(2*CB[msCb].K[msb]+1)+CB[msCb].J[msb]);
		MOL[msM].s[msS].r[msa].PJ[mscc]=CB[msCb].Jpop[msb];
		msInt=MOL[msM].s[msS].r[msa].PJ[mscc];
		msImax+=msInt;
		}
		}
	for(msb=0;msb<msKcountnum;msb++){ 
		for(msc=0;msc<((int)(2*msSpin+1));msc++){
			msJcurr=(double)msb-msSpin+(double)msc; 
			mscc=msb*((int)(2*msSpin+1))+msc;
			MOL[msM].s[msS].r[msa].EJ[mscc]+=msTeVib;
			MOL[msM].s[msS].r[msa].PJ[mscc]/=msImax;
			if(msJcurr<0){
				MOL[msM].s[msS].r[msa].EJ[mscc]=0;
				MOL[msM].s[msS].r[msa].PJ[mscc]=0;
				}
if(DEBUG>msdebugflag){
fprintf(DBG,"%d\t%f\t%d\t%f\t%f\n",msb,msJcurr,mscc,\
		MOL[msM].s[msS].r[msa].EJ[mscc],\
		MOL[msM].s[msS].r[msa].PJ[mscc]);
fflush(DBG);
}
			}
		} 
if(DEBUG>msdebugflag){
fprintf(DBG,"\n out of for loop \n");
fflush(DBG);
}
	} 
sprintf(msfnstr,"mkdir %s_molecules/%s",PREF,MOL[msM].Mol);
system(msfnstr); 
sprintf(msfnstr,"%s_molecules/%s/%s_%s_rot.dat",PREF,MOL[msM].Mol,\
		MOL[msM].Mol,MOL[msM].s[msS].Name);
MSROT=fopen(msfnstr,"w");
if(MSROT==NULL){
	printf("Error opening output file (singlet_JK) %s.  Exiting.\n",\
			msfnstr);
	exit(1);
	}
fprintf(MSROT,"## File generated by program %s.\n",PROGRAM_NAME);
fprintf(MSROT,"## This file contains ");
fprintf(MSROT,"rotational state information about\n");
fprintf(MSROT,"## state %s of molecule %s.\n",\
		MOL[msM].s[msS].Name,MOL[msM].Mol);
fprintf(MSROT,"## The columns are, in order:\n## K  J");
for(msa=0;msa<msVnum;msa++){
	msV=MOL[msM].s[msS].v[0].vlo+msa;
	fprintf(MSROT,"\tF(v=%d)\tPop(v=%d)",msV,msV); 
	}
fprintf(MSROT,"\n");
for(msb=0;msb<msKcountnum;msb++){ 
	for(msc=0;msc<((int)(2*msSpin+1));msc++){
		msJcurr=(double)msb-msSpin+(double)msc; 
		mscc=msb*((int)(2*msSpin+1))+msc;
		fprintf(MSROT,"%d\t%f",msb,msJcurr);
		for(msa=0;msa<msVnum;msa++){
			fprintf(MSROT,"\t%18.10e\t%18.10e",\
				MOL[msM].s[msS].r[msa].EJ[mscc],\
				MOL[msM].s[msS].r[msa].PJ[mscc]);
			}
		fprintf(MSROT,"\n");
		}
	}
fflush(MSROT);
fclose(MSROT);
if(DEBUG>msdebugflag){
fprintf(DBG,"\n out of msa(Vnum) loop.  msa=%d\n",msa);
fflush(DBG);
}
return;
}

/**************** singlet_JK *****************/

/* This function calculates energy levels and rotational distributions
   (if needed) for all singlet states (Hund's Cases a and b) */
void singlet_JK(int sM, int sS){

int sa=0,sb=0,sJck=1,sg=0,sL=0,sp=0,sV=0,sVnum=0,sItst=0;
int sCa=0,sCb=0,sC=0,sD=1,sf=0;
double sT=0,sDv=0,sBe=0,sAe=0,sJMAX=0,sJMAXc=0,sJMAXf=0,saA=0,sJcut=0;
double swe=0,swexe=0,sI=0,sInt=0,sImax=0,sOmm=0,sbeta=0,sFv=0,sBvDv=0;
double sIs=0,sIa=0,sTe=0,sDe=0,sBv=0,sTeVib=0,sJDmin=0,sJcountnum=0;
double sbJ=0;
char sfnstr[1000];
FILE *SROT;

if(DEBUG>sdebugflag){
fprintf(DBG,"\n\nTop of singlet_JK.  Molecule %s;  State %s.\n\n",\
		MOL[sM].Mol,MOL[sM].s[sS].Name);
fflush(DBG);
} 
/* Find lowest-level temperature designation */ 
if(MOL[sM].s[sS].T[0]!=0){
	sT=MOL[sM].s[sS].T[0];
if(DEBUG>sdebugflag){
fprintf(DBG,"Temperature defined at state level: %f",sT);
fflush(DBG);
}
	}
else{
	if(MOL[sM].T!=0){
		sT=MOL[sM].T;
if(DEBUG>sdebugflag){
fprintf(DBG,"Temperature defined at molecule level: %f",sT);
fflush(DBG);
}
		}
	else{
		sT=TEMP; 
if(DEBUG>sdebugflag){
fprintf(DBG,"Temperature defined at global level: %f",sT);
fflush(DBG);
}
		}
	}
if(sT==0){
	fprintf(PAR,"WARNING!! temperature defined as zero for state ");
	fprintf(PAR,"%s of molecule %s.\n\n",MOL[sM].s[sS].Name,MOL[sM].Mol);
	}

/* get the necessary info -- the sOmm variable is there because for singlet
   states, the only difference in the energy functions I'm using here is 
   that for Case a, Lambda occurs in Fv(J).  So, I will always put Lambda
   in the energy term -- it just gets multiplied by zero if the case is b. */ 
if(MOL[sM].s[sS].Ca>-1){ /* if this is Hund's Case a */
	sC=0;
	sCa=MOL[sM].s[sS].Ca;
	sBe=CA[sCa].Be;
	sAe=CA[sCa].ae; 
	sTe=CA[sCa].Te; 
	swe=CA[sCa].we; 
	swexe=CA[sCa].wexe; 
	sg=CA[sCa].g; 
	sI=CA[sCa].I; 
	sL=CA[sCa].L; 
/* This is the place to add in scan for the constant beta */
	sOmm=1;
if(DEBUG>sdebugflag){
fprintf(DBG,"State is Case (a).\n");
fprintf(DBG,"Be=%f; ae=%f; Te=%f; we=%f; wexe=%f; \n",sBe,sAe,sTe,swe,swexe);
fprintf(DBG,"sg=%d;  sI=%f;  sL=%d;  sOmm=%f ",sg,sI,sL,sOmm);
fflush(DBG);
}
	}
if(MOL[sM].s[sS].Cb>-1){ /* if this is Hund's Case b */
	sC=1;
	sCb=MOL[sM].s[sS].Cb;
	sBe=CB[sCb].Be;
	sAe=CB[sCb].ae; 
	sTe=CB[sCb].Te; 
	swe=CB[sCb].we; 
	swexe=CB[sCb].wexe; 
	sg=CB[sCb].g; 
	sI=CB[sCb].I; 
	sL=CB[sCb].L; 
	sp=CB[sCb].p; 
	/* add in scan for beta if that gets included */
	sOmm=0;
if(DEBUG>sdebugflag){
fprintf(DBG,"State is Case (b).\n");
fprintf(DBG,"Be=%f; ae=%f; Te=%f; we=%f; wexe=%f; \n",sBe,sAe,sTe,swe,swexe);
fprintf(DBG,"sg=%d;  sI=%f;  sL=%d; sp=%d;  sOmm=%f ",sg,sI,sL,sp,sOmm);
fflush(DBG);
}
	} 

/* Calculate some stuff that will be useful later */
if(swe!=0) sDe=4*pow(sBe,3)/(swe*swe); /* Centrifugal distortion */
if(DEBUG>sdebugflag){
fprintf(DBG,"sDe is %f;\n",sDe);
fflush(DBG);
}
if(sg!=0){
	sIs=1; /* Equivalent to:  sIs=(2*sI+1)*(sI+1); */
	sIa=sI/(sI+1); /* Equiv. to:  sIa=(2*sI+1)*sI; */
if(DEBUG>sdebugflag){
fprintf(DBG,"sIs=%f;  sIa=%f; \n",sIs,sIa);
fflush(DBG);
}
	}
if(sL==0){ /* if this is a Sigma state */
	fprintf(PAR,"STATE %s is a  Sigma",MOL[sM].s[sS].Name);
if(DEBUG>sdebugflag){
fprintf(DBG,"STATE %s is a Sigma",MOL[sM].s[sS].Name);
fflush(DBG);
}
	if(sg!=0){ /* and this is a homonuclear molecule */
		if(sp==+1){ /* and this is a Sigma+ state */
			fprintf(PAR,"+ ");
if(DEBUG>sdebugflag){
fprintf(DBG,"+ ");
fflush(DBG);
}
			if(sg==+1){ /* and the state is gerade */
			fprintf(PAR,"gerade state and even K's are (+), ");
if(DEBUG>sdebugflag){
fprintf(DBG,"gerade state and even K's are (+), ");
fflush(DBG);
}
				MOL[sM].s[sS].pmsymm=+1; 
				if(fmod(sI,1.0)==0){
					MOL[sM].s[sS].Isymm=+1;
					fprintf(PAR,"s (boson) \n");
if(DEBUG>sdebugflag){
fprintf(DBG,"s (boson)\n");
fflush(DBG);
}
					}
				if(fmod(sI,1.0)==0.5){
					MOL[sM].s[sS].Isymm=-1;
					fprintf(PAR,"a (fermion) \n");
if(DEBUG>sdebugflag){
fprintf(DBG,"a (fermion)\n");
fflush(DBG);
}
					}
				}
			else{ /* and the state is ungerade */
			fprintf(PAR,"ungerade state and even K's are (+), ");
if(DEBUG>sdebugflag){
fprintf(DBG,"ungerade state and even K's are (+), ");
fflush(DBG);
}
				MOL[sM].s[sS].pmsymm=+1; 
				if(fmod(sI,1.0)==0){
					MOL[sM].s[sS].Isymm=-1;
					fprintf(PAR,"a (boson) \n");
if(DEBUG>sdebugflag){
fprintf(DBG,"a (boson) \n");
fflush(DBG);
}
					}
				if(fmod(sI,1.0)==0.5){
					MOL[sM].s[sS].Isymm=+1;
					fprintf(PAR,"s (fermion) \n");
if(DEBUG>sdebugflag){
fprintf(DBG,"s (fermion)\n");
fflush(DBG);
}
					}
				}
			}
		else{ /* and this is a Sigma- state */
			fprintf(PAR,"- ");
if(DEBUG>sdebugflag){
fprintf(DBG,"- ");
fflush(DBG);
}
			if(sg==+1){ /* and the state is gerade */
			fprintf(PAR,"gerade state and even K's are (-), ");
if(DEBUG>sdebugflag){
fprintf(DBG,"gerade state and even K's are (-), ");
fflush(DBG);
}
				MOL[sM].s[sS].pmsymm=-1; 
				if(fmod(sI,1.0)==0){
					MOL[sM].s[sS].Isymm=-1;
					fprintf(PAR,"a (boson) \n");
if(DEBUG>sdebugflag){
fprintf(DBG,"a (boson) \n");
fflush(DBG);
}
					}
				if(fmod(sI,1.0)==0.5){
					MOL[sM].s[sS].Isymm=+1;
					fprintf(PAR,"s (fermion) \n");
if(DEBUG>sdebugflag){
fprintf(DBG,"s (fermion)\n");
fflush(DBG);
}
					}
				}
			else{ /* and the state is ungerade */
			fprintf(PAR,"ungerade state and even K's are (-), ");
if(DEBUG>sdebugflag){
fprintf(DBG,"ungerade state and even K's are (-), ");
fflush(DBG);
}
				MOL[sM].s[sS].pmsymm=-1;
				if(fmod(sI,1.0)==0){
					MOL[sM].s[sS].Isymm=+1;
					fprintf(PAR,"s (boson) \n");
if(DEBUG>sdebugflag){
fprintf(DBG,"s (boson)\n");
fflush(DBG);
}
					}
				if(fmod(sI,1.0)==0.5){
					MOL[sM].s[sS].Isymm=-1;
					fprintf(PAR,"a (fermion) \n");
if(DEBUG>sdebugflag){
fprintf(DBG,"a (fermion)\n");
fflush(DBG);
}
					}
				} 
			}
		}
	else{  /* and this a heteronuclear molecule */
		if(sp==+1){
			fprintf(PAR,"+ state and even K's are (+)\n");
			MOL[sM].s[sS].pmsymm=+1; /* and a Sigma+ state */
if(DEBUG>sdebugflag){
fprintf(DBG,"+ state and even K's are (+)\n");
fflush(DBG);
}
			}
		else{
			fprintf(PAR,"- state and even K's are (-)\n");
			MOL[sM].s[sS].pmsymm=-1; /* and a Sigma- state */
if(DEBUG>sdebugflag){
fprintf(DBG,"- state and even K's are (-)\n");
fflush(DBG);
}
			}
		}
	}
else{ /* if not a Sigma state, print info to the parameter file */
	fprintf(PAR,"STATE %s is a  ",MOL[sM].s[sS].Name);
if(DEBUG>sdebugflag){
fprintf(DBG,"STATE %s is a  ",MOL[sM].s[sS].Name);
fflush(DBG);
}
	switch(sL){
		case 1:
if(DEBUG>sdebugflag){
fprintf(DBG,"Pi state\n ");
fflush(DBG);
}
			fprintf(PAR,"Pi state\n");
			break;
		case 2:
if(DEBUG>sdebugflag){
fprintf(DBG,"Delta state\n ");
fflush(DBG);
}
			fprintf(PAR,"Delta state\n");
			break;
		case 3:
if(DEBUG>sdebugflag){
fprintf(DBG,"Phi state\n ");
fflush(DBG);
}
			fprintf(PAR,"Phi state\n");
			break;
		default:
fprintf(PAR,"state with Lambda>3 -- an exotic state\n");
fprintf(PAR,"Interpret the results from this program with care.\n");
			break; 
		}
	} 

if(MOL[sM].s[sS].Dist!=-1){
/* Estimate J/K max.  This comes from taking the derivative of the 
   population distribution and setting it equal to zero. */ 
	sBv=sBe-sAe/2; 
	saA=sBv/(kB*sT);
	sJMAX=1/sqrt(2*saA) - 0.5;
	sJMAXc=ceil(sJMAX);
	sJMAXf=floor(sJMAX);
/* Note the Omega(=Lambda for singlet states) in the next calculation. 
   The prefactor, called sOmm, is set to zero for Case b states and
   to one for Case a states.  So, the energy should be calculated 
   properly for either state */
	saA=(sBv*(sJMAXc*(sJMAXc+1)-sOmm*sL*sL))/(kB*sT);
	sImax=(2*sJMAXc+1)*exp(-saA); 
	saA=(sBv*(sJMAXf*(sJMAXf+1)-sOmm*sL*sL))/(kB*sT);
	sInt=(2*sJMAXf+1)*exp(-saA); 
if(DEBUG>sdebugflag){
fprintf(DBG,"sBv=%f; saA=%f; sJMAX=%f; sJMAXc=%f; aJMAXf=%f\n",\
		sBv,saA,sJMAX,sJMAXc,sJMAXf);
fprintf(DBG,"sImax is %f; sInt is %f ---  ",sImax,sInt);
fflush(DBG);
}
	if((sImax)<(sInt)){sJMAX=sJMAXf;}
	else{sJMAX=sJMAXc;}
if(DEBUG>sdebugflag){
fprintf(DBG,"sJMAX is now %f\n ",sJMAX);
fflush(DBG);
}

/* Find an upper limit for J/K corresponding to the user specification.  To
   find out where this equation comes from, see the documentation, 
   particularly the documentation found in files or directories with the
   word "trick" in the name.  Note that this trick is only good down to
   a cutoff intensity of about 1/10,000 of the maximum value.  */ 
	saA=sBv/(kB*sT);
	sJcut=sqrt((-2/saA)/(JKm + JKb/(log(JKCUT)))) + 1/(2*saA) - 0.5;
	sJcut=ceil(sJcut);
if(DEBUG>sdebugflag){
fprintf(DBG,"JKCUT is %f; sJcut is %f\n",JKCUT,sJcut);
fflush(DBG);
}
/* check this number and chastise programmer if not good...  Also, if not 
   good, find a better number using a less elegant method. */ 
	saA=(sBv*(sJcut*(sJcut+1)-sOmm*sL*sL))/(kB*sT);
	sInt=(2*sJcut+1)*exp(-saA);
	saA=(sBv*(sJMAX*(sJMAX+1)-sOmm*sL*sL))/(kB*sT);
	sImax=(2*sJMAX+1)*exp(-saA);
if(DEBUG>sdebugflag){
fprintf(DBG,"sInt=%f;  sImax=%f \n ",sInt,sImax);
fflush(DBG);
}
	if((sInt/sImax)>JKCUT){
		if(JKCUT>0){ /* put back to 0.0001 */
			printf("Dammit... Go fix the J/K business... \n");
			printf("The calculated sJcut is %5.4f ",sJcut);
			}
		sJck=1;
		while(sJck==1){
			sJcut++; 
			saA=(sBv*(sJcut*(sJcut+1)-sOmm*sL*sL))/(kB*sT);
			sInt=(2*sJcut+1)*exp(-saA); 
			if((sInt/sImax)<JKCUT){
				sJck=0;
				}
			} /*vvvvvvvvv put back to 0.0001 vvvvvvvv*/
		if(JKCUT>0) printf("The refined one is %5.4f \n",sJcut);
		}
	}

/* There can't be more than one omega here -- check for the distribution
   to be defined in a file.   */ 
if(MOL[sM].s[sS].Dist==-1){ 
	sD=0;
	if(sC==0){
		sJcut=CA[sCa].Jmax;
		sJDmin=CA[sCa].Jmin;
		}
	if(sC==1){
		sJcut=CB[sCb].Kmax;
		sJDmin=CB[sCb].Kmin;
		}
if(DEBUG>sdebugflag){
fprintf(DBG,"Distribution is by file.\n ");
fflush(DBG);
}
	} 
/* Allocate memory for J/K array.  Assign s/a and +/- flags. */ 
sVnum=MOL[sM].s[sS].v[0].vnum; /* in other functions, where there might be
multiple Omegas, this is vnum*numomega (so far as memory allocation goes) */
MOL[sM].s[sS].nV=sVnum;
MOL[sM].s[sS].nr=sVnum; 
MOL[sM].s[sS].r=(rotset*)calloc(sVnum,sizeof(rotset));
MOL[sM].s[sS].rc=(rotset*)calloc(30,sizeof(rotset));
for(sa=0;sa<sVnum;sa++){
	MOL[sM].s[sS].r[sa].k=sJcut+1;
	MOL[sM].s[sS].r[sa].j=sJcut+1;
	MOL[sM].s[sS].r[sa].J=(double*)calloc((sJcut+1),sizeof(double));
	MOL[sM].s[sS].r[sa].EJ=(double*)calloc((sJcut+1),sizeof(double));
	MOL[sM].s[sS].r[sa].PJ=(double*)calloc((sJcut+1),sizeof(double));
	MOL[sM].s[sS].r[sa].CJ=(double*)calloc((sJcut+1),sizeof(double)); 
/* I'm not allocating for K, even for Case b, because for singlet Case
   b states J=K */
	} 
/* Loop through J/K values, calculating energies and populations (include 
   nuclear effects if applicable).  Only one energy calculation is needed 
   here since the energy equation I'm using here is the same for singlet 
   states of Hund's Cases a and b.*/
for(sa=0;sa<sVnum;sa++){
	sV=MOL[sM].s[sS].v[0].vlo+sa;
	sImax=0;
	sTeVib=sTe+(sV+0.5)*swe-(sV+0.5)*(sV+0.5)*swexe; 
	if(swe!=0){
		sBv=sBe-sAe*(sV+0.5);
		sDv=sDe+sbeta*(sV+0.5); 
		sBvDv=sBv/(2*sDv); 
		MOL[sM].s[sS].r[sa].Jdissoc=floor((sqrt(1+4*sBvDv)-1)/2);
		MOL[sM].s[sS].r[sa].Kdissoc=MOL[sM].s[sS].r[sa].Jdissoc;
		}
	else{
		MOL[sM].s[sS].r[sa].Jdissoc=RAND_MAX;
		MOL[sM].s[sS].r[sa].Kdissoc=RAND_MAX;
		}
if(DEBUG>sdebugflag){
fprintf(DBG,"sVnum is %d; sTeVib=%f \n",sVnum,sTeVib);
fprintf(DBG,"sBv=%12.6e, sDv=%12.6e, sBvDv=%12.6e, ",sBv,sDv,sBvDv);
fprintf(DBG,"MOL[%d].s[%d].r[%d].Kdissoc=%12.6e\n",sM,sS,sa,\
		MOL[sM].s[sS].r[sa].Kdissoc);
fflush(DBG);
}
	if(MOL[sM].s[sS].r[sa].Jdissoc>((double)sJcut)){
		sJcountnum=sJcut+1; 
		} 
	else{
		sJcountnum=(int)(MOL[sM].s[sS].r[sa].Jdissoc);
fprintf(PAR,"\nWARNING!!  MOLECULE DISSOCIATION (see below)\nAccording to");
fprintf(PAR,"the spectroscopic constants in the state file,\n state");
fprintf(PAR,"%s of Molecule %s in vibration level %d dissociates\n",\
		MOL[sM].s[sS].Name,MOL[sM].Mol,sa);
fprintf(PAR,"at rotational level J=%.1f,",MOL[sM].s[sS].r[sa].Jdissoc);
fprintf(PAR," which is significantly populated at temperature %.1f.\n",sT);
fprintf(PAR,"Interpret results from this simulation carefully\n");
		} 
if(DEBUG>sdebugflag){
fprintf(DBG,"sVnum=%d; sTeVib=%f sBv=%12.6e, sDv=%12.6e\n",\
		sVnum,sTeVib,sBv,sDv);
fflush(DBG);
}
	for(sb=0;sb<sJcountnum;sb++){ 
		sbJ=sb+sOmm*sL;
		MOL[sM].s[sS].r[sa].J[sb]=sbJ;
		sFv=(sBv*(sbJ*(sbJ+1)-sOmm*sL*sL))-sDv*sbJ*sbJ*(sbJ+1)*(sbJ+1);
		MOL[sM].s[sS].r[sa].EJ[sb]=sFv;
		if(sD!=0){ 
			sInt=(2*sbJ+1)*exp(-sFv/(kB*sT)); 
if(DEBUG>sdebugflag){
fprintf(DBG,"sb is %d;  sFv=%f;  sInt=%f\n",sb,sFv,sInt);
fflush(DBG);
} 
				}
		if((sg!=0)&&(sD!=0)){
/* Note the following bit of trickery.  It works even if the state in
   question has Lambda>0 and doesn't need to have alternating spin-stats.
   The reason for this is that the value of Isymm is zero for those states.
   If Isymm is -1, even K's are a.  If Isymm is +1, even K's are s.  So,
   there are six possibilities:

 	 sb%2-->	0 (K is even)		1 (K is odd) 

	\Isymm/		

	  -1		-1 (a)			0 (s)

	  +1		+1 (s)			2 (a)

	   0		 0 (s-equiv.)		1 (s-equiv)

    The factor for the a state is relative to the weight of the s state 
    (weight with s = 1).  So, only when the result of the calculation is 
    -1 or 2 does sIa need to be multiplied. (but I'm lazy so the switch 
    still has all the cases in it...) */
			sItst=sb%2 + MOL[sM].s[sS].Isymm;
if(DEBUG>sdebugflag){
fprintf(DBG,"sg is 0\n");
fflush(DBG);
}
			switch(sItst){ 
				case 0:
					sInt*=sIs;
if(DEBUG>sdebugflag){
fprintf(DBG,"sInt=%f",sInt);
fflush(DBG);
}
					break;
				case 1:
					sInt*=sIs;
if(DEBUG>sdebugflag){
fprintf(DBG,"sInt=%f",sInt);
fflush(DBG);
}
					break; 
				case -1:
					sInt*=sIa;
if(DEBUG>sdebugflag){
fprintf(DBG,"sInt=%f",sInt);
fflush(DBG);
}
					break;
				case 2:
					sInt*=sIa;
if(DEBUG>sdebugflag){
fprintf(DBG,"sInt=%f",sInt);
fflush(DBG);
}
					break;
				default:
printf("problem with nuclear stats in singlet_JK. Exiting. \n");
					exit(1);
				}
			}
		if((sD==0)&&(sb>=sJDmin)){ 
			if(sC==0){
				sf=(CA[sCa].Jmax-CA[sCa].Jmin+1)*sa+sb-sJDmin;
				MOL[sM].s[sS].r[sa].PJ[sb]=\
					CA[sCa].Jpop[sf];
				}
			if(sC==1){
				sf=(CB[sCb].Jnum)*sa+sb-sJDmin;
				MOL[sM].s[sS].r[sa].PJ[sb]=\
					CB[sCb].Jpop[sf];
				}
			sInt=MOL[sM].s[sS].r[sa].PJ[sb];
if(DEBUG>sdebugflag){
fprintf(DBG,"sb is %d;  sFv=%f;  sInt=%f\n",sb,sFv,sInt);
fflush(DBG);
} 
			}
		if(sD!=0) MOL[sM].s[sS].r[sa].PJ[sb]=sInt;
		sImax+=sInt;
if(DEBUG>sdebugflag){
fprintf(DBG,"\n sImax is %f\n\n",sImax);
fflush(DBG);
}
		}
	for(sb=0;sb<sJcountnum;sb++){ 
		MOL[sM].s[sS].r[sa].EJ[sb]+=sTeVib;
		MOL[sM].s[sS].r[sa].PJ[sb]/=sImax;
if(DEBUG>sdebugflag){
fprintf(DBG,"%d\t%f\t%f\n",sb,MOL[sM].s[sS].r[sa].EJ[sb],\
		MOL[sM].s[sS].r[sa].PJ[sb]);
fflush(DBG);
}
		} 
if(DEBUG>sdebugflag){
fprintf(DBG,"\n out of for loop \n");
fflush(DBG);
}
	} 
sprintf(sfnstr,"mkdir %s_molecules/%s",PREF,MOL[sM].Mol);
system(sfnstr); 
sprintf(sfnstr,"%s_molecules/%s/%s_%s_rot.dat",PREF,MOL[sM].Mol,\
		MOL[sM].Mol,MOL[sM].s[sS].Name);
SROT=fopen(sfnstr,"w");
if(SROT==NULL){
	printf("Error opening output file (singlet_JK) %s.  Exiting.\n",\
			sfnstr);
	exit(1);
	}
fprintf(SROT,"## File generated by program %s.\n",PROGRAM_NAME);
fprintf(SROT,"## This file contains ");
fprintf(SROT,"rotational state information about\n");
fprintf(SROT,"## state %s of molecule %s.\n",MOL[sM].s[sS].Name,MOL[sM].Mol);
fprintf(SROT,"## The columns are, in order:\n## J");
for(sa=0;sa<sVnum;sa++){
	sV=MOL[sM].s[sS].v[0].vlo+sa;
	fprintf(SROT,"\tF(v=%d)\tPop(v=%d)",sV,sV); 
	}
fprintf(SROT,"\n");
for(sb=0;sb<sJcountnum;sb++){ 
	fprintf(SROT,"%d",sb);
	for(sa=0;sa<sVnum;sa++){
		fprintf(SROT,"\t%18.10e\t%18.10e",\
			MOL[sM].s[sS].r[sa].EJ[sb],\
			MOL[sM].s[sS].r[sa].PJ[sb]);
		}
	fprintf(SROT,"\n");
	}
fflush(SROT);
fclose(SROT);
if(DEBUG>sdebugflag){
fprintf(DBG,"\n out of sa(Vnum) loop.  sa=%d\n",sa);
fflush(DBG);
}
return;
} 

