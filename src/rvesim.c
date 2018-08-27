/* File rvesim.c extracted from RVESIM.c on 20080128 by BLFoley
 * Purpose: main control for rvesim program
 * See changelog for details
 */
#include "rvesim.h"

/************************** main **************************/

/* main opens and reads the input file, assigning values to variables as
	needed and handles other file opening/closing and function calls */
int main(int argc, char *argv[]){

/* generic counters */
int ma=0,mb=0,mc=0,mOnuma=0,mOnumb=0;
/* strings for various things */
char mstrdum[1000],mconf[100]=".rvesimconfig",mtmp[200];
FILE *MCONF;
/* for scaling the transition probabilites to one */
double mtrprobmax=0;

// Set some initial values
DEBUG=NUMTR=NUMAT=NUMOT=NUMEX=NUMMOL=NUMST=0;
/* A note on the DEBUG variable.  If this is equal to -1, no functions
  will print any messages. Adjusting values in .rvesimconfig can cause
 messages about the execution of functions to print.  Use this feature 
 cautiously as it can produce enormous files. */
Case_a_num=Case_b_num=SCANTYPE=UNITTYPE=ROVIB=0;
DETECTTYPE=3;
EFFTYPE=INTERACT=1;
SIMMAX=1;
TEMP=0;
JKCUT=0.001;
MINV=100;
MAXV=3500;
PSTEP=0.5;
MAXINT=5000;
RESN=2.0;
/* debugging-level flags */
mdebugflag=tdebugflag=mrdebugflag=mnadebugflag=mnbdebugflag=msdebugflag=sdebugflag=0;
mtdebugflag=aadebugflag=cacdebugflag=bbdebugflag=cbcdebugflag=bbSSdebugflag=bbSOdebugflag=0;
bbOOdebugflag=ddebugflag=rpdebugflag=radebugflag=rddebugflag=XYdebugflag=HLswitch=0; 
bbSSSATT=bbSOSATT=bbOOSATT=0;
PROGRAM_NAME=strdup("RVESIM");

// Read in the configuration file
MCONF=fopen(mconf,"r");
if(MCONF==NULL){
	printf("Configuration file .rvesimconfig not found.\n");
	printf("Using internal defaults.\n");
	}
else{
	fscanf(MCONF,"%s %s",mstrdum,PROGRAM_NAME);
	fscanf(MCONF,"%s %d",mstrdum,&DEBUG);
	fscanf(MCONF,"%s %d",mstrdum,&mdebugflag);
	fscanf(MCONF,"%s %d",mstrdum,&tdebugflag);
	fscanf(MCONF,"%s %d",mstrdum,&mrdebugflag);
	fscanf(MCONF,"%s %d",mstrdum,&mnadebugflag);
	fscanf(MCONF,"%s %d",mstrdum,&mnbdebugflag);
	fscanf(MCONF,"%s %d",mstrdum,&msdebugflag);
	fscanf(MCONF,"%s %d",mstrdum,&sdebugflag);
	fscanf(MCONF,"%s %d",mstrdum,&mtdebugflag);
	fscanf(MCONF,"%s %d",mstrdum,&cacdebugflag);
	fscanf(MCONF,"%s %d",mstrdum,&cbcdebugflag);
	fscanf(MCONF,"%s %d",mstrdum,&aadebugflag);
	fscanf(MCONF,"%s %d",mstrdum,&bbdebugflag);
	fscanf(MCONF,"%s %d",mstrdum,&bbSSdebugflag);
	fscanf(MCONF,"%s %d",mstrdum,&bbSOdebugflag);
	fscanf(MCONF,"%s %d",mstrdum,&bbOOdebugflag);
	fscanf(MCONF,"%s %d",mstrdum,&ddebugflag);
	fscanf(MCONF,"%s %d",mstrdum,&rpdebugflag);
	fscanf(MCONF,"%s %d",mstrdum,&radebugflag);
	fscanf(MCONF,"%s %d",mstrdum,&rddebugflag);
	fscanf(MCONF,"%s %d",mstrdum,&XYdebugflag);
	fscanf(MCONF,"%s %d",mstrdum,&HLswitch);
	fscanf(MCONF,"%s %lf",mstrdum,&bbSSSATT);
	fscanf(MCONF,"%s %lf",mstrdum,&bbSOSATT);
	fscanf(MCONF,"%s %lf",mstrdum,&bbOOSATT); 
	} 

/* inform user of usage if there is no input file on the command line */
if(argv[1]==NULL){
	printf("Usage:  %s input_file\n",PROGRAM_NAME); 
	printf("See the documentation files for more info.\n");
	exit(1); 
	}
/* otherwise, open the file from the command line */
else{
	INMAIN=fopen(argv[1],"r");
	if(INMAIN==NULL){
		printf("Input file open error.  Exiting.  \n");
		exit(1);
		}
	}
/* read the contents of the main input file */
fscanf(INMAIN,"%s %s",mstrdum,PREF);
if(strcmp(mstrdum,"OUTPREF")!=0){
	printf("Input file error.\n");
	printf("Entry %s should say \"OUTPREF\".  Exiting.\n",mstrdum);
	exit(1);
	}
if(DEBUG>-1){
	sprintf(mstrdum,"%s_debug.txt",PREF);
	DBG=fopen(mstrdum,"w");
	if(DBG==NULL){
		printf("Error opening debugging output file.  Exiting.\n");
		exit(1);
		}
	}
sprintf(mstrdum,"%s_parameter.txt",PREF);
PAR=fopen(mstrdum,"w");
if(PAR==NULL){
	printf("Error opening parameter file.  Exiting.\n");
	exit(1);
	}
strcpy(TMPFILE,".rvesim_temp");
fscanf(INMAIN,"%s %s",mstrdum,mtmp);
if(strcmp(mstrdum,"UNITS")!=0){
	printf("Input file error.\n");
	printf("Entry %s should say \"UNITS\".  Exiting.\n",mstrdum);
	exit(1);
	}
else{
	switch (mtmp[0]){
		case 'c':
			UNITTYPE=1;
			break;
		case 'n':
			UNITTYPE=0;
			break;
		default:
			printf("Unexpected entry for \"Units\". Exiting.\n");
			exit(1);
		}
	}
fscanf(INMAIN,"%s %lf",mstrdum,&RESN);
if(strcmp(mstrdum,"RESOLUTION")!=0){
	printf("Input file error.\n");
	printf("Entry %s should say \"RESOLUTION\".  Exiting.\n",mstrdum);
	exit(1);
	}
else RESN/=2;
fscanf(INMAIN,"%s %s",mstrdum,mtmp);
if(strcmp(mstrdum,"EFFICIENCY")!=0){
	printf("Input file error.\n");
	printf("Entry %s should say \"EFFICIENCY\".  Exiting.\n",mstrdum);
	exit(1);
	}
else{
	switch (mtmp[0]){
		case '1':
			EFFTYPE=1;
			break;
		default:
			EFFTYPE=0;
			EF=(XYlist*)calloc(1,sizeof(XYlist));
			strcpy(EF[0].f.f,mtmp);
			read_XY_file(EF);
		}
	}
fscanf(INMAIN,"%s %s",mstrdum,mtmp);
if(strcmp(mstrdum,"DETECTION")!=0){
	printf("Input file error.\n");
	printf("Entry %s should say \"DETECTION\". Exiting.\n",mstrdum);
	exit(1);
	}
else{
	switch (mtmp[0]){
		case 'c':
			DETECTTYPE=4;
			break;
		case 'p':
			DETECTTYPE=3;
			break;
		default:
		printf("Unexpected entry for \"DETECTION\". Exiting.\n");
			exit(1);
		}
	}
fscanf(INMAIN,"%s %s",mstrdum,mtmp);
if(strcmp(mstrdum,"SCAN")!=0){
	printf("Input file error.\n");
	printf("Entry %s should say \"SCAN\". Exiting.\n",mstrdum);
	exit(1);
	}
else{
	switch (mtmp[0]){
		case 's':
			SCANTYPE=0;
			break;
		case 'j':
			SCANTYPE=1;
			break;
		default:
	printf("Unexpected entry for \"SCAN\". Exiting.\n");
			exit(1);
		}
	} 
fscanf(INMAIN,"%s %lf",mstrdum,&JKCUT);
if(strcmp(mstrdum,"CUTOFF")!=0){
	printf("Input file error.\n");
	printf("Entry %s should say \"CUTOFF\". Exiting.\n",mstrdum);
	exit(1);
	}
if(JKCUT<=0){
fprintf(PAR,"The low intensity cutoff was set to a value too low for\n");
fprintf(PAR,"\tthis program to understand.  Resetting to default 0.001.\n");
	JKCUT=0.001;
	}
fscanf(INMAIN,"%s %s",mstrdum,mtmp);
if(strcmp(mstrdum,"INTERACTIVE")!=0){
	printf("Input file error.\n");
	printf("Entry %s should say \"INTERACTIVE\". Exiting.\n",mstrdum);
	exit(1);
	}
else{
	switch (mtmp[0]){
		case 'i':
			INTERACT=0;
			break;
		case 'n':
			INTERACT=1;
			break;
		default:
	printf("Unexpected entry for \"INTERACTIVE\". Exiting.\n");
			exit(1);
		}
	}
fscanf(INMAIN,"%s %lf",mstrdum,&MINV);
if(strcmp(mstrdum,"MINPOINT")!=0){
	printf("Input file error.\n");
	printf("Entry %s should say \"MINPOINT\". Exiting.\n",mstrdum);
	exit(1);
	}
fscanf(INMAIN,"%s %lf",mstrdum,&MAXV);
if(strcmp(mstrdum,"MAXPOINT")!=0){
	printf("Input file error.\n");
	printf("Entry %s should say \"MAXPOINT\". Exiting.\n",mstrdum);
	exit(1);
	} 
fscanf(INMAIN,"%s %lf",mstrdum,&PSTEP);
if(strcmp(mstrdum,"POINTSTEP")!=0){
	printf("Input file error.\n");
	printf("Entry %s should say \"POINTSTEP\". Exiting.\n",mstrdum);
	exit(1);
	}
MINV-=0.5*PSTEP;
MAXV-=0.5*PSTEP;
fscanf(INMAIN,"%s %lf",mstrdum,&MAXINT);
if(strcmp(mstrdum,"MAXIMUM")!=0){
	printf("Input file error.\n");
	printf("Entry %s should say \"MAXIMUM\". Exiting.\n",mstrdum);
	exit(1);
	}
fscanf(INMAIN,"%s %s",mstrdum,mtmp);
if(strcmp(mstrdum,"ROVIBSIM")!=0){
	printf("Input file error.\n");
	printf("Entry %s should say \"ROVIBSIM\". Exiting.\n",mstrdum);
	exit(1);
	}
else{
	switch (mtmp[0]){
		case 'Y':
			ROVIB=0;
			break;
		case 'N':
			ROVIB=1;
			break;
		case 'y':
			ROVIB=0;
			break;
		case 'n':
			ROVIB=1;
			break;
		default:
	printf("Unexpected entry for \"ROVIBSIM\". Exiting.\n");
			exit(1);
		}
	}
if(ROVIB==0){
	fscanf(INMAIN,"%s %s",INST.f,INTR.f);
	INST.F=fopen(INST.f,"r");
	if(INST.F==NULL){
		printf("Error opening state file. Exiting.\n");
		exit(1);
		}
	fclose(INST.F);
	INTR.F=fopen(INTR.f,"r");
	if(INTR.F==NULL){
		printf("Error opening transition file. Exiting.\n");
		exit(1);
		}
	fclose(INTR.F);
	}
fscanf(INMAIN,"%s %d",mstrdum,&NUMAT);
if(strcmp(mstrdum,"NUMAT")!=0){
	printf("Input file error.\n");
	printf("Entry %s should say \"NUMAT\". Exiting.\n",mstrdum);
	exit(1);
	}
if(NUMAT!=0){
	AT=(Ainfo*)calloc(NUMAT,sizeof(Ainfo));
	for(ma=0;ma<NUMAT;ma++){
		fscanf(INMAIN,"%s %lf",AT[ma].f.f,&AT[ma].p);
		AT[ma].f.F=fopen(AT[ma].f.f,"r");
		if(AT[ma].f.F==NULL){
			printf("Error opening %s. Exiting.\n",AT[ma].f.f);
			exit(1);
			}
		fclose(AT[ma].f.F);
		}
	}
fscanf(INMAIN,"%s %d",mstrdum,&NUMOT);
if(strcmp(mstrdum,"NUMOT")!=0){
	printf("Input file error.\n");
	printf("Entry %s should say \"NUMOT\". Exiting.\n",mstrdum);
	exit(1);
	}
if(NUMOT!=0){
	OT=(XYlist*)calloc(NUMOT,sizeof(XYlist));
	for(ma=0;ma<NUMOT;ma++){
		fscanf(INMAIN,"%s %lf",OT[ma].f.f,&OT[ma].p);
		read_XY_file(&OT[ma]); 
		}
	}
fscanf(INMAIN,"%s %d",mstrdum,&NUMEX);
if(strcmp(mstrdum,"NUMEX")!=0){
	printf("Input file error.\n");
	printf("Entry %s should say \"NUMEX\". Exiting.\n",mstrdum);
	exit(1);
	}
if(NUMEX!=0){
	EX=(XYlist*)calloc(NUMEX,sizeof(XYlist));
	for(ma=0;ma<NUMEX;ma++){
		fscanf(INMAIN,"%s %lf",EX[ma].f.f,&EX[ma].p);
		read_XY_file(&EX[ma]);
		}
	}
// Uncomment if ever absorption is coded in for diatomic transitions
//fscanf(INMAIN,"%s %s",mstrdum,mtmp);
//if(strcmp(mstrdum,"ABSEMISS")!=0){
	//printf("Input file error.\n");
	//printf("Entry %s should say \"ABSEMISS\". Exiting.\n",mstrdum);
	//exit(1);
	//}
//else{
        //switch (mtmp[0]){
                //case 'a':
                //case 'A':
                        //DETECTTYPE=1;
                        //break;
                //default:
			//break;
                //} 
	//}
fclose(INMAIN);

/*  Create directory for output */
sprintf(mstrdum,"mkdir %s_molecules",PREF); 
system(mstrdum); 
/*  FUNCTIONS for reading some of the input files */
if(ROVIB==0) read_state_and_transition_files(); 
if(ROVIB==0) read_FCF_file(); 
if(NUMAT>0) read_atomic_files(); 
/*  Loop through all the transition probabilities and make sure that they
    are all expressed in terms of fractions of one.  The only place where
    this is really necessary is in the calculation of cascade intensities.
    But, it won't hurt elsewhere, either.  If other files containing 
    transition probabilites are added to the program, don't forget to put 
    them in these loops, too.  
 	First:  find the maximum transition probability */
if(mdebugflag<DEBUG){
fprintf(DBG,"NUMMOL=%d, NUMAT=%d\n",NUMMOL,NUMTR);
fflush(DBG);
	}
for(ma=0;ma<NUMMOL;ma++){
	for(mb=0;mb<MOL[ma].trans;mb++){
		mOnuma=MOL[ma].t[mb].Ohi;
		if(mOnuma==0) mOnuma=1;
		mOnumb=MOL[ma].t[mb].Olo;
		if(mOnumb==0) mOnumb=1;
		mOnuma*=mOnumb; 
if(mdebugflag<DEBUG){
fprintf(DBG,"First trprob scan: ma=%d, mb=%d, mOnumb=%d, mOnuma=%d\n",\
		ma,mb,mOnumb,mOnuma);
fflush(DBG);
	}
		for(mc=0;mc<mOnuma;mc++){
			if(MOL[ma].t[mb].P[mc]>mtrprobmax){
				mtrprobmax=MOL[ma].t[mb].P[mc];
if(mdebugflag<DEBUG){
fprintf(DBG,"MOL[%d].t[%d].P[%d]=mtrprobmax=%f\n",ma,mb,mc,mtrprobmax);
fflush(DBG);
	}
				}
			}
		}
	}
for(ma=0;ma<NUMAT;ma++){
	for(mb=0;mb<AT[ma].n;mb++){
		if(AT[ma].A[mb]>mtrprobmax) mtrprobmax=AT[ma].A[mb];
if(mdebugflag<DEBUG){
fprintf(DBG,"AT[%d].A[%d]=%f; mtrprobmax=%f\n",ma,mb,AT[ma].A[mb],mtrprobmax);
fflush(DBG);
	}
		}
	}
/* Now, scale all the transition probabilites to that one */
for(ma=0;ma<NUMMOL;ma++){
	for(mb=0;mb<MOL[ma].trans;mb++){
		mOnuma=MOL[ma].t[mb].Ohi;
		if(mOnuma==0) mOnuma=1;
		mOnumb=MOL[ma].t[mb].Olo;
		if(mOnumb==0) mOnumb=1;
		mOnuma*=mOnumb;
if(mdebugflag<DEBUG){
fprintf(DBG,"Second trprob scan: ma=%d, mb=%d, mOnumb=%d, mOnuma=%d\n",\
		ma,mb,mOnumb,mOnuma);
fflush(DBG);
	}
		for(mc=0;mc<mOnuma;mc++){
if(mdebugflag<DEBUG){
fprintf(DBG,"MOL[%d].t[%d].P[%d]=%f\t",ma,mb,mc,MOL[ma].t[mb].P[mc]);
fflush(DBG);
	}
			MOL[ma].t[mb].P[mc]/=mtrprobmax;
if(mdebugflag<DEBUG){
fprintf(DBG,"MOL[%d].t[%d].P[%d]=%f\n",ma,mb,mc,MOL[ma].t[mb].P[mc]);
fflush(DBG);
	}
				}
			}
		}
for(ma=0;ma<NUMAT;ma++){
	for(mb=0;mb<AT[ma].n;mb++){
if(mdebugflag<DEBUG){
fprintf(DBG,"AT[%d].A[%d]=%f\t",ma,mb,AT[ma].A[mb]);
fflush(DBG);
	}
		AT[ma].A[mb]/=mtrprobmax;
if(mdebugflag<DEBUG){
fprintf(DBG,"AT[%d].A[%d]=%f\n",ma,mb,AT[ma].A[mb]);
fflush(DBG);
	}
		}
	}
/*  FUNCTION to calculate energies and relative intensities for
    rotational states */
if(ROVIB==0) manage_rotations(); 
/*  FUNCTION to step through transitions, figure out what sort of transitions
   they are and then call the right functions to simulate them */
if(ROVIB==0) manage_transitions(); 
/*  FUNCTION to add together all the individual pieces and make the
    simulation -- this function will call others that write the data */
do_simulation(); 
/*  FUNCTION to graph the simulated spectrum, with experimental if 
    requested */
/*  graph_results(); */ 
/*  FUNCTION to ask user if something should be changed and to make the
    change as needed */
/*  ask_for_changes(); */

if(ROVIB==0) fclose(INTR.F); 
if(NUMMOL>0) free(MOL);
if(Case_a_num>0) free(CA);
if(Case_b_num>0) free(CB);
if(NUMAT>0) free(AT);
if(NUMEX>0) free(EX);
if(NUMOT>0) free(OT);
if(EFFTYPE==0) free(EF);
sprintf(mstrdum,"rm %s",TMPFILE);
system(mstrdum);
if(DEBUG>-1) fclose(DBG);
fclose(PAR);
return 0;
}

