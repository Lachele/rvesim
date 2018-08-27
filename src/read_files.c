#include "rvesim.h"


/****************  new_mol_ck *****************/

/* this function checks entries to see if there is a "MOL" entry where the 
   program doesn't expect one.  If it finds such an entry, it writes a 
   message to stdout and to the parameter file and terminates the program. */
char * new_mol_ck(char *curr_mol, char *curr_st){

char nmol_ck[201]; 
fscanf(INTR.F,"%s",nmol_ck);
if(strcmp("MOL",nmol_ck)==0){
printf("Unexpected MOL entry in state %s of molecule %s.\n",\
	curr_st,curr_mol);
printf("Please edit transition file and restart program.  Exiting.\n");
fprintf(PAR,"Unexpected MOL entry in state %s of molecule %s.\n",\
	curr_st,curr_mol);
fprintf(PAR,"Please edit transition file and restart program.  Exiting.\n");
	exit(1);
	} 
return nmol_ck;
}


/**************** read_FCF_file() *****************/

/* This function reads the FCF file. */
void read_FCF_file(){

int rpC=0,rpa=0,rpb=0,rpc=0,rpd=0,rpz=0,rpy=0,rpw=1,rpv=1,rpvs=0,rpck=0,rpcc=0;
char rpsys[500],rpdum[500],rpfn[400],rpstrtoss[100];
FILE *RPF;
double rptoss=0,rpsumf=0,rpsumfvn=0,rpsumfvc=0,rpcma=0,rpcm=0;
double rpTeH=0,rpTeL=0,rpweH=0,rpweL=0,rpwexeH=0,rpwexeL=0;

sprintf(rpfn,".%s.temp1",PREF);
sprintf(rpdum,".%s.temp2",PREF); 
for(rpa=0;rpa<NUMMOL;rpa++){
	for(rpb=0;rpb<MOL[rpa].trans;rpb++){
/* Get the state and determine the relevant spectroscopic constants.  Also
   check for presence of multiple Omegas <<<  This doesn't work at the
   moment.  Below, only the first set of fcf's is assigned.*/
	if(MOL[rpa].s[MOL[rpa].t[rpb].Hi].Ca!=-1){
		rpC=MOL[rpa].s[MOL[rpa].t[rpb].Hi].Ca;
		rpw=CA[rpC].O;
		MOL[rpa].t[rpb].Ohi=rpw;
		if(rpw<1){rpw=1;} 
		rpTeH=CA[rpC].Te;
		rpweH=CA[rpC].we;
		rpwexeH=CA[rpC].wexe; 
		}
	if(MOL[rpa].s[MOL[rpa].t[rpb].Hi].Cb!=-1){
		rpC=MOL[rpa].s[MOL[rpa].t[rpb].Hi].Cb;
		rpw=1; 
		rpTeH=CB[rpC].Te;
		rpweH=CB[rpC].we;
		rpwexeH=CB[rpC].wexe; 
		}
/* change this if support for Hund's Case C or D is added */	
	if((MOL[rpa].s[MOL[rpa].t[rpb].Hi].Ca==-1)&&\
		(MOL[rpa].s[MOL[rpa].t[rpb].Hi].Cb==-1)){
printf("State %s of molecule ",MOL[rpa].s[MOL[rpa].t[rpb].Hi].Name);
printf("%s is neither Hund's Case a nor b.  Exiting.\n",MOL[rpa].Mol);
		exit(1);
		} 
	if(MOL[rpa].s[MOL[rpa].t[rpb].Lo].Ca!=-1){
		rpC=MOL[rpa].s[MOL[rpa].t[rpb].Lo].Ca;
		rpv=CA[rpC].O;
		MOL[rpa].t[rpb].Olo=rpv;
		if(rpv<1){rpv=1;}
		rpTeL=CA[rpC].Te;
		rpweL=CA[rpC].we;
		rpwexeL=CA[rpC].wexe; 
		}
	if(MOL[rpa].s[MOL[rpa].t[rpb].Lo].Cb!=-1){
		rpC=MOL[rpa].s[MOL[rpa].t[rpb].Lo].Cb;
		rpw=1; 
		rpTeH=CB[rpC].Te;
		rpweH=CB[rpC].we;
		rpwexeH=CB[rpC].wexe; 
		}
/* change this if support for Hund's Case C or D is added */	
	if((MOL[rpa].s[MOL[rpa].t[rpb].Lo].Ca==-1)&&\
		(MOL[rpa].s[MOL[rpa].t[rpb].Lo].Cb==-1)){
printf("State %s of molecule ",MOL[rpa].s[MOL[rpa].t[rpb].Lo].Name);
printf("%s is neither Hund's Case a nor b.  Exiting.\n",MOL[rpa].Mol);
		exit(1);
		} 
MOL[rpa].t[rpb].v=(vfcf*)calloc(30*rpv*rpw,sizeof(vfcf));
if(DEBUG>rpdebugflag){
fprintf(DBG,"rpv=%d, rpw=%d\n",rpv,rpw);
}

for(rpz=0;rpz<rpw;rpz++){
	for(rpy=0;rpy<rpv;rpy++){
		for(rpc=0;rpc<30;rpc++){
			rpcc=rpz*rpv*30+rpy*30+rpc; /* This is the 3D part 
	of a 4D data set.  This (3D) set is indexed so that the primary set
       	is the set of upper-level vibrational numbers.  The next level is
	the number of lower-state Omegas.  There are upper-state-Omega-number
	of sets of (30 upper vib levels)x(lower-state-Omega-number) pieces 
	of data.  later, each one of the lower vib levels gets its own set
	of 30 lower vib levels. */

/* In output from Ervin's program, the second quantum number is that of
   the higher state.  So, this grep is asking for all the entries that
   go to the high state rpc */
			sprintf(rpsys,"grep \"\\->%3d\" %s > %s",\
				rpc,MOL[rpa].t[rpb].f[rpz*rpv+rpy].f,rpfn); 
if(DEBUG>rpdebugflag){
	fprintf(DBG,"grep rpsys is %s\n",rpsys);
}
			system(rpsys);
			sprintf(rpsys,"wc -l %s > %s",rpfn,rpdum);
if(DEBUG>rpdebugflag){
	fprintf(DBG,"wc rpsys is %s\n",rpsys);
}
			system(rpsys);
			RPF=fopen(rpdum,"r");
			if(RPF==NULL){
printf("Error opening temporary file 2 in read FCF for rpa=%d, ",rpa);
printf("rpb=%d, and rpc=%d; filename %s. Exiting.\n",\
		rpb,rpc,MOL[rpa].t[rpb].f[rpz*rpv+rpy].f);
				exit(1);
				}
			fscanf(RPF,"%d",&rpvs);
			fclose(RPF); 
			RPF=fopen(rpfn,"r");
			if(RPF==NULL){
printf("Error opening temporary file 1 in read FCF for rpa=%d, ",rpa);
printf("rpb=%d, and rpc=%d; filename %s. Exiting.\n",\
		rpb,rpc,MOL[rpa].t[rpb].f[rpz*rpv+rpy].f);
				exit(1);
				}
if(DEBUG>rpdebugflag){
	fprintf(DBG,"The following are the entries from the FCF file.\n");
	fprintf(DBG,"Entries that are saved to an array are in brackets[].\n");
	fprintf(DBG,"rpvs is %d\n",rpvs);
}
	rpcma=rpTeH+rpweH*(rpc+0.5)-rpwexeH*(rpc+0.5)*(rpc+0.5)-rpTeL;
			rpsumf=rpsumfvn=rpsumfvc=0;	
			for(rpd=0;rpd<rpvs;rpd++){
		rpcm=rpcma-rpweL*(rpd+0.5)+rpwexeL*(rpd+0.5)*(rpd+0.5);
		fscanf(RPF,"%d",&MOL[rpa].t[rpb].v[rpc].vqn[rpd]); 
				fscanf(RPF,"%s",rpstrtoss); /* comma */
				fscanf(RPF,"%d",&rpck); /* vhi */
if(DEBUG>rpdebugflag){
fprintf(DBG,"[%d] %s %d ",MOL[rpa].t[rpb].v[rpc].vqn[rpd],rpstrtoss,rpck);
}
				if(rpck!=rpc){
printf("rpc not equal to rpck(%d) for rpa=%d, rpb=%d, rpc=%d, rpd=%d\n",\
						rpck,rpa,rpb,rpc,rpd);
				exit(1);
				}
				fscanf(RPF,"%s ",rpstrtoss); /* comma */ 
				fscanf(RPF,"%lf",&rptoss); /* wavenumbers */
if(DEBUG>rpdebugflag){ 
fprintf(DBG,"%s %f ",rpstrtoss,rptoss); 
fflush(DBG);
}
				fscanf(RPF,"%s ",rpstrtoss); /* comma */
		fscanf(RPF,"%lf",&MOL[rpa].t[rpb].v[rpc].fcfn[rpd]); 
				if(rptoss<=0){ /* negative trans. frequency */
					MOL[rpa].t[rpb].v[rpc].fcfn[rpd]=0; 
					} 
				MOL[rpa].t[rpb].v[rpc].fcfc[rpd]=\
					MOL[rpa].t[rpb].v[rpc].fcfn[rpd]; 
				rpsumf+=MOL[rpa].t[rpb].v[rpc].fcfc[rpd];
if(DEBUG>rpdebugflag){
fprintf(DBG,"%s [%f](n)\n",rpstrtoss,MOL[rpa].t[rpb].v[rpc].fcfn[rpd]);
fprintf(DBG,"MOL[%d].t[%d].v[%d].fcfc[%d]=%f\n",rpa,rpb,rpc,rpd,\
		MOL[rpa].t[rpb].v[rpc].fcfc[rpd]);
fflush(DBG);
}
				MOL[rpa].t[rpb].v[rpc].fcfc[rpd]*=pow(rpcm,3);
				MOL[rpa].t[rpb].v[rpc].fcfn[rpd]*=\
					pow(rpcm,DETECTTYPE);
				rpsumfvn+=MOL[rpa].t[rpb].v[rpc].fcfn[rpd];
				rpsumfvc+=MOL[rpa].t[rpb].v[rpc].fcfc[rpd];
if(DEBUG>rpdebugflag){
fprintf(DBG,"rpsumf=%f [%f](n)\n",rpsumf,MOL[rpa].t[rpb].v[rpc].fcfn[rpd]);
fprintf(DBG,"MOL[%d].t[%d].v[%d].fcfc[%d]=%f\n",rpa,rpb,rpc,rpd,\
		MOL[rpa].t[rpb].v[rpc].fcfc[rpd]);
fprintf(DBG,"rpsumfvn=%13.6e ; rpsumfvc=%13.6e\n",rpsumfvn,rpsumfvc);
fflush(DBG);
}
				fscanf(RPF,"%s ",rpstrtoss); /* comma */
				fscanf(RPF,"%lf",&rptoss); 
if(DEBUG>rpdebugflag){ 
fprintf(DBG,"%s %f ",rpstrtoss,rptoss);
fflush(DBG);
}
				fscanf(RPF,"%s ",rpstrtoss); /* comma */ 
				fscanf(RPF,"%lf",&rptoss); /* pec3 */ 
if(DEBUG>rpdebugflag){ 
fprintf(DBG,"%s %f\n",rpstrtoss,rptoss); 
fflush(DBG);
} 
				}  /* close read fcf & nm for loop */ 
/*  This piece of code scales FCF's for the frequency factor.  This isn't the
    best way to do this for two reasons.  One, this should be done at the 
    rot-->rot transition level.  Two, if the franck-condon factors don't all
    add to one, then this becomes an approximation to an approximation. */
			if((rpsumf-1)>0.01){
printf("rpsumf for state %s of molecule",MOL[rpa].s[MOL[rpa].t[rpb].Hi].Name);
printf(" %s is %f (greater than one).  Exiting.\n",MOL[rpa].Mol,rpsumf);
				exit(1);
				}
			if((1-rpsumf)>0.01){
fprintf(PAR,"WARNING!!:  The sum of the Franck-Condon factors (rpsumf) for");
fprintf(PAR,"\n\tstate %s ",MOL[rpa].s[MOL[rpa].t[rpb].Hi].Name);
fprintf(PAR,"of molecule %s is %f (v-hi=%d).\n",MOL[rpa].Mol,rpsumf,rpc); 
				}
			for(rpd=0;rpd<rpvs;rpd++){ 
			MOL[rpa].t[rpb].v[rpc].fcfc[rpd]*=rpsumf/rpsumfvc;
			MOL[rpa].t[rpb].v[rpc].fcfn[rpd]*=rpsumf/rpsumfvn; 
if(DEBUG>rpdebugflag){
fprintf(DBG,"MOL[%d].t[%d].v[%d].fcfn[%d]=%f ; ",rpa,rpb,rpc,rpd,\
		MOL[rpa].t[rpb].v[rpc].fcfn[rpd]);
fprintf(DBG,"MOL[%d].t[%d].v[%d].fcfc[%d]=%f\n",rpa,rpb,rpc,rpd,\
		MOL[rpa].t[rpb].v[rpc].fcfc[rpd]);
fflush(DBG);
}
				}/* close scale fcf for frequency loop */ 
			}  /* close vib levels for loop */
		} /* close rpy (num low Omegas) loop */
	} /* close rpz (num high Omegas) loop */
		}  /* close transitions for loop */
	} /* close NUMMOL for loop */ 

sprintf(rpsys,"rm %s %s",rpfn,rpdum);
system(rpsys);
return;
} 


/**************** read_atomic_file() *****************/

/* This function reads atomic files. */ 
void read_atomic_files(){

int raa=0,rab=0,raz=0;
char rtmp[500];
FILE *ATF;

if(DEBUG>radebugflag){
fprintf(DBG,"Top of read_atomic_files.  NUMAT=%d.\n",NUMAT); 
fflush(DBG);
} 
for(raa=0;raa<NUMAT;raa++){
	sprintf(rtmp,"wc -l %s > %s",AT[raa].f.f,TMPFILE);
if(DEBUG>radebugflag){
fprintf(DBG,"System call string is %s.\n",rtmp);
fflush(DBG);
}
	system(rtmp);
	SYS=fopen(TMPFILE,"r");
	if(SYS==NULL){ 
printf("Error opening atomic temporary file #%d.  Exiting.\n",raa);
		exit(1);
		}
	fscanf(SYS,"%d",&raz);
if(DEBUG>radebugflag){
fprintf(DBG,"raz is %d.\n",raz);
fflush(DBG);
}
	AT[raa].n=raz;
	fclose(SYS);
	AT[raa].x=(double*)calloc(raz,sizeof(double));
	AT[raa].A=(double*)calloc(raz,sizeof(double));
	AT[raa].pop=(double*)calloc(raz,sizeof(double));
	ATF=fopen(AT[raa].f.f,"r");
	if(ATF==NULL){
printf("Error opening atomic info file #%d.  Exiting.\n",raa);
		exit(1);
		} 
	for(rab=0;rab<raz;rab++){
		fscanf(ATF,"%lf %lf %lf",&AT[raa].x[rab],&AT[raa].A[rab],\
				&AT[raa].pop[rab]); 
if(DEBUG>radebugflag){
fprintf(DBG,"xposition is %f, A is %f\n",AT[raa].x[rab],AT[raa].A[rab]); 
fflush(DBG);
}
		}
	fclose(ATF);
	} 
return; 
}


/**************** read_rot_dist_file() *****************/

/* This function reads rotational distribution files. 

   A note on the J(K) array:  positions in the array are given by:

   rdc*rdJK*rdv  +  rdd*rdv  +  rde

   where rdc is the current Omega number (not the value of this Omega);
	rdJK is the total number of lines in the J(K) distribution file
	as determined by "wc -l filename";  rdv is the number of vibration
	levels to be considered; rdd is the current J(K) number; rde is
	the current vibrational level number.  */ 

void read_rot_dist_files(){

int rdnumf=1,rda=0,rdb=0,rdc=0,rdd=0,rde=0,rdv=0,rdJK=0,rdsz=0,rdcs=0;
int rdcc=0,rddd=0,rdee=0,rdnum=0;
double rdJmax=0,rdJprev=0,rdJdum=0;
char rdsys[500];

sprintf(TMPFILE,".%s.temporary",PREF); 
/* start loop over molecules and states.  Reinitialize rdnumf */
for(rda=0;rda<NUMMOL;rda++){
if(DEBUG>rddebugflag){
fprintf(DBG,"Top of rda loop for rda=%d.  NUMMOL is %d\n",rda,NUMMOL);
fflush(DBG);
}
	for(rdb=0;rdb<MOL[rda].states;rdb++){ 
if(DEBUG>rddebugflag){
fprintf(DBG,"Top of rdb loop for rdb=%d\n",rdb);
fflush(DBG);
}
/* For each state, determine if there are files to read and, if so, 
   is the state Case a or Case b, and the position in the Case array */
if(MOL[rda].s[rdb].Dist==-1){
/* determine the number of vibrational levels to consider */	
	rdv=MOL[rda].s[rdb].v[0].vnum;
if(DEBUG>rddebugflag){
fprintf(DBG,"State %s (molecule %s) has distribution file(s) ",\
		MOL[rda].s[rdb].Name,MOL[rda].Mol);
fprintf(DBG,"and vnum is %d\n",MOL[rda].s[rdb].v[0].vnum);
fflush(DBG);
}
/* if Case a, figure out what rdnumf really is */
	if(MOL[rda].s[rdb].Ca>-1){
		rdcs=MOL[rda].s[rdb].Ca;
		rdnumf=CA[rdcs].O;
		rdJK=0;
		for(rdc=0;rdc<rdnumf;rdc++){
			sprintf(rdsys,"wc -l %s > %s",\
				MOL[rda].s[rdb].f[rdc].f,TMPFILE);
			system(rdsys);
if(DEBUG>rddebugflag){
fprintf(DBG,"System string is %s\n",rdsys);
fflush(DBG);
}
			SYS=fopen(TMPFILE,"r");
			if(SYS==NULL){
printf("Error opening temporary file for file %s. Exiting.\n",\
			MOL[rda].s[rdc].f[rdc].f);
				exit(1);
				}
			fscanf(SYS,"%d",&rdd);
			if(rdd>rdJK) rdJK=rdd;
			fclose(SYS);
			} 
		rdsz=rdnumf*rdv*rdJK;
		CA[rdcs].Jpop=(double*)calloc(rdsz,sizeof(double));
if(DEBUG>rddebugflag){
fprintf(DBG,"MOL[%d].Mol (%s) MOL[%d].s[%d].Name (%s) ",rda,MOL[rda].Mol,\
		rda,rdb,MOL[rda].s[rdb].Name);
fprintf(DBG,"is defined as Case a.\nThere are %d lines in the file.\n",rdJK);
fflush(DBG);
}
		}
	else{
		if(MOL[rda].s[rdb].Cb==-1){
printf("Molecule Name %s State %s has no defined case.  Exiting.\n",\
			MOL[rda].Mol,MOL[rda].s[rdb].Name);
			exit(1);
			}
		rdcs=MOL[rda].s[rdb].Cb;
		rdnum=1;
		sprintf(rdsys,"wc -l %s > %s",\
				MOL[rda].s[rdb].f[rdc].f,TMPFILE);
		system(rdsys);
		SYS=fopen(TMPFILE,"r");
		if(SYS==NULL){
printf("Error opening temporary file for file %s. Exiting.\n",\
		MOL[rda].s[rdc].f[rdc].f);
			exit(1);
			}
		fscanf(SYS,"%d",&rdJK);
		fclose(SYS);
		rdsz=rdv*rdJK;
		CB[rdcs].Jnum=rdJK;
		CB[rdcs].Jpop=(double*)calloc(rdsz,sizeof(double));
		CB[rdcs].K=(double*)calloc(rdJK,sizeof(double));
		CB[rdcs].J=(double*)calloc(rdJK,sizeof(double));
if(DEBUG>rddebugflag){
fprintf(DBG,"MOL[%d].Mol (%s) MOL[%d].s[%d].Name (%s) ",rda,MOL[rda].Mol,\
		rda,rdb,MOL[rda].s[rdb].Name);
fprintf(DBG,"is defined as Case b.\nThere are %d lines in the file.\n",rdJK);
fflush(DBG);
}
		}

/* start loop over files [rdc] -- there will only be one file for Case b or
   for Case a where a list of Omegas aren't specified */ 
	for(rdc=0;rdc<rdnumf;rdc++){
		rdcc=rdc*rdJK*rdv;
		rdJdum=rdJprev=rdJmax=0;
if(DEBUG>rddebugflag){
fprintf(DBG,"rdcc is %d, rdc is %d\n",rdcc,rdc);
fflush(DBG);
}
	MOL[rda].s[rdb].f[rdc].F=fopen(MOL[rda].s[rdb].f[rdc].f,"r"); 
		if(MOL[rda].s[rdb].f[rdc].F==NULL){
printf("Error opening distribution file %s.  Exiting.\n",\
		MOL[rda].s[rdb].f[rdc].f);
			exit(1);
			}
/* start loop over lines in file [rdd] and read in J/K values */
		for(rdd=0;rdd<rdJK;rdd++){ 
			rddd=rdd*rdv;
			fscanf(MOL[rda].s[rdb].f[rdc].F,"%lf",&rdJdum);
if(DEBUG>rddebugflag){
fprintf(DBG,"rddd id %d, rdd is %d, rdJdum is %f\n",rddd,rdd,rdJdum);
fflush(DBG);
}
/* check to make sure that the J(K) values in the file are sequential */
			if(rdd==0){
				rdJprev=rdJdum;
				if(MOL[rda].s[rdb].Ca>-1){
					CA[rdcs].Jmin=rdJdum;
if(DEBUG>rddebugflag){
fprintf(DBG,"CA[%d].Jmin is %f\n",rdcs,CA[rdcs].Jmin);
fflush(DBG);
}
					}
				else{ 
					CB[rdcs].Kmin=rdJdum;
					CB[rdcs].K[rdd]=rdJdum;
	fscanf(MOL[rda].s[rdb].f[rdc].F,"%lf",&CB[rdcs].J[rdd]);
if(DEBUG>rddebugflag){
fprintf(DBG,"CB[%d].Kmin is %f\n",rdcs,CB[rdcs].Kmin);
fprintf(DBG,"CB[%d].K[%d] is %f\n",rdcs,rdd,CB[rdcs].K[rdd]);
fprintf(DBG,"CB[%d].J[%d] is %f\n",rdcs,rdd,CB[rdcs].J[rdd]);
fflush(DBG);
}
					}
				}
			else{
				if(rdJdum!=(rdJprev+1)){
printf("Non-sequential J(or K) values found in file %s ",\
		MOL[rda].s[rdb].f[rdc].f);
printf("between lines numbered %d and %d.\n",rdd,(rdd-1));
printf("The two J-values read are %f and %f.\n",rdJprev,rdJdum);
printf("Exiting.  Please edit distribution file and restart program.\n");
fprintf(PAR,"Non-sequential J(or K) values found in file %s ",\
		MOL[rda].s[rdb].f[rdc].f);
fprintf(PAR,"between lines numbered %d and %d.\n",rdd,(rdd-1));
fprintf(PAR,"The two J-values read are %f and %f.\n",rdJprev,rdJdum);
fprintf(PAR,"Exiting.  Please edit distribution file and restart program.\n");
				exit(1);
					} 
				if(MOL[rda].s[rdb].Cb>-1){
					CB[rdcs].K[rdd]=rdJdum;
	fscanf(MOL[rda].s[rdb].f[rdc].F,"%lf",&CB[rdcs].J[rdd]);
if(DEBUG>rddebugflag){
fprintf(DBG,"CB[%d].K[%d] is %f\n",rdcs,rdd,CB[rdcs].K[rdd]);
fprintf(DBG,"CB[%d].J[%d] is %f\n",rdcs,rdd,CB[rdcs].J[rdd]);
fflush(DBG);
}
					}
				rdJprev=rdJdum;
/* check for maximum J(K) value. */
				if(rdJdum>rdJmax) rdJmax=rdJdum;
if(DEBUG>rddebugflag){
fprintf(DBG,"rdJmax is %f\n",rdJmax);
fflush(DBG);
}
				}
/* read in the populations for each vib level [rde] */
			for(rde=0;rde<rdv;rde++){
				rdee=rdcc+rddd+rde; 
				if(MOL[rda].s[rdb].Ca>-1){
					fscanf(MOL[rda].s[rdb].f[rdc].F,\
						"%lf",&CA[rdcs].Jpop[rdee]);
if(DEBUG>rddebugflag){
fprintf(DBG,"rdee is %d, rde is %d, ",rdee,rde);
fprintf(DBG,"CA[%d].Jpop[%d] is %f\n",rdcs,rdee,CA[rdcs].Jpop[rdee]);
fflush(DBG);
}
					}
				if(MOL[rda].s[rdb].Cb>-1){
					fscanf(MOL[rda].s[rdb].f[rdc].F,\
						"%lf",&CB[rdcs].Jpop[rdee]);
if(DEBUG>rddebugflag){
fprintf(DBG,"rdee is %d, rde is %d, ",rdee,rde);
fprintf(DBG,"CB[%d].Jpop[%d] is %f\n",rdcs,rdee,CB[rdcs].Jpop[rdee]);
fflush(DBG);
}
					} 
				} /* close read over v-levels */
			} /* close loop over lines in file (J or K values) */
		fclose(MOL[rda].s[rdb].f[rdc].F);
		if(MOL[rda].s[rdb].Ca>-1){ 
			CA[rdcs].Jmax=rdJmax;
if(DEBUG>rddebugflag){
fprintf(DBG,"rdJmax is %f CA[%d].Jmax is %f\n",rdJmax,rdcs,CA[rdcs].Jmax);
fflush(DBG);
}			}
		if(MOL[rda].s[rdb].Cb>-1){
			CB[rdcs].Kmax=rdJmax;
if(DEBUG>rddebugflag){
fprintf(DBG,"rdJmax is %f CB[%d].Kmax is %f\n",rdJmax,rdcs,CB[rdcs].Kmax);
fflush(DBG);
}			
			} 
		} /* close loop over rdnumf (num Omegas if Case a, else 1) */ 
} /* close if Dist==-1 condition */
		} /* close loop over number of states */
	} /* close loop over number of molecules */ 
return;
}


/**************** read_state_and_transition_files() *****************/

/* This function reads the transition file and the state file.  */
void read_state_and_transition_files(){

char tdum[2000],tdumm[8],tsys[2000],tstring[2000],thom[15]="N"; 
int ta=0,tb=0,tc=0,td=0,te=0,tz=0,ty=0,tx=0,tw=0,tv=0,tu=0,tt=0;
int tisomega=1,tnumomega=1,tf=0;
char tmolck[201],tmoldum[201],tstdum[201],tsysdir[1000];
double tmult=0,tomeglo=0,tomeghi=0,tmodO=1,tmodS=1,tvmaxpop=0;
FILE *TSYSF;

INTR.F=fopen(INTR.f,"r");
if(INTR.F==NULL){
printf("Error opening transition file.  Exiting.\n");
exit(1);
} 
fscanf(INTR.F,"%lf",&TEMP); 
if(DEBUG>tdebugflag){
fprintf(DBG,"%f\n",TEMP);
fflush(DBG);
} 
if(TEMP==0){ /* write warning to parameter file */
fprintf(PAR,"WARNING!!!  Temperature NOT defined globally.\n");
fprintf(PAR,"If there is not a temperature designation for any molecule or\n");
fprintf(PAR,"state, that molecule or state will be at T=0 (brrr....).\n\n");
	}
else{ /* write global temperature to parameter file */
fprintf(PAR,"Temperature defined globally as %f.\n\n",TEMP);
	}
	fflush(PAR);

fscanf(INTR.F,"%d",&NUMMOL); 
if(DEBUG>tdebugflag){
fprintf(DBG,"NUMMOL is %d\n",NUMMOL);
fflush(DBG);
}

/******* BEGIN LOOP OVER MOLECULES **************/
MOL=(Molecule *)calloc(NUMMOL, sizeof(Molecule)); 
for(ta=0;ta<NUMMOL;ta++){ 
	sprintf(tmoldum,"Molecule number %d",ta);
	strcpy(tstdum,"No state specified yet");
	fscanf(INTR.F,"%s",tmolck);
	if(strcmp(tmolck,"MOL")!=0){
printf("MOL entry not found at beginning of data for molecule");
printf("number %d.  Exiting. \nPlease edit transition file and",(ta+1));
printf(" restart program.\n"); 
fprintf(PAR,"MOL entry not found at beginning of data for molecule");
fprintf(PAR,"number %d.  Exiting. \nPlease edit transition file and",(ta+1));
fprintf(PAR," restart program.\n"); 
		exit(1);
		}
	strcpy(tmolck,new_mol_ck(tmoldum,tstdum));
	sscanf(tmolck,"%s",MOL[ta].Mol);
	strcpy(tmolck,new_mol_ck(MOL[ta].Mol,tstdum));
	sscanf(tmolck,"%lf",&MOL[ta].pop); 
	strcpy(tmolck,new_mol_ck(MOL[ta].Mol,tstdum));
	sscanf(tmolck,"%lf",&MOL[ta].T); 
	sprintf(tsysdir,"mkdir %s_molecules/%s",PREF,MOL[ta].Mol);
	system(tsysdir);
if(DEBUG>tdebugflag){
fprintf(DBG,"The molecule is %s with pop %f and temp %f\n",MOL[ta].Mol,\
		MOL[ta].pop,MOL[ta].T); 
fflush(DBG);
}
	if(MOL[ta].T==0){
fprintf(PAR,"Temperature NOT specified for molecule %s.\n\n",MOL[ta].Mol);
		}
	else{
fprintf(PAR,"Temperature specified as %f for molecule %s.\n\n",\
		MOL[ta].T,MOL[ta].Mol);
		} 
	fflush(PAR);
	strcpy(tmolck,new_mol_ck(MOL[ta].Mol,tstdum));
	sscanf(tmolck,"%d",&MOL[ta].states); 
	NUMST+=MOL[ta].states; 
if(DEBUG>tdebugflag){
fprintf(DBG,"MOL[%d].states is %d\n",ta,MOL[ta].states);
fflush(DBG);
} 
	MOL[ta].s=(State *)calloc(MOL[ta].states, sizeof(State));
	for(tb=0;tb<MOL[ta].states;tb++){
		MOL[ta].s[tb].Ca=-1;
		MOL[ta].s[tb].Cb=-1;
		MOL[ta].s[tb].Cc=-1;
		MOL[ta].s[tb].Cd=-1;
		} 
if(DEBUG>tdebugflag){
fprintf(DBG,"Mol. %s Rel. Pop. is %f at Temp. %f\tMOL[%d].states is %d\n",\
		MOL[ta].Mol,MOL[ta].pop,MOL[ta].T,ta,MOL[ta].states);
fflush(DBG);
}

/*** OPEN LOOP OVER STATES ***/
	for(tb=0;tb<MOL[ta].states;tb++){
		tisomega=1;
		tnumomega=1;
		sprintf(tstdum,"State number %d",(tb+1));
		strcpy(tmolck,new_mol_ck(MOL[ta].Mol,tstdum));
		sscanf(tmolck,"%s",MOL[ta].s[tb].Name); 
		strcpy(tmolck,new_mol_ck(MOL[ta].Mol,MOL[ta].s[tb].Name));
		sscanf(tmolck,"%lf",&MOL[ta].s[tb].pop); 
/* read state info before getting to the Omega issue */
/*sprintf(tsys,"grep \"%s %s[ |\t]\" %s  > .%s.temporary",MOL[ta].Mol,\
   <<- this is an alternate form of the next half-line.  */
sprintf(tsys,"grep \"%s %s[[:space:]]\" %s  > .%s.temporary",MOL[ta].Mol,\
	MOL[ta].s[tb].Name,INST.f,PREF); 
		sprintf(tstring,".%s.temporary",PREF);
		system(tsys); 
if(DEBUG>tdebugflag){
fprintf(DBG,"tsys is %s; tstring is %s\n",tsys, tstring);
fflush(DBG);
} 
		TSYSF=fopen(tstring,"r");
		if(TSYSF==NULL){
printf("Error opening temporary file # %d for molecule %s in read_trans.  \
Exiting.\n",tb,MOL[ta].Mol);
			exit(1);
			}
		fscanf(TSYSF,"%s",tdumm);
		if(strcmp(tdumm,MOL[ta].Mol)!=0){
printf("Error reading MOL in temporary file for tdumm=%s and MOL[%d].Mol=%s \
and MOL[%d].s[%d].Name=%s.  Exiting.\n",tdumm,ta,MOL[ta].Mol,ta,tb,\
		MOL[ta].s[tb].Name);
			exit(1);
			}
		fscanf(TSYSF,"%s",tdumm);
		if(strcmp(tdumm,MOL[ta].s[tb].Name)!=0){
printf("Error reading state in temporary file for tdumm=%s and MOL[%d].Mol=%s\
 and MOL[%d].s[%d].Name=%s.  Exiting.\n",tdumm,ta,MOL[ta].Mol,ta,tb,\
		MOL[ta].s[tb].Name);
			exit(1);
			}
		fscanf(TSYSF,"%s",MOL[ta].s[tb].Case);
		if(strcmp(MOL[ta].s[tb].Case,"a")==0){
			tc=Case_a_num;
			MOL[ta].s[tb].Ca=Case_a_num; 
if(DEBUG>tdebugflag){
fprintf(DBG,"MOL[%d].s[%d].Ca is %d\n",ta,tb,MOL[ta].s[tb].Ca);
fflush(DBG);
}
			CA=(Case_a_stateinfo*)realloc(CA,\
					((tc+1)*sizeof(Case_a_stateinfo)));
			fscanf(TSYSF,"%d",&CA[tc].L); 
if(DEBUG>tdebugflag){
fprintf(DBG,"State %s of molecule %s has Lambda=%d and is \n",\
		MOL[ta].s[tb].Name,MOL[ta].Mol,CA[tc].L); 
fprintf(DBG,"defined as Case a.\n"); 
fflush(DBG);
} 
			if(CA[tc].L==0){
fprintf(PAR,"State %s of molecules %s ",MOL[ta].s[tb].Name,MOL[ta].Mol);
fprintf(PAR,"has Lambda=0\n and is defined as Case a.\n");
fprintf(PAR,"This program isn't set up to call Sigma states Case a.\n");
fprintf(PAR,"Redefining this state as Case b.\n");
fflush(PAR);
				strcpy(MOL[ta].s[tb].Case,"b");
				tc=Case_b_num;
				MOL[ta].s[tb].Cb=Case_b_num;
				MOL[ta].s[tb].Ca=-1;
				CB=(Case_b_stateinfo*)realloc(CB,\
					((tc+1)*(sizeof(Case_b_stateinfo))));
				CB[tc].L=0;
				fscanf(TSYSF,"%d",&CB[tc].p);
				fscanf(TSYSF,"%lf",&tmult);
				CB[tc].S=(tmult-1)/2;
				CB[tc].sflag=fmod(CB[tc].S,1);
				fscanf(TSYSF,"%s",thom);
				if(strcmp(thom,"Y")==0){
					fscanf(TSYSF,"%d",&CB[tc].g);
					fscanf(TSYSF,"%lf",&CB[tc].I);
					}
				fscanf(TSYSF,"%lf",&CB[tc].Te);
				fscanf(TSYSF,"%lf",&CB[tc].we);
				fscanf(TSYSF,"%lf",&CB[tc].wexe);
				fscanf(TSYSF,"%lf",&CB[tc].Be);
				fscanf(TSYSF,"%lf",&CB[tc].ae); 
				tnumomega=1;  /* Case b, tnumomega = 1 */
				Case_b_num++; 
if(DEBUG>tdebugflag){
fprintf(DBG,"Redefined as Case (b):  MOL[%d].s[%d].Cb is %d\n",\
		ta,tb,MOL[ta].s[tb].Cb);
fprintf(DBG,"tc: %d; Mult: %f; S: %f\n",tc,tmult,CB[tc].S); 
fprintf(DBG,"Homonuclear? %s; g: %d;  I: %f; \n",thom,CB[tc].g,CB[tc].I);
fprintf(DBG,"Te=%f; we=%f; wexe=%f; Be=%f; ae=%f\n",CB[tc].Te,CB[tc].we,\
		CB[tc].wexe,CB[tc].Be,CB[tc].ae);
fflush(DBG);
}
if(CB[tc].Be==0){
	printf("It isn't possible to simulate a rotational spectrum if\n");
	printf("the rotational constant (Be) is zero.  Call me picky. :-)\n"); 
	printf("The program will exit.  Either edit the state file for\n");
	printf("state %s of molecule %s or delete its transition(s).\n",\
		MOL[ta].s[tb].Name,MOL[ta].Mol); 	
	exit(1);
	}
				}
			else{
				fscanf(TSYSF,"%lf",&tmult); 
				CA[tc].S=(tmult-1)/2;
				CA[tc].sflag=fmod(CA[tc].S,1);
				fscanf(TSYSF,"%s",thom);
				if(strcmp(thom,"Y")==0){
					fscanf(TSYSF,"%d",&CA[tc].g);
					fscanf(TSYSF,"%lf",&CA[tc].I);
					}
				fscanf(TSYSF,"%lf",&CA[tc].Te);
				fscanf(TSYSF,"%lf",&CA[tc].we);
				fscanf(TSYSF,"%lf",&CA[tc].wexe);
				fscanf(TSYSF,"%lf",&CA[tc].Be);
				fscanf(TSYSF,"%lf",&CA[tc].ae); 
if(DEBUG>tdebugflag){
fprintf(DBG,"tc: %d; Mult: %f; S: %f\n",tc,tmult,CA[tc].S); 
fprintf(DBG,"Homonuclear? %s; g: %d;  I: %f; \n",thom,CA[tc].g,CA[tc].I);
fprintf(DBG,"Te=%f; we=%f; wexe=%f; Be=%f; ae=%f\n",CA[tc].Te,CA[tc].we,\
		CA[tc].wexe,CA[tc].Be,CA[tc].ae);
fflush(DBG);
}
				tnumomega=tmult;  /* Case a, tnumomega = the
multiplicity unless the user specifies otherwise (see below) */
				Case_a_num++; 
if(CA[tc].Be==0){
	printf("It isn't possible to simulate a rotational spectrum if\n");
	printf("the rotational constant (Be) is zero.  Call me picky. :-)\n"); 
	printf("The program will exit.  Either edit the state file for\n");
	printf("state %s of molecule %s or delete its transition(s).\n",\
		MOL[ta].s[tb].Name,MOL[ta].Mol); 	
	exit(1);
	}
				}
			} 
if((strcmp(MOL[ta].s[tb].Case,"b")==0)&&(MOL[ta].s[tb].Cb==-1)){ 
			tc=Case_b_num;
			MOL[ta].s[tb].Cb=Case_b_num;
CB=(Case_b_stateinfo*)realloc(CB,((tc+1)*sizeof(Case_b_stateinfo)));
			fscanf(TSYSF,"%d",&CB[tc].L);
			if(CB[tc].L==0){
				fscanf(TSYSF,"%d",&CB[tc].p);
				} 
if(DEBUG>tdebugflag){
fprintf(DBG,"State %s of molecules %s has Lambda=%d and is \n",\
		MOL[ta].s[tb].Name,MOL[ta].Mol,CB[tc].L); 
fprintf(DBG,"defined as Case b.\n"); 
fflush(DBG);
}
			fscanf(TSYSF,"%lf",&tmult);
			tmodS=(float)fmod(tmult,1); 
			if(tmodS!=0){
printf("Multiplicities must be integers.  For state %s ",MOL[ta].s[tb].Name);
printf("of molecule %s, ",MOL[ta].Mol);
printf("the multiplicity is %6.2f. Exiting.\n",tmult);
printf("Please edit transtion file or state file and restart program.\n");
fprintf(PAR,"Multiplicities must be integers.  ");
fprintf(PAR,"For state %s ",MOL[ta].s[tb].Name);
fprintf(PAR,"of molecule %s, ",MOL[ta].Mol);
fprintf(PAR,"the multiplicity is %6.2f. Exiting.\n",tmult);
fprintf(PAR,"Please edit transtion file or state file and restart program.\n");
				exit(1);
				}
			CB[tc].S=(tmult-1)/2;
			CB[tc].sflag=fmod(CB[tc].S,1);
			fscanf(TSYSF,"%s",thom);
			if(strcmp(thom,"Y")==0){
				fscanf(TSYSF,"%d",&CB[tc].g);
				fscanf(TSYSF,"%lf",&CB[tc].I);
				}
			fscanf(TSYSF,"%lf",&CB[tc].Te);
			fscanf(TSYSF,"%lf",&CB[tc].we);
			fscanf(TSYSF,"%lf",&CB[tc].wexe);
			fscanf(TSYSF,"%lf",&CB[tc].Be);
			fscanf(TSYSF,"%lf",&CB[tc].ae); 
if(DEBUG>tdebugflag){
fprintf(DBG,"tc: %d; Mult: %f; S: %f\n",tc,tmult,CB[tc].S); 
fprintf(DBG,"Homonuclear? %s; g: %d;  I: %f; \n",thom,CB[tc].g,CB[tc].I);
fprintf(DBG,"Te=%f; we=%f; wexe=%f; Be=%f; ae=%f\n",CB[tc].Te,CB[tc].we,\
		CB[tc].wexe,CB[tc].Be,CB[tc].ae);
fflush(DBG);
} 
			tnumomega=1;  /* for Case b, tnumomega is always 1 */
			Case_b_num++;
if(CB[tc].Be==0){
	printf("It isn't possible to simulate a rotational spectrum if\n");
	printf("the rotational constant (Be) is zero.  Call me picky. :-)\n"); 
	printf("The program will exit.  Either edit the state file for\n");
	printf("state %s of molecule %s or delete its transition(s).\n",\
		MOL[ta].s[tb].Name,MOL[ta].Mol); 	
	exit(1);
	}
			} 
		fclose(TSYSF);
		if((MOL[ta].s[tb].Cb==-1)&&(MOL[ta].s[tb].Ca==-1)){
printf("Case neither a nor b in temporary file for MOL[ta].s[tb].Case=%s and \
MOL[%d].Mol=%s and MOL[%d].s[%d].Name=%s.  Exiting.\n",MOL[ta].s[tb].Case,ta,
MOL[ta].Mol,ta,tb,MOL[ta].s[tb].Name);
			exit(1);
			} 

/* Read in statement about Omega Specifications.  Then, read in 
   vibrational level information and other info as appropriate */ 
		strcpy(tmolck,new_mol_ck(MOL[ta].Mol,MOL[ta].s[tb].Name));
		sscanf(tmolck,"%s",tdum);
		if(strcmp(tdum,"Y")==0){
			if(MOL[ta].s[tb].Ca==-1){
printf("Omega restriction requested for a state that isn't Case a\n");
printf("This program doesn't understand that. Exiting.\n");
printf("Please edit transtion (or state) file and restart program.\n\n");
fprintf(PAR,"Omega restriction requested for a state that isn't Case a\n");
fprintf(PAR,"This program doesn't understand that. Exiting.\n");
fprintf(PAR,"Please edit transtion (or state) file and restart program.\n\n"); 
fflush(PAR); 
				exit(1);
				} 
			tisomega=0; 
			tv=MOL[ta].s[tb].Ca;
			strcpy(tmolck,new_mol_ck(MOL[ta].Mol,\
				MOL[ta].s[tb].Name));
			sscanf(tmolck,"%d",&tw); 
				CA[tv].O=tw;
			if(tw<1){
				tw=2*CA[tv].S+1;
				tisomega=1;
fprintf(PAR,"WARNING!!! Omega restriction requested, but with the default\n");
fprintf(PAR,"!! number of omegas.  The program will not read for a set of\n");
fprintf(PAR,"!! user-specified omegas and populations.  If the program\n");
fprintf(PAR,"!! exits after complaining of a MOL entry in the wrong place,\n");
fprintf(PAR,"!! this might be the reason why.  \n\n");
				}
			tnumomega=tw; 
			CA[tv].lO=(double*)calloc(tw,sizeof(double));
			CA[tv].pO=(double*)calloc(tw,sizeof(double));
			tomeglo=fabs((double)CA[tv].L-CA[tv].S);
			tomeghi=(double)CA[tv].L+CA[tv].S;
			tmodS=(float)fmod(CA[tv].S,1);
			} 
MOL[ta].s[tb].nO=tnumomega;
MOL[ta].s[tb].f=(fileset *)calloc(tnumomega,sizeof(fileset));
MOL[ta].s[tb].T=(double *)calloc(tnumomega,sizeof(double));
MOL[ta].s[tb].v=(vset *)calloc(tnumomega,sizeof(vset)); 
if(DEBUG>tdebugflag){
fprintf(DBG,"tnumomega is %d.  Just allocated f, T, v.\n",tnumomega);
fflush(DBG);
	} 
	if(tisomega!=0){
		tnumomega=1;
		}
		for(te=0;te<tnumomega;te++){ 
			if(tisomega==0){
				strcpy(tmolck,new_mol_ck(MOL[ta].Mol,\
					MOL[ta].s[tb].Name));
					sscanf(tmolck,"%lf",&CA[tv].lO[te]); 
					tmodO=(float)fmod(CA[tv].lO[te],1);
	if((CA[tv].lO[te]<tomeglo)||(CA[tv].lO[te]>tomeghi)){
printf("For state %s of molecule %s, ",MOL[ta].s[tb].Name,MOL[ta].Mol);
printf("Omega values must lie between %f and %f.\n",tomeglo,tomeghi);
printf("You requested an Omega of %f.  Exiting.\n",CA[tv].lO[te]);
printf("Please edit the transition file (or the state file) and ");
printf("restart the program.\n\n"); 
fprintf(PAR,"For state %s of molecule %s, ",MOL[ta].s[tb].Name,MOL[ta].Mol);
fprintf(PAR,"Omega values must lie between %f and %f.\n",tomeglo,tomeghi);
fprintf(PAR,"You requested an Omega of %f.  Exiting.\n",CA[tv].lO[te]);
fprintf(PAR,"Please edit the transition file (or the state file) and ");
fprintf(PAR,"restart the program.\n\n"); 
		exit(1);
		}
	if(tmodS!=tmodO){
printf("For state %s of molecule %s, ",MOL[ta].s[tb].Name,MOL[ta].Mol);
printf("Omega values for spins equal to %3.1f ",CA[tv].S);
if(tmodS==0.5){printf("must be half-integral.\n");}
if(tmodS==0){printf("must be integral.\n");}
printf("You requested an Omega of %f.  Exiting.\n",CA[tv].lO[te]);
printf("Please edit the transition file (or the state file) and ");
printf("restart the program.\n\n"); 
fprintf(PAR,"For state %s of molecule %s, ",MOL[ta].s[tb].Name,MOL[ta].Mol);
fprintf(PAR,"Omega values for spins equal to %3.1f ",CA[tv].S);
if(tmodS==0.5){fprintf(PAR,"must be half-integral.\n");}
if(tmodS==0){fprintf(PAR,"must be integral.\n");} 
fprintf(PAR,"You requested an Omega of %f.  Exiting.\n",CA[tv].lO[te]);
fprintf(PAR,"Please edit the transition file (or the state file) and ");
fprintf(PAR,"restart the program.\n\n"); 
		exit(1);
		}
				strcpy(tmolck,new_mol_ck(MOL[ta].Mol,\
					MOL[ta].s[tb].Name));
				sscanf(tmolck,"%lf",&CA[tv].pO[te]);
				} 
                strcpy(tmolck,new_mol_ck(MOL[ta].Mol,MOL[ta].s[tb].Name)); 
		sscanf(tmolck,"%s",tdum); 
		tz=ty=tx=0;
		for(td=0;td<strlen(tdum);td++){
			ty=isdigit((int)tdum[td]); 
			if((int)tdum[td]==46) tx++;
			if((ty==0)&&((int)tdum[td]!=46)) tz=1;
			}
		if((tz==0)&&(tx<=1)){ 
			sscanf(tdum,"%lf",&MOL[ta].s[tb].T[te]);
			MOL[ta].s[tb].Dist=1;
fprintf(PAR,"Temp for state %s of Molecule %s is %f\n",MOL[ta].s[tb].Name,\
		MOL[ta].Mol,MOL[ta].s[tb].T[te]);
fflush(PAR); 
if(DEBUG>tdebugflag){
fprintf(DBG,"Temp for state %s of Molecule %s is %f\n",MOL[ta].s[tb].Name,\
		MOL[ta].Mol,MOL[ta].s[tb].T[te]);
fflush(DBG);
} 
			} /* close if tdum is a number (read temperature) */ 
		else{
			sscanf(tdum,"%s",MOL[ta].s[tb].f[te].f);
			MOL[ta].s[tb].Dist=-1;
			MOL[ta].s[tb].f[te].F=fopen(MOL[ta].s[tb].f[te].f,"r");
			if(MOL[ta].s[tb].f[te].F==NULL){
printf("Cannot open input file %s.  Exiting.\n",MOL[ta].s[tb].f[te].f);
				exit(1);
				} 
			fclose(MOL[ta].s[tb].f[te].F);
fprintf(PAR,"Distribution file for state %s of Molecule %s is %s\n",\
		MOL[ta].s[tb].Name,MOL[ta].Mol,MOL[ta].s[tb].f[te].f);
fflush(PAR); 
if(DEBUG>tdebugflag){
fprintf(DBG,"Dist (%d) file for state %s of Molecule %s is %s\n",\
	MOL[ta].s[tb].Dist,MOL[ta].s[tb].Name,\
	MOL[ta].Mol,MOL[ta].s[tb].f[te].f);
fflush(DBG);
} 
			} /* close if tdum isn't a number (open file) */ 
if(DEBUG>tdebugflag){
fprintf(DBG,"%s\t%f\t%s\t(%d)\n",MOL[ta].s[tb].Name,MOL[ta].s[tb].pop,tdum,\
		MOL[ta].s[tb].Dist); 
fflush(DBG);
} 
		strcpy(tmolck,new_mol_ck(MOL[ta].Mol,MOL[ta].s[tb].Name));
		sscanf(tmolck,"%d",&MOL[ta].s[tb].v[te].vnum);	
if(DEBUG>tdebugflag){
fprintf(DBG,"vnum is %d\n",MOL[ta].s[tb].v[te].vnum);
fflush(DBG);
} 
		if(MOL[ta].s[tb].v[te].vnum>0){
			strcpy(tmolck,new_mol_ck(MOL[ta].Mol,\
						MOL[ta].s[tb].Name));
			sscanf(tmolck,"%d",&MOL[ta].s[tb].v[te].vlo);	
if(DEBUG>tdebugflag){
fprintf(DBG,"vlo = %d  ",MOL[ta].s[tb].v[te].vlo);
fflush(DBG);
} 
MOL[ta].s[tb].v[te].p=(double*)calloc(MOL[ta].s[tb].v[te].vnum,\
		sizeof(double));
MOL[ta].s[tb].v[te].c=(int*)calloc(MOL[ta].s[tb].v[te].vnum,sizeof(int));
			} 
		tt=MOL[ta].s[tb].v[te].vlo;
		tvmaxpop=0;
		for(td=0;td<MOL[ta].s[tb].v[te].vnum;td++){
			strcpy(tmolck,new_mol_ck(MOL[ta].Mol,\
						MOL[ta].s[tb].Name));
			sscanf(tmolck,"%lf",&MOL[ta].s[tb].v[te].p[td]); 
			tvmaxpop+=MOL[ta].s[tb].v[te].p[td]; 
if(DEBUG>tdebugflag){
fprintf(DBG,"MOL[%d].s[%d].v[%d].p[%d]=%f\n",\
		ta,tb,te,td,MOL[ta].s[tb].v[te].p[td]);
fflush(DBG);
} 
			} 
/* Make sure the vibrational populations sum to one in the simulation */
		for(td=0;td<MOL[ta].s[tb].v[te].vnum;td++){
			MOL[ta].s[tb].v[te].p[td]/=tvmaxpop; 
			}
/* Beginning of scan for file containing rotational distributions */ 
		if(MOL[ta].s[tb].Dist==-1){ 
			if((MOL[ta].s[tb].Ca!=-1)&&(MOL[ta].s[tb].Cb!=-1)){
printf("MOL[%d].s[%d].Ca is %d\n",ta,tb,MOL[ta].s[tb].Ca);
printf("MOL[%d].s[%d].Cb is %d\n",ta,tb,MOL[ta].s[tb].Cb);
				printf("fix Ca/Cb assignment problem\n");
				exit(1);
				} 
			sprintf(tsys,"head -1 %s > .%s.temporary", \
					MOL[ta].s[tb].f[te].f,PREF);
			system(tsys);
			TSYSF=fopen(tstring,"r");
			if(TSYSF==NULL){
				printf("fix tstring 1\n");
				exit(1);
				}
			fscanf(TSYSF,"%lf",&tmult);
			if(MOL[ta].s[tb].Ca!=-1){
				CA[tc].Jmin=tmult;
				}
			if(MOL[ta].s[tb].Cb!=-1){
				CB[tc].Kmin=tmult;
				}
			fclose(TSYSF);
if(DEBUG>tdebugflag){
fprintf(DBG,"string for MOL[%d].s[%d].f[%d].f is %s\n",\
		ta,tb,te,MOL[ta].s[tb].f[te].f); 
if(MOL[ta].s[tb].Ca!=-1){
fprintf(DBG,"Jmin is %f, ",CA[tc].Jmin);
}
if(MOL[ta].s[tb].Cb!=-1){
fprintf(DBG,"Kmin is %f, ",CB[tc].Kmin);
}
fflush(DBG);
} 
			sprintf(tsys,"tail -1 %s > .%s.temporary", \
					MOL[ta].s[tb].f[te].f,PREF);
			system(tsys);
			TSYSF=fopen(tstring,"r");
			if(TSYSF==NULL){
				printf("fix tstring 2\n");
				exit(1);
				}
			fscanf(TSYSF,"%lf",&tmult);
			if(MOL[ta].s[tb].Ca!=-1){
				CA[tc].Jmax=tmult;
				}
			if(MOL[ta].s[tb].Cb!=-1){
				CB[tc].Kmax=tmult;
				}
			fclose(TSYSF);
if(DEBUG>tdebugflag){
fprintf(DBG,"string for MOL[%d].s[%d].f[%d].f is %s\n",\
		ta,tb,te,MOL[ta].s[tb].f[te].f); 
if(MOL[ta].s[tb].Ca!=-1){ fprintf(DBG,"Jmax is %f, ",CA[tc].Jmax); }
if(MOL[ta].s[tb].Cb!=-1){ fprintf(DBG,"Kmax is %f, ",CB[tc].Kmax); }
fflush(DBG);
} 
			sprintf(tsys,"wc -l %s > .%s.temporary", \
					MOL[ta].s[tb].f[te].f,PREF);
			system(tsys);
			TSYSF=fopen(tstring,"r");
			if(TSYSF==NULL){
				printf("fix tstring 3\n");
				exit(1);
				}
			fscanf(TSYSF,"%lf",&tmult);
			if(MOL[ta].s[tb].Ca!=-1){
				if(tmult!=(CA[tc].Jmax-CA[tc].Jmin+1)){
					printf("Jmax-Jmin range not wc -l\n");
					printf("This might cause trouble\n");
					}
				}
			if(MOL[ta].s[tb].Cb!=-1){
				if(tmult!=(CB[tc].Kmax-CB[tc].Kmin+1)){
					printf("Kmax-Kmin range not wc -l\n");
					printf("This might cause trouble\n");
					}
				}
			fclose(TSYSF); 
if(DEBUG>tdebugflag){
fprintf(DBG,"tmult is %f\n",tmult); 
fflush(DBG);
} 
			}  /* close read for rotation distribution file */ 
			} /* close loop over number of specified omegas --
			     or over the one set of vibrational states */ 
		} /* close states for loop */ 

/* Preliminary stuff for loop over transitions for this molecule */ 
	sprintf(tstdum,"(not a state, top of transition list)");
	strcpy(tmolck,new_mol_ck(MOL[ta].Mol,tstdum));
	sscanf(tmolck,"%d",&MOL[ta].trans); 
if(DEBUG>tdebugflag){
fprintf(DBG,"MOL[%d].trans is %d\n",ta,MOL[ta].trans);
fflush(DBG);
}
	MOL[ta].t=(Trans*)calloc(MOL[ta].trans,sizeof(Trans));
/* Start loop over transitions for this molecule */ 
	for(tb=0;tb<MOL[ta].trans;tb++){
		sprintf(tstdum,"(not a state, transition number %d)",(tb+1));
		strcpy(tmolck,new_mol_ck(MOL[ta].Mol,tstdum));
		sscanf(tmolck,"%s",tdum); 
/* Separate the two states in the transition */ 
		tw=tv=0;
		while(tw==0){
			MOL[ta].t[tb].Nhi[tv]=tdum[tv];
			if(tdum[tv]=='-'){
				MOL[ta].t[tb].Nhi[tv]='\0';
				tw=1;
				}
			tv++;
			}
		tu=tv;
		tw=0;
		while(tw==0){ 
			MOL[ta].t[tb].Nlo[tv-tu]=tdum[tv];
			if(tdum[tv]=='\0') tw=1;
			tv++;
			}
/* Find the position in the state array for each state in the transition */
		tw=tv=tu=tt=0;
		while((tw==0)&&(tv<MOL[ta].states)){ 
			if(strcmp(MOL[ta].s[tv].Name,MOL[ta].t[tb].Nhi)==0){
				MOL[ta].t[tb].Hi=tv;
				MOL[ta].s[tv].no+=1;
MOL[ta].s[tv].o=(int*)realloc(MOL[ta].s[tv].o,MOL[ta].s[tv].no*sizeof(int));
				MOL[ta].s[tv].o[MOL[ta].s[tv].no-1]=tb;
				tu=1;
if(DEBUG>tdebugflag){
fprintf(DBG,"MOL[%d].t[%d].Nhi %s matches ",ta,tb,MOL[ta].t[tb].Nhi);
fprintf(DBG,"MOL[%d].s[%d].Name %s and ",ta,tv,MOL[ta].s[tv].Name);
fprintf(DBG,"MOL[%d].t[%d].Hi is %d.\n",ta,tb,MOL[ta].t[tb].Hi);
fprintf(DBG,"MOL[%d].s[%d].no is %d.\n",ta,tv,MOL[ta].s[tv].no);
for(tf=0;tf<MOL[ta].s[tv].no;tf++){
fprintf(DBG,"State %s is an origination state in ",MOL[ta].s[tv].Name);
fprintf(DBG,"transition number %d.\n\n",tb); 
	}
fflush(DBG);
}
				}
			if(strcmp(MOL[ta].s[tv].Name,MOL[ta].t[tb].Nlo)==0){
				MOL[ta].t[tb].Lo=tv;
				MOL[ta].s[tv].nd+=1;
MOL[ta].s[tv].d=(int*)realloc(MOL[ta].s[tv].d,MOL[ta].s[tv].nd*sizeof(int));
				MOL[ta].s[tv].d[MOL[ta].s[tv].nd-1]=tb;
				tt=1;
if(DEBUG>tdebugflag){
fprintf(DBG,"MOL[%d].t[%d].Nlo %s matches ",ta,tb,MOL[ta].t[tb].Nlo);
fprintf(DBG,"MOL[%d].s[%d].Name %s and ",ta,tv,MOL[ta].s[tv].Name);
fprintf(DBG,"MOL[%d].t[%d].Lo is %d.\n",ta,tb,MOL[ta].t[tb].Lo);
fprintf(DBG,"MOL[%d].s[%d].nd is %d.\n",ta,tv,MOL[ta].s[tv].nd);
for(tf=0;tf<MOL[ta].s[tv].nd;tf++){
fprintf(DBG,"State %s is an destination state in ",MOL[ta].s[tv].Name);
fprintf(DBG,"transition number %d.\n\n",tb); 
	}
fflush(DBG);
}
				}
			if((tu==1)&&(tt==1)){
				tw=1;
				}
			tv++;
			} 
if(DEBUG>tdebugflag){
fprintf(DBG,"tdum is %s\n",tdum);
fflush(DBG);
}	
		if(tw==0){
printf("Missing definition for at least one state in transition %s.",tdum);
printf("  Exiting.\nPlease edit transition file and restart program\n");
fprintf(PAR,"Missing definition for at least one state in transition %s.",tdum);
fprintf(PAR,"  Exiting.\nPlease edit transition file and restart program.\n");
			exit(1);
			} 
		if(MOL[ta].s[MOL[ta].t[tb].Hi].Ca!=-1){
			tw=CA[MOL[ta].s[MOL[ta].t[tb].Hi].Ca].O;
			MOL[ta].t[tb].Ohi=tw;
			if(tw==0){tw=1;}
			}
		else{
			tw=1;
			}
		if(MOL[ta].s[MOL[ta].t[tb].Lo].Ca!=-1){
			tv=CA[MOL[ta].s[MOL[ta].t[tb].Lo].Ca].O;
			MOL[ta].t[tb].Olo=tv;
			if(tv==0){tv=1;}
			}
		else{
			tv=1;
			}
		MOL[ta].t[tb].f=(fileset*)calloc((tw*tv),sizeof(fileset));
		MOL[ta].t[tb].P=(double*)calloc((tw*tv),sizeof(double));
		MOL[ta].t[tb].cP=(int*)calloc((tw*tv),sizeof(int)); 
		for(te=0;te<tw;te++){
		for(td=0;td<tv;td++){
			strcpy(tmolck,new_mol_ck(MOL[ta].Mol,tstdum));
			sscanf(tmolck,"%lf",&MOL[ta].t[tb].P[te*tv+td]);
			strcpy(tmolck,new_mol_ck(MOL[ta].Mol,tstdum));
			sscanf(tmolck,"%s",MOL[ta].t[tb].f[te*tv+td].f);
if(DEBUG>tdebugflag){
fprintf(DBG,"Prob is %f; file is %s\n",MOL[ta].t[tb].P[te*tv+td],\
		MOL[ta].t[tb].f[te*tv+td].f);
fflush(DBG);
}	
if((MOL[ta].t[tb].P[te*tv+td]==0)&&(MOL[ta].t[tb].f[te*tv+td].f!="NULL")){
fprintf(PAR,"WARNING!!!  FCF-FILE specified for zero transition ");
fprintf(PAR,"probability for transtion %s-%s \n\tand specified Omegas ",\
		MOL[ta].t[tb].Nhi,MOL[ta].t[tb].Nlo);
fprintf(PAR,"numbers %d (Hi) and %d (Lo).\n",(te+1),(td+1));
fprintf(PAR,"\t(This is a non-fatal error, the program will continue)\n");
	}
if((MOL[ta].t[tb].P[te*tv+td]!=0)&&(MOL[ta].t[tb].f[te*tv+td].f=="NULL")){
printf("FCF-FILE specified as NULL for non-zero transition ");
printf("probability\nfor transtion %s-%s and specified Omegas ",\
		MOL[ta].t[tb].Nhi,MOL[ta].t[tb].Nlo);
printf("numbers %d (Hi) and %d (Lo).\n",(te+1),(td+1));
printf("Exiting.  Please edit transition file and restart program.\n");
fprintf(PAR,"FCF-FILE specified as NULL for non-zero transition ");
fprintf(PAR,"probability\nfor transtion %s-%s and specified Omegas ",\
		MOL[ta].t[tb].Nhi,MOL[ta].t[tb].Nlo);
fprintf(PAR,"numbers %d (Hi) and %d (Lo).\n",(te+1),(td+1));
fprintf(PAR,"Exiting.  Please edit transition file and restart program.\n");
	exit(1);
	} 
		if(MOL[ta].t[tb].P[te*tv+td]!=0){
			TSYSF=fopen(MOL[ta].t[tb].f[te*tv+td].f,"r");
if(DEBUG>tdebugflag){
fprintf(DBG,"FCF file for MOL[%d].t[%d].f[%d].f is %s.\n\n",\
		ta,tb,(te*tv+td),MOL[ta].t[tb].f[te*tv+td].f);
fflush(DBG);
	}
			if(TSYSF==NULL){
printf("Error opening FCF file MOL[%d].t[%d].f[%d].f=%s.  Exiting.\n\n",\
		ta,tb,(te*tv+td),MOL[ta].t[tb].f[te*tv+td].f);
				exit(1);
				}
			fclose(TSYSF); 
			} /* close if trans. prob. not zero file check */
			} /* close td loop over omegas-in-transitions */ 
			} /* close te loop over omegas-in-transitions */
		} /* close transitions for loop */ 
	} /* end NUMMOL loop */ 
fclose(INTR.F);
return;
} 


/**************** read_XY_file *****************/

/* This function can be used to read in any file whose contents consist
   solely of two columns of numbers. */
void read_XY_file(XYlist *XY){

int XYa=0;
char XYsys[500];
FILE *XYF;

sprintf(XYsys,"wc -l %s > %s",XY[0].f.f,TMPFILE);
system(XYsys); 
if(DEBUG>XYdebugflag){
fprintf(DBG,"XYsys is %s.\n",XYsys);
fflush(DBG);
} 
SYS=fopen(TMPFILE,"r");
if(SYS==NULL){
        printf("Error opening system temporary file.  Exiting.\n");
        exit(1);
        } 
fscanf(SYS,"%d",&XY[0].n);
fclose(SYS);

if(DEBUG>XYdebugflag){
fprintf(DBG,"Here 1\n");
fflush(DBG);
} 
XY[0].x=(double*)calloc(XY[0].n,sizeof(double));
if(DEBUG>XYdebugflag){
fprintf(DBG,"Here 2\n");
fflush(DBG);
}
XY[0].y=(double*)calloc(XY[0].n,sizeof(double));
if(DEBUG>XYdebugflag){
fprintf(DBG,"Here 3\n");
fflush(DBG);
} 
XYF=fopen(XY[0].f.f,"r");
if(XYF==NULL){
        printf("Error opening XY input file %s.  Exiting.\n",XY[0].f.f);
        exit(1);
        } 
if(DEBUG>XYdebugflag){
fprintf(DBG,"Here 4\n");
fflush(DBG);
} 
for(XYa=0;XYa<XY[0].n;XYa++){
	fscanf(XYF,"%lf",&XY[0].x[XYa]);
if(DEBUG>XYdebugflag){
fprintf(DBG,"X is %f  ",XY[0].x[XYa]);
fflush(DBG);
}
	fscanf(XYF,"%lf",&XY[0].y[XYa]);
if(DEBUG>XYdebugflag){
fprintf(DBG,"Y is %f\n",XY[0].y[XYa]);
fflush(DBG);
} 
	}
fclose(XYF);
return;
} 
