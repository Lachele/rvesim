#include "rvesim.h"

/***************** a-a *****************/
  
/* This function calculates transitions between Hund's Cases (a) and (a) */

void a_a(int aaa, int aab){
/* indexes and dummy variables: */
int aac=0,aacc=0,aad=0,aadd=0,aae=0,aaee=0,aaf=0,aaff=0,aag=0;
int aah=0,aaj=0,aafvl=0,aafvh=0,aafrnh=0,aafrcl=0,aafrch=0;
int aadum=0;
/* for variables below: n="native"; c="cascade"; St="state"; Ca="Case a";
   v="vibrational"; num="number"; des="designation"; hi="high"; lo="low";
   O="Omega"; max="maximum"; pec="pseudo-Einstein coefficient"; 
   comp="comparison"; a,b,<etc> are indexes */
/* variables about vibrational levels */
int aahvnlo=0,aahvnnum=0,aahvclo=0,aahvcnum=0,aahivlo=0,aahivnum=0,aaHV=0;
int aalvnlo=0,aalvnnum=0,aalvclo=0,aalvcnum=0,aalovlo=0,aalovnum=0,aaLV=0;
/* variables for info about the two states */
int aaStlo=0,aaSthi=0,aaStnlovlo=30,aaStnlovhi=0,aaStclovlo=30,aaStclovhi=0;
int aaCahi=0,aaCalo=0,aaJsh=0;
double aahipop=0,aatrprob=0,aavhpop=0,aaBvDv=0,aaDe=0,aavfcf=0,aavfcfr=0;
/* variables about Omegas and J-values */
int aaOhnum=0, aaOlnum=0,aaOvhnum=0,aaOvlnum=0,aamaxhiJ=0;
double aahiOlo=0,aaloOlo=0,aamaxOmJ=0,aaOmJ=0,aaOlocurr=0,aaOhicurr=0;
/* variables related to which v-v transitions need to be included */
double aapecnmax=0,aapecncompb=0,aapeccmax=0,aapecccompb=0; 
/* variables for the spectrscopic constants */
double aawexehi=0,aawexelo=0,aaTehi=0,aaTelo=0,aawehi=0,aawelo=0;
/* some local pointers to use as abbreviations for the longer addresses */
Case_a_stateinfo *aaHI,*aaLO;
Trans *aaT;
rotset *aaRh,*aaRl; 
vset *aaVh;
/* a few flags for various purposes */
int aaflg=0,aafrnhf=0,aafrchf=0;
/* variables for creating and writing to output files */
char aafile[1000];
FILE *AAFILE;
/* variables for calculating Holn-London factors */
int aaL=0,aaB=0,aadL=0;
double aaJ=0;

/* Check FCF's for lower-state vibrational levels.   Find highest and lowest
   lower-state levels that wil be populated (within the user-set lower 
   precision limit).  Check to see that there are 
   energies for all those levels.  If not, reallocate the state info to have
   room for the new levels.  Don't add "native" population -- only cascade.
   (the added population will be ignored if it's a destination state only) */ 
aaSthi=MOL[aaa].t[aab].Hi;
aaStlo=MOL[aaa].t[aab].Lo;
aaCahi=MOL[aaa].s[aaSthi].Ca;
aaCalo=MOL[aaa].s[aaStlo].Ca;
aaLO=&CA[aaCalo];
aaHI=&CA[aaCahi];
aahipop=MOL[aaa].s[aaSthi].pop; /* upper state relative population */
aahvnlo=MOL[aaa].s[aaSthi].v[0].vlo; /* user-specified low vib number */
aahvnnum=MOL[aaa].s[aaSthi].v[0].vnum; /* user-specified number of vib levels */
aahvclo=MOL[aaa].s[aaSthi].v[0].vclo; /* low vib level from cascade */
aahvcnum=MOL[aaa].s[aaSthi].v[0].vcnum; /* number of vib levels from cascade */
aalvnlo=MOL[aaa].s[aaStlo].v[0].vlo; /* user-specified low vib number */
aalvnnum=MOL[aaa].s[aaStlo].v[0].vnum; /* user-specified number of vib levels */
aalvclo=MOL[aaa].s[aaStlo].v[0].vclo; /* low vib level from cascade */
aalvcnum=MOL[aaa].s[aaStlo].v[0].vcnum; /* number of vib levels from cascade */ 
aaOhnum=CA[aaCahi].O;
aaOlnum=CA[aaCalo].O;
aaOvhnum=CA[aaCahi].O; /* these because the number of sets of v-v transitions */
aaOvlnum=CA[aaCalo].O; /* might not be the same as the number of O-O trans's */
aawexehi=aaHI[0].wexe;
aawexelo=aaLO[0].wexe;
aaTehi=aaHI[0].Te;
aaTelo=aaLO[0].Te;
aawehi=aaHI[0].we;
aawelo=aaLO[0].we;
if((aaHI[0].L-aaLO[0].L)==+1) aadL=+1;
if((aaHI[0].L==aaLO[0].L)) aadL=+0;
if((aaHI[0].L-aaLO[0].L)==-1) aadL=-1; 
aaL=aaHI[0].L;

if(aadebugflag<DEBUG){
fprintf(DBG,"aaSthi=%d, aaStlo=%d, aaCahi=%d, aaCalo=%d, aahipop=%f\n",\
		aaSthi,aaStlo,aaCahi,aaCalo,aahipop);
fprintf(DBG,"aahvnlo=%d, aahvnnum=%d, aahvclo=%d, aahvcnum=%d\n",\
		aahvnlo,aahvnnum,aahvclo,aahvcnum);
fprintf(DBG,"aalvnlo=%d, aalvnnum=%d, aalvclo=%d, aalvcnum=%d\n",\
		aalvnlo,aalvnnum,aalvclo,aalvcnum); 
fprintf(DBG,"aaOhnum=%d, aaOlnum=%d, aaOvhnum=%d, aaOvlnum=%d\n",\
		aaOhnum,aaOlnum,aaOvhnum,aaOvlnum);
fflush(DBG);
}
/* There needs to be either (a) all specified Omegas or (b) all possible 
   Omegas.  Find all possible if they weren't specified.  Also, if they 
   weren't specified, there will only be one set of v-v transitions */
if(aaOhnum==0){
	aaOvhnum=1;
	aaOhnum=(int)(2*CA[aaCahi].S+1.000000001);/* hail, hail, o great
			 Numera / may we not have a truncation error :-) */
	aahiOlo=(double)CA[aaCahi].L-CA[aaCahi].S;
if(aadebugflag<DEBUG){
fprintf(DBG,"aaOhnum is zero and aaOvhnum=%d, aaOhnum=%d, aahiOlo=%f\n",\
	       aaOvhnum,aaOhnum,aahiOlo);
fflush(DBG);
}
	if(aahiOlo<0){
		aahiOlo*=-1; /* make a positive number for clarity */
		aahiOlo-=0.0000001; /* avoid truncation mishaps */
		aahiOlo=ceil(aaloOlo); /* get the next largest integer */
		aaOhnum-=(int)aaloOlo; /* subtract that from 2S+1 */
		aahiOlo+=(double)CA[aaCahi].L-CA[aaCahi].S;
if(aadebugflag<DEBUG){
fprintf(DBG,"aahiOlo<0 and aaOhnum=%d, aahiOlo=%f\n",aaOhnum,aahiOlo);
fflush(DBG);
}
		} 
	}
if(aaOlnum==0){
	aaOvlnum=1;
	aaOlnum=(int)(2*CA[aaCalo].S+1.000000001);
	aaloOlo=(double)CA[aaCalo].L-CA[aaCalo].S;
if(aadebugflag<DEBUG){
fprintf(DBG,"aaOlnum is zero and aaOvlnum=%d, aaOlnum=%d, aaloOlo=%f\n",\
	       aaOvlnum,aaOlnum,aaloOlo);
fflush(DBG);
}
	if(aaloOlo<0){
		aaloOlo*=-1; /* make a positive number for clarity */
		aaloOlo-=0.0000001; /* avoid truncation mishaps */
		aaloOlo=ceil(aaloOlo); /* get the next largest integer */
		aaOlnum-=(int)aaloOlo; /* subtract that from 2S+1 */
		aaloOlo+=(double)CA[aaCalo].L-CA[aaCalo].S;
if(aadebugflag<DEBUG){
fprintf(DBG,"aaloOlo<0 and aaOlnum=%d, aaloOlo=%f\n",aaOlnum,aaloOlo);
fflush(DBG);
}
		} 
	} 

/* find total number of sub-transitions (per v-v and also Omega-Omega) for 
   this transition, then allocate memory for them */
if((aahvnlo<aahvclo)||(aahvcnum==0)) aahivlo=aahvnlo;
else aahivlo=aahvclo;
if((aahvnlo+aahvnnum)>(aahvclo+aahvcnum)){
	aahivnum=aahvnlo+aahvnnum-aahivlo;
	}
else{
	aahivnum=aahvclo+aahvcnum-aahivlo;
	}
MOL[aaa].t[aab].vnh=aahivnum;
MOL[aaa].t[aab].vhlo=aahivlo;
if(aadebugflag<DEBUG){
fprintf(DBG,"aahvnlo=%d, aahvclo=%d ",aahvnlo,aahvclo);
fprintf(DBG,"(aahvnlo+aahvnnum)=%d, (aahvclo+aahvcnum)=%d\n",\
		(aahvnlo+aahvnnum),(aahvclo+aahvcnum));
fprintf(DBG,"aahivnum=%d, aahivlo=%d\n",aahivnum,aahivlo); 
fprintf(DBG,"aahivnum is %d, aahivlo is %d\n",MOL[aaa].t[aab].vnh,\
		MOL[aaa].t[aab].vhlo);
fflush(DBG);
}

/* loop to find maximum PEC (=FCF times FREQ^DETECTTYPE) */
for(aac=0;aac<aaOvhnum;aac++){
for(aad=0;aad<aaOvlnum;aad++){ /* for user-specified populations */ 
	aaHV=aahvnnum+aahvnlo;
if(aadebugflag<DEBUG){
fprintf(DBG,"aac=%d, aad=%d, aaHV=%d\n",aac,aad,aaHV);
fflush(DBG);
}
	for(aae=aahvnlo;aae<aaHV;aae++){ /* find max pec for native v's */
		for(aaf=0;aaf<30;aaf++){ 
			if(MOL[aaa].t[aab].v[aaee].fcfn[aaf]>aapecnmax){
				aapecnmax=MOL[aaa].t[aab].v[aaee].fcfn[aaf];
				}
if(aadebugflag<(DEBUG-1)){
fprintf(DBG,"aapecnmax=%11.4e\n",aapecnmax);
fflush(DBG);
}
			}
		}
	aaHV=aahvcnum+aahvclo;
if(aadebugflag<(DEBUG-1)){
fprintf(DBG,"Out of find-max-pec for native v's.  aaHV=%d\n",aaHV);
fflush(DBG);
}
	for(aae=aahvclo;aae<aaHV;aae++){ /* find max pec for cascade v's */
		for(aaf=0;aaf<30;aaf++){ 
			if(MOL[aaa].t[aab].v[aaee].fcfc[aaf]>aapeccmax){
				aapeccmax=MOL[aaa].t[aab].v[aaee].fcfc[aaf];
				}
if(aadebugflag<(DEBUG-1)){
fprintf(DBG,"aapeccmax=%11.4e\n",aapeccmax);
fflush(DBG);
}
			}
		}
	}
	}
/* find minimum and maximum lower-state v's with PEC's that are above the 
   low-intensity cutoff (JKCUT) */ 
for(aac=0;aac<aaOvhnum;aac++){
for(aad=0;aad<aaOvlnum;aad++){ /* for user-specified populations */ 
	aaHV=aahvnnum+aahvnlo;
if(aadebugflag<DEBUG){
fprintf(DBG,"In cut-off loop:  aac=%d, aad=%d, aaHV=%d\n",aac,aad,aaHV);
fflush(DBG);
}
	for(aae=aahvnlo;aae<aaHV;aae++){ /* for native populations */
		aaee=aac*aaOvlnum*30+aad*30+aae;
		aaf=0;
if(aadebugflag<(DEBUG-1)){
fprintf(DBG,"aaee=%d\n",aaee);
fflush(DBG);
}
		while(aaflg==0){
			aapecncompb=MOL[aaa].t[aab].v[aaee].fcfn[aaf]/JKCUT;
			if(aapecncompb<aapecnmax) aadum=aaf;
			else aaflg=1;
			aaf++;
			}
if(aadebugflag<(DEBUG-1)){
fprintf(DBG,"aapecncompb=%11.4e, aaf=%d\n",aapecncompb,aaf);
fflush(DBG);
}
		if(aadum<aaStnlovlo) aaStnlovlo=aadum;
		aaf=29;
		aaflg=0;
if(aadebugflag<(DEBUG-1)){
fflush(DBG);
}
		while(aaflg==0){
			aapecncompb=MOL[aaa].t[aab].v[aaee].fcfn[aaf]/JKCUT;
			if(aapecncompb<aapecnmax) aadum=aaf;
			else aaflg=1;
			aaf--;
			}
if(aadebugflag<(DEBUG-1)){
fprintf(DBG,"aapecncompb=%11.4e, aaf=%d, aadum=%d\n",aapecncompb,aaf,aadum);
fflush(DBG);
}
		if(aadum>aaStnlovhi) aaStnlovhi=aadum;
		}
	aaHV=aahvcnum+aahvclo;
if(aadebugflag<(DEBUG-1)){
fprintf(DBG,"out of native population cutoff. aaHV=%d\n",aaHV);
fflush(DBG);
}
	for(aae=aahvclo;aae<aaHV;aae++){ /* for cascade populations */
		aaee=aac*aaOvlnum*30+aad*30+aae;
		aaf=0;
if(aadebugflag<(DEBUG-1)){
fflush(DBG);
}
		while(aaflg==0){
			aapecccompb=MOL[aaa].t[aab].v[aaee].fcfc[aaf]/JKCUT;
			if(aapecccompb<aapeccmax) aadum=aaf; 
			else aaflg=1;
			aaf++;
			}
if(aadebugflag<(DEBUG-1)){
fprintf(DBG,"aapecccompb=%11.4e, aaf=%d, aadum=%d\n",aapecccompb,aaf,aadum);
fflush(DBG);
}
		if(aadum<aaStclovlo) aaStclovlo=aadum;
		aaf=29;
		aaflg=0;
if(aadebugflag<(DEBUG-1)){
fflush(DBG);
}
		while(aaflg==0){
			aapecccompb=MOL[aaa].t[aab].v[aaee].fcfc[aaf]/JKCUT;
			if(aapecccompb<aapeccmax) aadum=aaf;
			else aaflg=1;
			aaf--;
			}
if(aadebugflag<(DEBUG-1)){
fprintf(DBG,"aapecccompb=%11.4e, aaf=%d, aadum=%d\n",aapecccompb,aaf,aadum);
fflush(DBG);
}
		if(aadum>aaStclovhi) aaStclovhi=aadum;
		}
	}
	}
if(aadebugflag<(DEBUG-1)){
fprintf(DBG,"out of cutoff loop\n\n");
fflush(DBG);
}

/* determine the lowest necessary destination-state vib level and the
   number of destination-state vibration levels */
if((aaStnlovlo)<(aaStclovlo)) aalovlo=aaStnlovlo;
else aalovlo=aaStclovlo;
if((aaStnlovhi)>(aaStclovhi)) aalovnum=aaStnlovhi+1-aalovlo;
else aalovnum=aaStclovhi+1-aalovlo;
MOL[aaa].t[aab].vnl=aalovnum;
MOL[aaa].t[aab].vllo=aalovlo;
/* allocate memory for the simsets in the transition structure */
MOL[aaa].t[aab].nS=aaOhnum*aaOlnum*aahivnum*aalovnum;
MOL[aaa].t[aab].S=(simset*)calloc(MOL[aaa].t[aab].nS,sizeof(simset)); 
MOL[aaa].t[aab].fS=(int*)calloc(MOL[aaa].t[aab].nS,sizeof(int)); 
if(aadebugflag<DEBUG){
fprintf(DBG,"MOL[%d].t[%d].vnl=%d, MOL[%d].t[%d].vllo=%d ",aaa,aab,\
		MOL[aaa].t[aab].vnl,aaa,aab,MOL[aaa].t[aab].vllo);
fprintf(DBG,"MOL[%d].t[%d].nS=%d\n",aaa,aab,MOL[aaa].t[aab].nS);
fflush(DBG);
}

/* get maximum J value present in any of the high-state Omegas */
for(aae=0;aae<aaOhnum;aae++){
if(aadebugflag<DEBUG){
fprintf(DBG,"At top of J-value loop. aae=%d aaLO[0].O=%d\n",aae,aaLO[0].O);
fflush(DBG);
}
	aaHV=aahvnnum+aahvnlo;
	if(aaLO[0].O!=0){ /* if Omegas are user-defined */
		aaOlocurr=aaLO[0].lO[aae]; 
if(aadebugflag<DEBUG){
fprintf(DBG,"aaLO[0].O!=0; aaOlocurr=%f\n",aaLO[0].lO[aae]); 
fflush(DBG);
}
		}
	else{ 
		aaOlocurr=aaloOlo + aae; 
if(aadebugflag<DEBUG){
fprintf(DBG,"else, and aaOlocurr=%f\n",aaOlocurr); 
fflush(DBG);
}
		}
	for(aaf=aahvnlo;aaf<aaHV;aaf++){ /* check native populations */ 
if(aadebugflag<(DEBUG)){
fprintf(DBG,"In native population loop \n");
fflush(DBG);
}
		aaj=MOL[aaa].s[aaSthi].r[aaf-aahvnlo].j;
		aaOmJ=aaOlocurr+aaj-1;
		if(aaj>aamaxhiJ){ 
			aamaxhiJ=aaj;
			}
		if(aaOmJ>aamaxOmJ){ 
			aamaxOmJ=aaOmJ;
			}
		}
if(aadebugflag<(DEBUG)){
fprintf(DBG,"Native:  aae=%d, aaf=%d, aamaxhiJ=%d\n",aae,aaf,aamaxhiJ);
fflush(DBG);
}
	aaHV=aahvcnum+aahvclo;
	for(aaf=aahvclo;aaf<aaHV;aaf++){ /* check cascade populations */
		aaj=MOL[aaa].s[aaSthi].rc[aaf-aahvclo].j;
		aaOmJ=aaOlocurr+aaj-1;
		if(aaj>aamaxhiJ){ 
			aamaxhiJ=aaj;
			}
		if(aaOmJ>aamaxOmJ){ 
			aamaxOmJ=aaOmJ;
			}
		}
if(aadebugflag<(DEBUG)){
fprintf(DBG,"Cascade:  aae=%d, aaf=%d, aamaxhiJ=%d\n",aae,aaf,aamaxhiJ);
fflush(DBG);
}
	}

/* loop through the list of low-state Omegas first.  See if all the needed
   destination-state v's are calculated.  Check also that enough J's are
   calculated.  Calculate any states that aren't already done.  The 
   following method isn't efficient.  It will calculate more lower-state 
   energy levels than are actually needed. But, this shouldn't be a problem 
   other than having a few unnecessary calculations happen. */ 
aamaxhiJ++; /* because of the Delta-J = +1 possibility */
aamaxOmJ++; /* '' */
if(aadebugflag<(DEBUG)){
fprintf(DBG,"Before new calculations loop aamaxhiJ=%d\n",aamaxhiJ);
fflush(DBG);
}
/* Calculate any newly needed energy levels. */
/* loop through each lower Omega */
for(aad=0;aad<aaOlnum;aad++){ /* loop through each lower-state Omega */
	if(aaOvlnum==1) aaee=0; /* Only one vset if auto-Omega */
	else aaee=aad; 
	if(aaLO[0].O!=0){ /* if Omegas are user-defined */
		aaOlocurr=aaLO[0].lO[aad]; 
if(aadebugflag<(DEBUG)){
fprintf(DBG,"Omegas User-defined and aaOlocurr=%15.8e.\n",aaOlocurr);
fflush(DBG);
}
		}
	else{
		aaOlocurr=aaloOlo + aad; 
if(aadebugflag<(DEBUG)){
fprintf(DBG,"Omegas Not user-defined and aaOlocurr=%15.8e.\n",aaOlocurr);
fprintf(DBG,"aaloOlo=%15.8e; aaOlnum=%d\n",aaloOlo,aaOlnum);
fflush(DBG);
}
		}
	aaLV=aalovlo+aalovnum;
if(aadebugflag<(DEBUG)){
fprintf(DBG,"aad=%d, aaee=%d, aaLV=%d\n",aad,aaee,aaLV);
fflush(DBG);
}
	for(aaf=aalovlo;aaf<aaLV;aaf++){/* check low-state v-levels */
		aaff=aad*30+aaf;/* Index for the array of rc rotational sets */ 
if(aadebugflag<(DEBUG)){
fprintf(DBG,"aaf=%d, aaff=%d, aadum=%d, ",aaf,aaff,aadum);
fflush(DBG);
}
		if(MOL[aaa].s[aaStlo].rc[aaff].jc==0){
/* this v has never been calculated, so call a function to calculate it. */ 
		if(aaLO[0].we!=0){
			aaDe=4*pow(aaLO[0].Be,3)/(aaLO[0].we*aaLO[0].we);
			aaBvDv=(aaLO[0].Be-aaLO[0].ae*(aaf+0.5))/\
				(2*(aaDe+aaLO[0].beta*(aaf+0.5)));
			MOL[aaa].s[aaStlo].rc[aaff].Jdissoc=aaLO[0].sflag+\
				floor((sqrt(1+4*aaBvDv)-1)/2);
			}
		else MOL[aaa].s[aaStlo].rc[aaff].Jdissoc=RAND_MAX;
			if(aamaxhiJ<MOL[aaa].s[aaStlo].rc[aaff].Jdissoc){
				case_a_cascade(&MOL[aaa].s[aaStlo].rc[aaff],\
					aaCalo,aaOlocurr,aaf,0,aamaxhiJ);
				}
			else{
				case_a_cascade(&MOL[aaa].s[aaStlo].rc[aaff],\
					aaCalo,aaOlocurr,aaf,0,\
			(int)(MOL[aaa].s[aaStlo].rc[aaff].Jdissoc+1));
				}
if(aadebugflag<(DEBUG)){
fprintf(DBG,"just called case_a_cascade, v never calc'd\n");
fprintf(DBG,"aaDe=%12.6e, aaBvDv=%12.6e, ",aaDe,aaBvDv);
fprintf(DBG,"MOL[%d].s[%d].rc[%d].Jdissoc=%12.6e\n",aaa,aaStlo,aaff,\
		MOL[aaa].s[aaStlo].rc[aaff].Jdissoc); 
fflush(DBG);
fprintf(DBG,"jwas zero; just called cascade \n");
fprintf(DBG,"MOL[%d].s[%d].rc[%d].jc=%d\n",aaa,aaStlo,aaff,\
		MOL[aaa].s[aaStlo].rc[aaff].jc);
fflush(DBG);
}
			} 
		else{
/* this v has been calculated.  See if enough J's are calculated */
			if(MOL[aaa].s[aaStlo].rc[aaff].jc<aamaxhiJ){
if(aadebugflag<(DEBUG)){
fprintf(DBG,"j<maxJ; about to call cascade \n");
fprintf(DBG,"MOL[%d].s[%d].rc[%d].jc=%d\n",aaa,aaStlo,aaff,\
		MOL[aaa].s[aaStlo].rc[aaff].jc);
fflush(DBG);
}
			if(aamaxhiJ<MOL[aaa].s[aaStlo].rc[aaff].Jdissoc){
				case_a_cascade(&MOL[aaa].s[aaStlo].rc[aaff],\
					aaCalo,aaOlocurr,aaf,\
				MOL[aaa].s[aaStlo].rc[aaff].jc,aamaxhiJ);
				}
			else{
				case_a_cascade(&MOL[aaa].s[aaStlo].rc[aaff],\
		aaCalo,aaOlocurr,aaf,MOL[aaa].s[aaStlo].rc[aaff].jc,\
			(int)(MOL[aaa].s[aaStlo].rc[aaff].Jdissoc+1));
				}
if(aadebugflag<(DEBUG)){
fprintf(DBG,"j was <maxJ; just called cascade \n");
fprintf(DBG,"MOL[%d].s[%d].rc[%d].jc=%d\n",aaa,aaStlo,aaff,\
		MOL[aaa].s[aaStlo].rc[aaff].jc);
fflush(DBG);
fprintf(DBG,"just called case_a_cascade for cascade only\n");
fflush(DBG);
}
				}
			}
		} 
/* Reset values of vclo and vcnum for the next transition */ 
		MOL[aaa].s[aaStlo].v[aaee].vclo=aalovlo;
		MOL[aaa].s[aaStlo].v[aaee].vcnum=aalovnum;
if(aadebugflag<(DEBUG)){
fprintf(DBG,"MOL[%d].s[%d].v[%d].vclo=%d ",aaa,aaStlo,aaee,\
		MOL[aaa].s[aaStlo].v[aaee].vclo);
fprintf(DBG,"MOL[%d].s[%d].v[%d].vcnum=%d\n",aaa,aaStlo,aaee,\
		MOL[aaa].s[aaStlo].v[aaee].vcnum);
fflush(DBG);
}
	} 
	
/* Check whether larger selection rules (Lambda, S, +/- in Sigma states) are
 * violated.  If they are, write a note to that effect to the DBG file. */ 
if(aadebugflag<(DEBUG-1)){
fprintf(DBG,"aaHI[0].L=%d, aaLO[0].L=%d\n",aaHI[0].L,aaLO[0].L);
fflush(DBG);
}
if((abs(aaHI[0].L-aaLO[0].L))>1){
	fprintf(PAR,"\nTRANSITION OMITTED!!! (see below)\n\n");
	fprintf(PAR,"Normal selection rule violated for Delta-");
	fprintf(PAR,"Lambda=0,+/-1 for\n\tTransition %s-->%s\n",\
			MOL[aaa].t[aab].Nhi,MOL[aaa].t[aab].Nlo);
	fprintf(PAR,"\tThis program doesn't understand that.\n\n");
if(aadebugflag<(DEBUG)){
fprintf(DBG,"Delta-L selection rule violated.\n");
fflush(DBG);
}
	return;
	}
if((fabs(aaHI[0].S-aaLO[0].S))>0.01){ /* <<< an overkill against round-off */
	fprintf(PAR,"WARNING!!! Normal selection rule violated for Delta-");
	fprintf(PAR,"S=0 for\n\tTransition %s-->%s\n",\
			MOL[aaa].t[aab].Nhi,MOL[aaa].t[aab].Nlo);
	fprintf(PAR,"\tThis is non-fatal; only a warning.\n");
if(aadebugflag<(DEBUG)){
fprintf(DBG,"Delta-S selection rule violated\n");
fflush(DBG);
}
	} 

/* Start with the largest Omega of the high state (in case someone writes
   in cascade between different Omegas).  Check for presence of lower-state 
   Omegas to which to transit.  If present,
   transit, if not move on to the next higher Omega. 
  
   If Delta-spinSigma is not zero (true if, for example, a pi-pi transition
   occurs and one Omega is 2 and the other is 1), flag this set of transitions
   to be left out of simulation and placed in separate file.  At the moment I
   don't recall whether thst separate file actually gets written. */ 
aaT=&MOL[aaa].t[aab]; /* an abbreviation */
if(aadebugflag<DEBUG){
fprintf(DBG,"aaT.Nhi=%s, aaT.Nlo=%s\n",aaT[0].Nhi,aaT[0].Nlo);
fprintf(DBG,"High: MOL[%d].s[%d].v[0].vclo=%d ",aaa,aaSthi,\
		MOL[aaa].s[aaSthi].v[0].vclo);
fprintf(DBG,"High: MOL[%d].s[%d].v[0].vcnum=%d\n",aaa,aaSthi,\
		MOL[aaa].s[aaSthi].v[0].vcnum);
fprintf(DBG,"High: MOL[%d].s[%d].v[0].vlo=%d ",aaa,aaSthi,\
		MOL[aaa].s[aaSthi].v[0].vlo);
fprintf(DBG,"High: MOL[%d].s[%d].v[0].vnum=%d\n",aaa,aaSthi,\
		MOL[aaa].s[aaSthi].v[0].vnum);
fprintf(DBG,"aahivnum=%d; aalovnum=%d; aaOhnum=%d, aaOlnum=%d\n",\
		aahivnum,aalovnum,aaOhnum,aaOlnum);
fflush(DBG);
}
for(aac=(aaOhnum-1);aac>=0;aac--){ /* down through high-state Omegas */
	if(aaHI[0].O!=0){ /* if Omegas are user-defined */
		aaOhicurr=aaHI[0].lO[aac]; 
		aafvh=aac*aahivnum; 
		aaVh=&MOL[aaa].s[aaSthi].v[aac];
		aatrprob=aaT[0].P[aac]; 
		}
	else{
		aaOhicurr=aahiOlo + aac; 
		aafvh=0;
		aaVh=&MOL[aaa].s[aaSthi].v[0];
		aatrprob=aaT[0].P[0];
		}
	aafrnh=aac*aahvnnum;/* indexes for high state rotsets */
	aacc=aac*aaOlnum*aahivnum*aalovnum; /* posn. in 4 dimensions */
if(aadebugflag<DEBUG){
fprintf(DBG,"aac=%d, aafvh=%d, aatrprob=%f ",\
		aac,aafvh,aatrprob);
fprintf(DBG,"aafrnh=%d, aafrch=%d, aacc=%d\n",aafrnh,aafrch,aacc); 
fprintf(DBG,"aaOhicurr=%11.4e\n",aaOhicurr);
fflush(DBG);
}
for(aad=(aaOlnum-1);aad>=0;aad--){  /* down through low-state Omegas */
	aafrnh+=aahivnum;
	if(aaLO[0].O!=0){ /* if Omegas are user-defined */
		aaOlocurr=aaLO[0].lO[aad]; 
		aafvl=aad*aalovnum;
		}
	else{
		aaOlocurr=aaloOlo + aad; 
		aafvl=0;
		}
	aaJsh=(int)(aaOhicurr-aaOlocurr+0.0000001); /* shift in J indexing 
				 due to difference in Omegas */
	if((fabs(aaOhicurr-aaOlocurr))>1.000000001) aaflg=1;
	else aaflg=0;
	aadd=aacc+aad*aahivnum*aalovnum; /* posn. in remaining 3 dimensions */ 
if(aadebugflag<DEBUG){
fprintf(DBG,"aad=%d, aafrcl=%d, aadd=%d\n",aad,aafrcl,aadd);
fprintf(DBG,"aaOlocurr=%11.3e\n",aaOlocurr);
fflush(DBG);
} 
	for(aae=(aahivlo+aahivnum-1);aae>=aahivlo;aae--){/*high-state viblev*/
		aaee=aadd+(aae-aahivlo)*aalovnum; /* posn. in remaining 
						     2 dimensions */
		aafrnh--;
/* these flags (aafrnhf and aafrchf) tell if this high vib state is natively
   populated, populated by cascade, or both.  They will be used later */ 
		if((aae>=aahvnlo)&&(aae<(aahvnlo+aahvnnum))){
			aafrnhf=0;
			aavhpop=aaVh[0].p[aae-aahvnlo]; 
			}
		else aafrnhf=-1; 
		if((aae>=aahvclo)&&(aae<(aahvclo+aahvcnum))){
			aafrchf=0;
			}
		else aafrchf=-1;
if(aadebugflag<DEBUG){
fprintf(DBG,"aafrnhf=%d, aafrnh=%d, aavhpop=%f, aafrchf=%d\n",\
		aafrnhf,aafrnh,aavhpop,aafrchf);
fprintf(DBG,"aae=%d, aaee=%d\n",aae,aaee);
fprintf(DBG,"aahivlo=%d, aahivnum=%d, aafrnh=%d, aafrch=%d\n",\
		aahivlo,aahivnum,aafrnh,aafrch);
fflush(DBG);
} 
	aafrcl=aad*30+aalovnum; /* indexes for low state rotsets */
	for(aaf=(aalovlo+aalovnum-1);aaf>=aalovlo;aaf--){/*low-state viblev*/ 
if(aadebugflag<DEBUG){
fprintf(DBG,"aaf=%d, aalovlo=%d, aalovnum=%d \n",aaf,aalovlo,aalovnum);
fflush(DBG);
}
/* Right now, all Omegas get the same set of FCF's. */
		aaff=aaee+aaf-aalovlo; /* posn. in final dimension */
		aafrcl--; /* position in low-state cascade rotset */
		aaT[0].fS[aaff]=aaflg; 
		aavfcf=aaT[0].v[aae].fcfn[aaf];
		if(aavfcf!=0){aavfcfr=aaT[0].v[aae].fcfc[aaf]/aavfcf;}
		if((aaOhicurr==0)&&(aaOlocurr==0)) aaT[0].S[aaff].n=2*aamaxhiJ;
		else aaT[0].S[aaff].n=3*aamaxhiJ;
		aaT[0].S[aaff].f=\
			(double*)calloc(aaT[0].S[aaff].n,sizeof(double));
		aaT[0].S[aaff].ni=\
			(double*)calloc(aaT[0].S[aaff].n,sizeof(double));
		aaT[0].S[aaff].ci=\
			(double*)calloc(aaT[0].S[aaff].n,sizeof(double));
		aaRl=&MOL[aaa].s[aaStlo].rc[aafrcl]; 
if(aadebugflag<DEBUG){
fprintf(DBG,"aaff=%d, aafrcl=%d, aaT[0].S[aaff].n=%d\n",\
		aaff,aafrcl,aaT[0].S[aaff].n);
fflush(DBG);
}
/* start with highest J value (same reason as before), and loop down looking 
   for lower J's to which to transit.  assign intensities.  If both states are 
   Omega=0, then assign zero intensity to transitions for J=0<->J=0.  */ 
		if(aafrnhf==0){/* if this v is natively populated.*/
sprintf(aafile,"%s_molecules/%s/%s-%.1f--%s-%.1f_v%d--v%d_NAT.dat",\
	PREF,MOL[aaa].Mol,aaT[0].Nhi,aaOhicurr,\
		aaT[0].Nlo,aaOlocurr,aae,aaf);
			AAFILE=fopen(aafile,"w");
			if(AAFILE==NULL){
printf("Error opening transition sub-file %s.  Exiting.\n",aafile);
				exit(1);
				}
fprintf(AAFILE,"# File created by %s.\n# File contains ",PROGRAM_NAME);
if((aaOhicurr==0)&&(aaOlocurr==0)) fprintf(AAFILE,"P and R");
else fprintf(AAFILE,"P, Q and R");
fprintf(AAFILE," branches from the simulation titled %s\n",PREF);
fprintf(AAFILE,"# THESE INTENSITIES ARE DUE TO USER-SPECIFIED POPULATIONS");
fprintf(AAFILE," ONLY -- NO CASCADE\n");
fprintf(AAFILE,"# Transition is Hund's Case (a) to Hund's Case (a).\n");
fprintf(AAFILE,"# Molecule %s, \tHigh state: %s, Omega=%.1f, v=%d\n",\
		MOL[aaa].Mol,aaT[0].Nhi,aaOhicurr,aae);
fprintf(AAFILE,"#\t\tLow state: %s, Omega=%.1f, v=%d\n# Columns are: ",\
		aaT[0].Nlo,aaOlocurr,aaf);
if((aaOhicurr==0)&&(aaOlocurr==0)){
fprintf(AAFILE,"J_hi P(J_hi+1)_v P_i  R(J_hi-1)_v R_i\n#\n");
	}
else fprintf(AAFILE,"J_hi P(Jhi+1)_v P_i Q(J_hi+0)_v Q_i R(J_hi-1)_v R_i\n#\n");

			aaRh=&MOL[aaa].s[aaSthi].r[aafrnh]; 
			if((aamaxhiJ)<MOL[aaa].s[aaSthi].r[aafrnh].j){
				aaj=aamaxhiJ-1;
				}
			else aaj=MOL[aaa].s[aaSthi].r[aafrnh].j-1; 
			if((aaRl[0].Jdissoc-aaOlocurr)<(aaj-1)){ 
fprintf(PAR,"\nWARNING!! State %s(%.1f) of molecule %s has significant ",\
		MOL[aaa].s[aaSthi].Name,aaOhicurr,MOL[aaa].Mol);
fprintf(PAR,"rotational population at level %.1f\n",aaRh[0].J[aaj-1]);
fprintf(PAR,"\tBut lower state %s(%.1f) dissociates at level %.1f\n",\
		MOL[aaa].s[aaStlo].Name,aaOlocurr,aaRl[0].Jdissoc);
fprintf(PAR,"It is likely that the bound upper state is transiting to (an)");
fprintf(PAR," unbound lower state(s).\n");
fprintf(PAR,"Look for continuum emissions trailing off to the low-energy ");
fprintf(PAR,"end of the v(hi)=%d to v(lo)=%d\nband in the real ",aae,aaf);
fprintf(PAR,"spectrum (not the simulated one)\nCurrently, this program will");
fprintf(PAR," not simulate bound-->unbound spectra.\n\n");
				aaj=(int)(aaRl[0].Jdissoc-aaOlocurr);
				}
if(aadebugflag<DEBUG){
fprintf(DBG,"aafrnh=%d, aaj=%d\n",aafrnh,aaj);
fprintf(DBG,"MOL[%d].s[%d].r[%d].j=%d\n",aaa,aaSthi,aafrnh,\
		MOL[aaa].s[aaSthi].r[aafrnh].j);
fflush(DBG);
}
			for(aag=(aaj-1);aag>=0;aag--){
				aah=aag+aaJsh; 
				aaJ=aaRh[0].J[aag];
if(aadebugflag<DEBUG){
fprintf(DBG,"aah=%d, aag=%d, aaJsh=%d, aaJ=%.1f \n",aah,aag,aaJsh,aaJ);
fflush(DBG);
}
				if((aaOhicurr==0)&&(aaOlocurr==0)){ 
				if(((aah+1)>=0)&&((aah+1)<=aaj)){
					aaB=+1;
					aaT[0].S[aaff].f[aag*2+1]=\
						aaRh[0].EJ[aag]-\
						aaRl[0].EJ[aah+1]; 
					aaT[0].S[aaff].ni[aag*2+1]=\
						aahipop*aatrprob*aavhpop*\
						aaRh[0].PJ[aag]*aavfcf*\
						HL(aaB,aaJ,aaL,aadL);
					aaRl[0].CJ[aah+1]+=aavfcfr*\
						aaT[0].S[aaff].ni[aag*2+1];
if(aadebugflag<DEBUG){
fprintf(DBG,"aaRh[0].J[%d]=%.1f  ",aag,aaRh[0].J[aag]);
fprintf(DBG,"aaRh[0].EJ[%d]=%13.6e  ",aag,aaRh[0].EJ[aag]);
fprintf(DBG,"aaRl[0].J[%d]=%.1f  ",aah+1,aaRl[0].J[aah+1]);
fprintf(DBG,"aaRl[0].EJ[%d]=%13.6e\n",aah+1,aaRl[0].EJ[aah+1]);
fprintf(DBG,"aaRh[0].PJ[%d]=%13.6e  ",aag,aaRh[0].PJ[aag]);
fprintf(DBG,"aaRl[0].CJ[%d]=%13.6e\n",aah+1,aaRl[0].CJ[aah+1]);
fprintf(DBG,"aahipop=%13.6e aatrprob=%13.6e aavfcf=%13.6e aavhpop=%13.6e\n",\
		aahipop,aatrprob,aavfcf,aavhpop);
fflush(DBG);
}
					} 
				if(((aah-1)>=0)&&((aah-1)<=aaj)){
					aaB=-1;
					aaT[0].S[aaff].f[aaf*2+0]=\
						aaRh[0].EJ[aag]-\
						aaRl[0].EJ[aah-1]; 
					aaT[0].S[aaff].ni[aag*2+0]=\
						aahipop*aatrprob*aavhpop*\
						aaRh[0].PJ[aag]*aavfcf*\
						HL(aaB,aaJ,aaL,aadL);
					aaRl[0].CJ[aah-1]+=aavfcfr*\
						aaT[0].S[aaff].ni[aag*2+0];
if(aadebugflag<DEBUG){
fprintf(DBG,"aaRh[0].J[%d]=%.1f  ",aag,aaRh[0].J[aag]);
fprintf(DBG,"aaRh[0].EJ[%d]=%13.6e  ",aag,aaRh[0].EJ[aag]);
fprintf(DBG,"aaRl[0].J[%d]=%.1f  ",aah+1,aaRl[0].J[aah-1]);
fprintf(DBG,"aaRl[0].EJ[%d]=%13.6e\n",aah-1,aaRl[0].EJ[aah-1]);
fprintf(DBG,"aaRh[0].PJ[%d]=%13.6e  ",aag,aaRh[0].PJ[aag]);
fprintf(DBG,"aaRl[0].CJ[%d]=%13.6e\n",aah-1,aaRl[0].CJ[aah-1]);
fprintf(DBG,"aahipop=%13.6e aatrprob=%13.6e aavfcf=%13.6e aavhpop=%13.6e\n",\
		aahipop,aatrprob,aavfcf,aavhpop);
fflush(DBG);
}
					} 
if(aadebugflag<DEBUG){
fprintf(DBG,"aag=%d, aaT[0].S[aaff].f[aag*2+1]=%13.6e, ",aag,\
		aaT[0].S[aaff].f[aag*2+1]); 
fprintf(DBG,"aaT[0].S[aaff].f[aag*2+0]=%13.6e\n",aaT[0].S[aaff].f[aag*2+0]);
fprintf(DBG,"aaT[0].S[aaff].ni[aag*2+1]=%13.6e, ",aaT[0].S[aaff].ni[aag*2+1]);
fprintf(DBG,"aaRl[0].CJ[aah+1]=%13.6e\naaT[0].S[aaff].ni[aag*2+0]=%13.6e, ",\
		aaRl[0].CJ[aah+1],aaT[0].S[aaff].ni[aag*2+0]);
fprintf(DBG,"aaRl[0].CJ[aah-1]=%13.6e\n",aaRl[0].CJ[aah-1]);
fflush(DBG);
}
fprintf(AAFILE,"%5.1f\t%18.10e\t%18.10e\t%18.10e\t%18.10e\n",aaRh[0].J[aag],\
	aaT[0].S[aaff].f[aag*2+1],aaT[0].S[aaff].ni[aag*2+1],\
	aaT[0].S[aaff].f[aag*2+0],aaT[0].S[aaff].ni[aag*2+0]);
					}
				else{
				if(((aah+1)>=0)&&((aah+1)<=aaj)){ 
					aaB=+1;
					aaT[0].S[aaff].f[aag*3+2]=\
						aaRh[0].EJ[aag]-\
						aaRl[0].EJ[aah+1]; 
					aaT[0].S[aaff].ni[aag*3+2]=\
						aahipop*aatrprob*aavhpop*\
						aaRh[0].PJ[aag]*aavfcf*\
						HL(aaB,aaJ,aaL,aadL);
					aaRl[0].CJ[aah+1]+=aavfcfr*\
						aaT[0].S[aaff].ni[aag*3+2];
if(aadebugflag<DEBUG){
fprintf(DBG,"aaRh[0].J[%d]=%.1f  ",aag,aaRh[0].J[aag]);
fprintf(DBG,"aaRh[0].EJ[%d]=%13.6e  ",aag,aaRh[0].EJ[aag]);
fprintf(DBG,"aaRl[0].J[%d]=%.1f  ",aah+1,aaRl[0].J[aah+1]);
fprintf(DBG,"aaRl[0].EJ[%d]=%13.6e\n",aah+1,aaRl[0].EJ[aah+1]);
fprintf(DBG,"aaRh[0].PJ[%d]=%13.6e  ",aag,aaRh[0].PJ[aag]);
fprintf(DBG,"aaRl[0].CJ[%d]=%13.6e\n",aah+1,aaRl[0].CJ[aah+1]);
fprintf(DBG,"aahipop=%13.6e aatrprob=%13.6e aavfcf=%13.6e aavhpop=%13.6e\n",\
		aahipop,aatrprob,aavfcf,aavhpop);
fflush(DBG);
}
					}
				if(((aah+0)>=0)&&((aah+0)<=aaj)){
					aaB=0;
					aaT[0].S[aaff].f[aag*3+1]=\
						aaRh[0].EJ[aag]-\
						aaRl[0].EJ[aah+0]; 
					aaT[0].S[aaff].ni[aag*3+1]=\
						aahipop*aatrprob*aavhpop*\
						aaRh[0].PJ[aag]*aavfcf*\
						HL(aaB,aaJ,aaL,aadL);
					aaRl[0].CJ[aah+0]+=aavfcfr*\
						aaT[0].S[aaff].ni[aag*3+1];
if(aadebugflag<DEBUG){
fprintf(DBG,"aaRh[0].J[%d]=%.1f  ",aag,aaRh[0].J[aag]);
fprintf(DBG,"aaRh[0].EJ[%d]=%13.6e  ",aag,aaRh[0].EJ[aag]);
fprintf(DBG,"aaRl[0].J[%d]=%.1f  ",aah,aaRl[0].J[aah+0]);
fprintf(DBG,"aaRl[0].EJ[%d]=%13.6e\n",aah,aaRl[0].EJ[aah+0]);
fprintf(DBG,"aaRh[0].PJ[%d]=%13.6e  ",aag,aaRh[0].PJ[aag]);
fprintf(DBG,"aaRl[0].CJ[%d]=%13.6e\n",aah,aaRl[0].CJ[aah+0]);
fprintf(DBG,"aahipop=%13.6e aatrprob=%13.6e aavfcf=%13.6e aavhpop=%13.6e\n",\
		aahipop,aatrprob,aavfcf,aavhpop);
fflush(DBG);
}
					}
				if(((aah-1)>=0)&&((aah-1)<=aaj)){
					aaB=-1;
					aaT[0].S[aaff].f[aag*3+0]=\
						aaRh[0].EJ[aag]-\
						aaRl[0].EJ[aah-1]; 
					aaT[0].S[aaff].ni[aag*3+0]=\
						aahipop*aatrprob*aavhpop*\
						aaRh[0].PJ[aag]*aavfcf*\
						HL(aaB,aaJ,aaL,aadL);
					aaRl[0].CJ[aah-1]+=aavfcfr*\
						aaT[0].S[aaff].ni[aag*3+0];
if(aadebugflag<DEBUG){
fprintf(DBG,"aaRh[0].J[%d]=%.1f  ",aag,aaRh[0].J[aag]);
fprintf(DBG,"aaRh[0].EJ[%d]=%13.6e  ",aag,aaRh[0].EJ[aag]);
fprintf(DBG,"aaRl[0].J[%d]=%.1f  ",aah-1,aaRl[0].J[aah-1]);
fprintf(DBG,"aaRl[0].EJ[%d]=%13.6e\n",aah-1,aaRl[0].EJ[aah-1]);
fprintf(DBG,"aaRh[0].PJ[%d]=%13.6e  ",aag,aaRh[0].PJ[aag]);
fprintf(DBG,"aaRl[0].CJ[%d]=%13.6e\n",aah-1,aaRl[0].CJ[aah-1]);
fprintf(DBG,"aahipop=%13.6e aatrprob=%13.6e aavfcf=%13.6e aavhpop=%13.6e\n",\
		aahipop,aatrprob,aavfcf,aavhpop);
fflush(DBG);
}
					} 
if(aadebugflag<DEBUG){
fprintf(DBG,"aaT[0].S[aaff].f[aag*3+2]=%13.6e ",aaT[0].S[aaff].f[aag*3+2]);
fprintf(DBG,"aaT[0].S[aaff].f[aag*3+1]=%13.6e ",aaT[0].S[aaff].f[aag*3+1]);
fprintf(DBG,"aaT[0].S[aaff].f[aag*3+0]=%13.6e\n",aaT[0].S[aaff].f[aag*3+0]);
fprintf(DBG,"aaT[0].S[aaff].ni[aag*3+2]=%13.6e ",aaT[0].S[aaff].ni[aag*3+2]);
fprintf(DBG,"aaT[0].S[aaff].ni[aag*3+1]=%13.6e ",aaT[0].S[aaff].ni[aag*3+1]);
fprintf(DBG,"aaT[0].S[aaff].ni[aag*3+0]=%13.6e\n",aaT[0].S[aaff].ni[aag*3+0]);
fprintf(DBG,"aaRl[0].CJ[aah+1]=%13.6e ",aaRl[0].CJ[aah+1]); 
fprintf(DBG,"aaRl[0].CJ[aah+0]=%13.6e ",aaRl[0].CJ[aah+0]); 
fprintf(DBG,"aaRl[0].CJ[aah-1] =%13.6e\n",aaRl[0].CJ[aah-1]);
fflush(DBG);
}
fprintf(AAFILE,\
	"%5.1f\t%18.10e\t%18.10e\t%18.10e\t%18.10e\t%18.10e\t%18.10e\n",\
	aaRh[0].J[aag],\
	aaT[0].S[aaff].f[aag*3+2],aaT[0].S[aaff].ni[aag*3+2],\
	aaT[0].S[aaff].f[aag*3+1],aaT[0].S[aaff].ni[aag*3+1],\
	aaT[0].S[aaff].f[aag*3+0],aaT[0].S[aaff].ni[aag*3+0]);
					}
				}
			fclose(AAFILE);
			}

		if(aafrchf==0){ /* if this v is cascade populated */ 
sprintf(aafile,"%s_molecules/%s/%s-%.1f--%s-%.1f_v%d--v%d_CAS.dat",\
	PREF,MOL[aaa].Mol,aaT[0].Nhi,aaOhicurr,\
		aaT[0].Nlo,aaOlocurr,aae,aaf);
			AAFILE=fopen(aafile,"w");
			if(AAFILE==NULL){
printf("Error opening transition sub-file %s.  Exiting.\n",aafile);
				exit(1);
				}
fprintf(AAFILE,"# File created by %s.\n# File contains ",PROGRAM_NAME);
if((aaOhicurr==0)&&(aaOlocurr==0)) fprintf(AAFILE,"P and R");
else fprintf(AAFILE,"P, Q and R");
fprintf(AAFILE," branches from the simulation titled %s\n",PREF);
fprintf(AAFILE,"# THESE INTENSITIES ARE DUE TO CASCADE ");
fprintf(AAFILE," ONLY -- NO USER-SPECIFIED POPULATIONS\n");
fprintf(AAFILE,"# Transition is Hund's Case (a) to Hund's Case (a).\n");
fprintf(AAFILE,"# Molecule %s, \tHigh state: %s, Omega=%.1f, v=%d\n",\
		MOL[aaa].Mol,aaT[0].Nhi,aaOhicurr,aae);
fprintf(AAFILE,"#\t\tLow state: %s, Omega=%.1f, v=%d\n# Columns are: ",\
		aaT[0].Nlo,aaOlocurr,aaf);
if((aaOhicurr==0)&&(aaOlocurr==0)){
fprintf(AAFILE,"J_hi P(J_hi+1)_v P_i  R(J_hi-1)_v R_i\n#\n");
	}
else fprintf(AAFILE,"J_hi P(Jhi+1)_v P_i Q(J_hi+0)_v Q_i R(J_hi-1)_v R_i\n#\n");
			aaRh=&MOL[aaa].s[aaSthi].rc[aae]; 
			if((aamaxhiJ-1)<MOL[aaa].s[aaSthi].rc[aae].jc){
				aaj=aamaxhiJ-1;
				}
			else aaj=MOL[aaa].s[aaSthi].rc[aae].jc;
			if((aaRl[0].Jdissoc-aaOlocurr)<(aaj-1)){ 
fprintf(PAR,"\nWARNING!! State %s(%.1f) of molecule %s has significant ",\
		MOL[aaa].s[aaSthi].Name,aaOhicurr,MOL[aaa].Mol);
fprintf(PAR,"CASCADE rotational population at level %.1f\n",aaRh[0].J[aaj-1]);
fprintf(PAR,"\tBut lower state %s(%.1f) dissociates at level %.1f\n",\
		MOL[aaa].s[aaStlo].Name,aaOlocurr,aaRl[0].Jdissoc);
fprintf(PAR,"It is likely that the bound upper state is transiting to (an)");
fprintf(PAR," unbound lower state(s).\n");
fprintf(PAR,"Look for continuum emissions trailing off to the low-energy ");
fprintf(PAR,"end of the v(hi)=%d to v(lo)=%d\nband in the real ",aae,aaf);
fprintf(PAR,"spectrum (not the simulated one)\nCurrently, this program will");
fprintf(PAR," not simulate bound-->unbound spectra.\n\n");
				aaj=(int)(aaRl[0].Jdissoc-aaOlocurr);
				}
			for(aag=(aaj-1);aag>=0;aag--){
				aah=aag+aaJsh;
				aaJ=aaRh[0].J[aag];
				if((aaOhicurr==0)&&(aaOlocurr==0)){ 
				if(((aah+1)>=0)&&((aah+1)<=aaj)){
					aaB=+1;
					aaT[0].S[aaff].f[aag*2+1]=\
						aaRh[0].EJ[aag]-\
						aaRl[0].EJ[aah+1]; 
					aaT[0].S[aaff].ci[aag*2+1]=aatrprob*\
						aaRh[0].CJ[aag]*aavfcf*\
						HL(aaB,aaJ,aaL,aadL);
					aaRl[0].CJ[aah+1]+=aavfcfr*\
						aaT[0].S[aaff].ci[aag*2+1];
if(aadebugflag<DEBUG){
fprintf(DBG,"aaRh[0].EJ[%d]=%13.6e  ",aag,aaRh[0].EJ[aag]);
fprintf(DBG,"aaRl[0].EJ[%d]=%13.6e  ",aah+1,aaRl[0].EJ[aah+1]);
fprintf(DBG,"aaRh[0].CJ[%d]=%13.6e  ",aag,aaRh[0].CJ[aag]);
fprintf(DBG,"aaRl[0].CJ[%d]=%13.6e\n",aah+1,aaRl[0].CJ[aah+1]);
fprintf(DBG,"aavfcf=%13.6e  aatrprob=%13.6e \n",\
		aavfcf,aatrprob);
fflush(DBG);
}
					}
				if(((aah-1)>=0)&&((aah-1)<=aaj)){
					aaB=-1;
					aaT[0].S[aaff].f[aag*2+0]=\
						aaRh[0].EJ[aag]-\
						aaRl[0].EJ[aah-1]; 
					aaT[0].S[aaff].ci[aag*2+0]=aatrprob*\
						aaRh[0].CJ[aag]*aavfcf*\
						HL(aaB,aaJ,aaL,aadL);
					aaRl[0].CJ[aah-1]+=aavfcfr*\
						aaT[0].S[aaff].ci[aag*2+0];
if(aadebugflag<DEBUG){
fprintf(DBG,"aaRh[0].EJ[%d]=%13.6e  ",aag,aaRh[0].EJ[aag]);
fprintf(DBG,"aaRl[0].EJ[%d]=%13.6e  ",aah-1,aaRl[0].EJ[aah-1]);
fprintf(DBG,"aaRh[0].CJ[%d]=%13.6e  ",aag,aaRh[0].CJ[aag]);
fprintf(DBG,"aaRl[0].CJ[%d]=%13.6e\n",aah-1,aaRl[0].CJ[aah-1]);
fprintf(DBG,"aavfcf=%13.6e  aatrprob=%13.6e \n",\
		aavfcf,aatrprob);
fflush(DBG);
}
					} 
if(aadebugflag<DEBUG){
fprintf(DBG,"aag=%d, aaT[0].S[aaff].f[aag*2+1]=%13.6e, ",aag,\
		aaT[0].S[aaff].f[aag*2+1]); 
fprintf(DBG,"aaT[0].S[aaff].f[aag*2+0]=%13.6e\n",aaT[0].S[aaff].f[aag*2+0]);
fprintf(DBG,"aaT[0].S[aaff].ci[aag*2+1]=%13.6e, ",aaT[0].S[aaff].ci[aag*2+1]);
fprintf(DBG,"aaRl[0].CJ[aah+1]=%13.6e\naaT[0].S[aaff].ci[aag*2+0]=%13.6e, ",\
		aaRl[0].CJ[aah+1],aaT[0].S[aaff].ci[aag*2+0]);
fprintf(DBG,"aaRl[0].CJ[aah-1]=%13.6e\n",aaRl[0].CJ[aah-1]);
fflush(DBG);
}
fprintf(AAFILE,"%5.1f\t%18.10e\t%18.10e\t%18.10e\t%18.10e\n",\
	aaRh[0].J[aag],\
	aaT[0].S[aaff].f[aag*2+1],aaT[0].S[aaff].ci[aag*2+1],\
	aaT[0].S[aaff].f[aag*2+0],aaT[0].S[aaff].ci[aag*2+0]);
					}
				else{ 
				if(((aah+1)>=0)&&((aah+1)<=aaj)){
					aaB=+1;
					aaT[0].S[aaff].f[aag*3+2]=\
					aaRh[0].EJ[aag]-aaRl[0].EJ[aah+1]; 
					aaT[0].S[aaff].ci[aag*3+2]=aatrprob*\
						aaRh[0].CJ[aag]*aavfcf*\
						HL(aaB,aaJ,aaL,aadL);
					aaRl[0].CJ[aah+1]+=aavfcfr*\
						aaT[0].S[aaff].ci[aag*3+2];
if(aadebugflag<DEBUG){
fprintf(DBG,"aaRh[0].EJ[%d]=%13.6e  ",aag,aaRh[0].EJ[aag]);
fprintf(DBG,"aaRl[0].EJ[%d]=%13.6e  ",aah+1,aaRl[0].EJ[aah+1]);
fprintf(DBG,"aaRh[0].CJ[%d]=%13.6e  ",aag,aaRh[0].CJ[aag]);
fprintf(DBG,"aaRl[0].CJ[%d]=%13.6e\n",aah+1,aaRl[0].CJ[aah+1]);
fprintf(DBG,"aavfcf=%13.6e  aatrprob=%13.6e \n",\
		aavfcf,aatrprob);
fflush(DBG);
}
					}
				if(((aah+0)>=0)&&((aah+0)<=aaj)){ 
					aaB=0;
					aaT[0].S[aaff].f[aag*3+1]=\
					aaRh[0].EJ[aag]-aaRl[0].EJ[aah+0]; 
					aaT[0].S[aaff].ci[aag*3+1]=aatrprob\
						*aavfcf*aaRh[0].CJ[aag]*\
						HL(aaB,aaJ,aaL,aadL);
					aaRl[0].CJ[aah+0]+=aavfcfr*\
						aaT[0].S[aaff].ci[aag*3+1];
if(aadebugflag<DEBUG){
fprintf(DBG,"aaRh[0].EJ[%d]=%13.6e  ",aag,aaRh[0].EJ[aag]);
fprintf(DBG,"aaRl[0].EJ[%d]=%13.6e  ",aah,aaRl[0].EJ[aah+0]);
fprintf(DBG,"aaRh[0].CJ[%d]=%13.6e  ",aag,aaRh[0].CJ[aag]);
fprintf(DBG,"aaRl[0].CJ[%d]=%13.6e\n",aah,aaRl[0].CJ[aah+0]);
fprintf(DBG,"Here 20\n");
fprintf(DBG,"aavfcf=%13.6e  aatrprob=%13.6e\n",\
		aavfcf,aatrprob);
fflush(DBG);
}
					}
				if(((aah-1)>=0)&&((aah-1)<=aaj)){
					aaB=-1;
					aaT[0].S[aaff].f[aag*3+0]=\
					aaRh[0].EJ[aag]-aaRl[0].EJ[aah-1]; 
					aaT[0].S[aaff].ci[aag*3+0]=aatrprob\
						*aavfcf*aaRh[0].CJ[aag]*\
						HL(aaB,aaJ,aaL,aadL);
					aaRl[0].CJ[aah-1]+=aavfcfr*\
						aaT[0].S[aaff].ci[aag*3+0];
if(aadebugflag<DEBUG){
fprintf(DBG,"aaRh[0].EJ[%d]=%13.6e  ",aag,aaRh[0].EJ[aag]);
fprintf(DBG,"aaRl[0].EJ[%d]=%13.6e  ",aah-1,aaRl[0].EJ[aah-1]);
fprintf(DBG,"aaRh[0].CJ[%d]=%13.6e  ",aag,aaRh[0].CJ[aag]);
fprintf(DBG,"aaRl[0].CJ[%d]=%13.6e\n",aah-1,aaRl[0].CJ[aah-1]);
fprintf(DBG,"aavfcf=%13.6e  aatrprob=%13.6e\n",\
		aavfcf,aatrprob);
fflush(DBG);
}
					} 
if(aadebugflag<DEBUG){
fprintf(DBG,"aaT[0].S[aaff].f[aag*3+2]=%13.6e ",aaT[0].S[aaff].f[aag*3+2]);
fprintf(DBG,"aaT[0].S[aaff].f[aag*3+1]=%13.6e ",aaT[0].S[aaff].f[aag*3+1]);
fprintf(DBG,"aaT[0].S[aaff].f[aag*3+0]=%13.6e\n",aaT[0].S[aaff].f[aag*3+0]);
fprintf(DBG,"aaT[0].S[aaff].ci[aag*3+2]=%13.6e ",aaT[0].S[aaff].ci[aag*3+2]);
fprintf(DBG,"aaT[0].S[aaff].ci[aag*3+1]=%13.6e ",aaT[0].S[aaff].ci[aag*3+1]);
fprintf(DBG,"aaT[0].S[aaff].ci[aag*3+0]=%13.6e\n",aaT[0].S[aaff].ci[aag*3+0]);
fprintf(DBG,"aaRl[0].CJ[aah+1]=%13.6e ",aaRl[0].CJ[aah+1]); 
fprintf(DBG,"aaRl[0].CJ[aah+0]=%13.6e ",aaRl[0].CJ[aah+0]); 
fprintf(DBG,"aaRl[0].CJ[aah-1] =%13.6e\n",aaRl[0].CJ[aah-1]);
fflush(DBG);
}
fprintf(AAFILE,\
	"%5.1f\t%18.10e\t%18.10e\t%18.10e\t%18.10e\t%18.10e\t%18.10e\n",\
	aaRh[0].J[aag],\
	aaT[0].S[aaff].f[aag*3+2],aaT[0].S[aaff].ci[aag*3+2],\
	aaT[0].S[aaff].f[aag*3+1],aaT[0].S[aaff].ci[aag*3+1],\
	aaT[0].S[aaff].f[aag*3+0],aaT[0].S[aaff].ci[aag*3+0]);
					}
				}
			} 
		}
		}
	}
	} 
return; 
} 

