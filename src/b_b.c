#include "rvesim.h"

/***************** b-b *****************/
  
/* This function calls other functions that calculate transitions 
   between Hund's Case (b) and Hund's Case (b) */

void b_b(int bba, int bbb){
/* indexes and dummy variables: */
int bbc=0,bbd=0,bbe=0,bbf=0,bbff=0,bbk=0,bbdum=0;
/* for variables below: n="native"; c="cascade"; St="state"; Cb="Case b";
   v="vibrational"; num="number"; des="designation"; hi="high"; lo="low";
   max="maximum"; pec="pseudo-Einstein coefficient"; 
   comp="comparison"; a,b,<etc> are indexes */
/* variables about vibrational levels */
int bbhvnlo=0,bbhvnnum=0,bbhvclo=0,bbhvcnum=0,bbhivlo=0,bbhivnum=0,bbHV=0;
int bblvnlo=0,bblvnnum=0,bblvclo=0,bblvcnum=0,bblovlo=0,bblovnum=0,bbLV=0;
/* variables for info about the two states */
int bbStlo=0,bbSthi=0,bbStnlovlo=30,bbStnlovhi=0,bbStclovlo=30,bbStclovhi=0;
int bbCbhi=0,bbCblo=0;
double bbhipop=0,bbDe=0,bbBvDv=0;
/* variables related to which v-v transitions need to be included */
double bbpecnmax=0,bbpecncompb=0,bbpeccmax=0,bbpecccompb=0; 
/* variables to use for relating J and K */
int bbHiJnum=0,bbLoJnum=0,bbmaxhiK=0,bbtransperK=0;
double bbSh=0,bbSl=0;
/* variables for the spectrscopic constants */
double bbwexehi=0,bbwexelo=0,bbTehi=0,bbTelo=0,bbwehi=0,bbwelo=0;
/* some local pointers to use as abbreviations for the longer addresses */
Case_b_stateinfo *bbHI,*bbLO;
Trans *bbT;
/* a few flags for various purposes */
int bbflg=0;
/* a set of info for transferring to the sub-transition functions */
BBTinfo bbinfo;

/* Check FCF's for lower-state vibrational levels.   Find highest and lowest
   lower-state levels that wil be populated (within the user-set lower 
   precision limit).  Check to see that there are 
   energies for all those levels.  If not, reallocate the state info to have
   room for the new levels.  Don't add "native" population -- only cascade.
   (the added population will be ignored if it's a destination state only) */

/* a whole lot of the stuff below is redundant.  That's because I copied
   this function from a-a.c and it takes less time to be redundant than
   to make sure I changed everything correctly */
bbSthi=MOL[bba].t[bbb].Hi;
bbStlo=MOL[bba].t[bbb].Lo;
bbinfo.T=&MOL[bba].t[bbb];
bbCbhi=MOL[bba].s[bbSthi].Cb;
bbCblo=MOL[bba].s[bbStlo].Cb;
bbinfo.Sh=&MOL[bba].s[bbSthi];
bbinfo.Sl=&MOL[bba].s[bbStlo];
bbLO=&CB[bbCblo]; 
bbHI=&CB[bbCbhi];
bbinfo.Cl=bbLO;
bbinfo.Ch=bbHI; 
bbinfo.m=bba;
bbhipop=MOL[bba].s[bbSthi].pop; /* upper state relative population */
bbhvnlo=MOL[bba].s[bbSthi].v[0].vlo; /* user-specified low vib number */
bbhvnnum=MOL[bba].s[bbSthi].v[0].vnum; /* user-specified number of vib levels */
bbhvclo=MOL[bba].s[bbSthi].v[0].vclo; /* low vib level from cascade */
bbhvcnum=MOL[bba].s[bbSthi].v[0].vcnum; /* number of vib levels from cascade */
bblvnlo=MOL[bba].s[bbStlo].v[0].vlo; /* user-specified low vib number */
bblvnnum=MOL[bba].s[bbStlo].v[0].vnum; /* user-specified number of vib levels */
bblvclo=MOL[bba].s[bbStlo].v[0].vclo; /* low vib level from cascade */
bblvcnum=MOL[bba].s[bbStlo].v[0].vcnum; /* number of vib levels from cascade */ 
bbinfo.Vh=MOL[bba].s[bbSthi].v; 
bbinfo.hvnlo=bbhvnlo;
bbinfo.hvnnum=bbhvnnum;
bbinfo.hvclo=bbhvclo;
bbinfo.hvcnum=bbhvcnum;
bbinfo.lvnlo=bblvnlo;
bbinfo.lvnnum=bblvnnum;
bbinfo.lvclo=bblvclo;
bbinfo.lvcnum=bblvcnum; 
bbwexehi=bbHI[0].wexe;
bbwexelo=bbLO[0].wexe;
bbTehi=bbHI[0].Te;
bbTelo=bbLO[0].Te;
bbwehi=bbHI[0].we;
bbwelo=bbLO[0].we; 
if(bbdebugflag<DEBUG){
fprintf(DBG,"bbSthi=%d, bbStlo=%d, bbCbhi=%d, bbCblo=%d, bbhipop=%f\n",bbSthi,\
		bbStlo,bbCbhi,bbCblo,bbhipop);
fprintf(DBG,"bbhvnlo=%d, bbhvnnum=%d, bbhvclo=%d, bbhvcnum=%d\n",bbhvnlo,\
		bbhvnnum,bbhvclo,bbhvcnum);
fprintf(DBG,"bblvnlo=%d, bblvnnum=%d, bblvclo=%d, bblvcnum=%d\n",bblvnlo,\
		bblvnnum,bblvclo,bblvcnum); 
fprintf(DBG,"bbinfo.Sh[0].Cb=%d, bbinfo.Sl[0].Cb=%d\n",\
		bbinfo.Sh[0].Cb,bbinfo.Sl[0].Cb);
fprintf(DBG,"Molecule is %s;  Transition %s --> %s\n",MOL[bba].Mol,\
		MOL[bba].s[bbSthi].Name,MOL[bba].s[bbStlo].Name);
fprintf(DBG,"MOL[%d].t[%d].Nhi=%s, MOL[%d].t[%d].Nlo=%s\n",bba,bbb,\
		MOL[bba].t[bbb].Nhi,bba,bbb,MOL[bba].t[bbb].Nlo);
fflush(DBG);
}

/* Find the number of J's per (most) K(s) */
bbHiJnum=(int)(2*CB[bbCbhi].S+1.000000001);/* hail, hail, o great Numera 
	      / may we not have a truncation error :-) */
bbLoJnum=(int)(2*CB[bbCblo].S+1.000000001);
bbinfo.nhJ=bbHiJnum;
bbinfo.nlJ=bbLoJnum;
bbSh=CB[bbCbhi].S;
bbSl=CB[bbCblo].S;
if(bbdebugflag<DEBUG){
fprintf(DBG,"bbSh=%.1f, bbsL=%.1f\n",bbSh,bbSl);
fflush(DBG);
}
/* Find the maximum number of transitions to expect for each value of K.
   There must be some general form to use for this, but I spent 20 minutes
   working on it and decided I could just let the program calculate it
   by brute force before I could work out the general form.  This is 
   necessary to get the relative intensities right, but will be crucial if
   anyone ever adds any fine splitting effects.  See the documentation on
   transitions for more info. */
/* NOTE:::  This is a neat idea, but I'm running low on time and need
   to get this done.  So, I'll be being lazy a little bit further down. 
   Scan for "Laziness:". */
bbtransperK=0;
for(bbc=0;bbc<bbHiJnum;bbc++){
/* for K to K-1 transitions: */
	if(((-1-bbSl)<(-bbSh+bbc))&&((bbSl-1)>(-bbSh+bbc))) bbtransperK+=3;
	if(((-1-bbSl)==(-bbSh+bbc))&&((bbSl-1)>(-bbSh+bbc))) bbtransperK+=2;
	if(((-1-bbSl)<(-bbSh+bbc))&&((bbSl-1)==(-bbSh+bbc))) bbtransperK+=2;
	if((bbSl==0)&&((-bbSh+bbc)>=-2)&&((-bbSh+bbc)<=0)) bbtransperK+=1; 
if(bbdebugflag<DEBUG){
fprintf(DBG,"after Delta-K=-1 for bbc=%d, bbtransperK=%d\n",bbc,bbtransperK);
fflush(DBG);
} 
/* for K to K transitions: */
if((bbHI[0].L!=0)||(bbLO[0].L!=0)){ /* only if not Sigma-Sigma */
	if(((-bbSl)<(-bbSh+bbc))&&((bbSl)>(-bbSh+bbc))) bbtransperK+=3;
	if(((-bbSl)==(-bbSh+bbc))&&((bbSl)>(-bbSh+bbc))) bbtransperK+=2;
	if(((-bbSl)<(-bbSh+bbc))&&((bbSl)==(-bbSh+bbc))) bbtransperK+=2;
	if((bbSl==0)&&((-bbSh+bbc)>=-1)&&((-bbSh+bbc)<=1)) bbtransperK+=1; 
	}
if(bbdebugflag<DEBUG){
fprintf(DBG,"after Delta-K=0 for bbc=%d, bbtransperK=%d\n",bbc,bbtransperK);
fflush(DBG);
}
/* for K to K+1 transitions: */
	if(((+1-bbSl)<(-bbSh+bbc))&&((bbSl+1)>(-bbSh+bbc))) bbtransperK+=3;
	if(((+1-bbSl)==(-bbSh+bbc))&&((bbSl+1)>(-bbSh+bbc))) bbtransperK+=2;
	if(((+1-bbSl)<(-bbSh+bbc))&&((bbSl+1)==(-bbSh+bbc))) bbtransperK+=2;
	if((bbSl==0)&&((-bbSh+bbc)>=0)&&((-bbSh+bbc)<=2)) bbtransperK+=1; 
if(bbdebugflag<DEBUG){
fprintf(DBG,"after Delta-K=+1 for bbc=%d, bbtransperK=%d\n",bbc,bbtransperK);
fflush(DBG);
}
	}
bbinfo.nT=bbtransperK;
/* Laziness:  If the business above worked in all cases, the following 
   line wouldn't be necessary.  But, this way, there is definitely 
   enough space in the array.  In some cases, there will be far too much. */
bbinfo.nT=bbtransperK=(int)((2*bbSh+1)*3*(2*bbSl+1));
if(bbdebugflag<DEBUG){
fprintf(DBG,"The Lazy bbtransperK (bbinfo.nT) is %d\n",bbtransperK);
fflush(DBG);
}

/* find total number of sub-transitions (v-v) for this transition.  This
   will take some doing.  After, allocate memory for them */
if((bbhvnlo<bbhvclo)||(bbhvcnum==0)) bbhivlo=bbhvnlo;
else bbhivlo=bbhvclo;
if((bbhvnlo+bbhvnnum)>(bbhvclo+bbhvcnum)){
	bbhivnum=bbhvnlo+bbhvnnum-bbhivlo;
	}
else{
	bbhivnum=bbhvclo+bbhvcnum-bbhivlo;
	}
MOL[bba].t[bbb].vnh=bbhivnum;
MOL[bba].t[bbb].vhlo=bbhivlo;
bbinfo.hivnum=bbhivnum;
bbinfo.hivlo=bbhivlo;
if(bbdebugflag<DEBUG){
fprintf(DBG,"bbhvnlo=%d, bbhvclo=%d ",bbhvnlo,bbhvclo);
fprintf(DBG,"(bbhvnlo+bbhvnnum)=%d, (bbhvclo+bbhvcnum)=%d\n",\
		(bbhvnlo+bbhvnnum),(bbhvclo+bbhvcnum));
fprintf(DBG,"bbhivnum=%d, bbhivlo=%d\n",bbhivnum,bbhivlo); 
fprintf(DBG,"bbhivnum is %d, bbhivlo is %d\n",MOL[bba].t[bbb].vnh,\
		MOL[bba].t[bbb].vhlo);
fprintf(DBG,"MOL[%d].t[%d].Nhi=%s, MOL[%d].t[%d].Nlo=%s\n",bba,bbb,\
		MOL[bba].t[bbb].Nhi,bba,bbb,MOL[bba].t[bbb].Nlo);
fflush(DBG);
} 
/* loop to find maximum PEC (=FCF times FREQ^DETECTTYPE) */
bbHV=bbhvnnum+bbhvnlo;
if(bbdebugflag<DEBUG){
fprintf(DBG,"bbc=%d, bbd=%d, bbHV=%d\n",bbc,bbd,bbHV);
fflush(DBG);
}
for(bbe=bbhvnlo;bbe<bbHV;bbe++){ /* find max pec for native v's */
	for(bbf=0;bbf<30;bbf++){ 
		if(MOL[bba].t[bbb].v[bbe].fcfn[bbf]>bbpecnmax){
			bbpecnmax=MOL[bba].t[bbb].v[bbe].fcfn[bbf];
			}
if(bbdebugflag<(DEBUG-1)){
fprintf(DBG,"bbpecnmax=%11.4e\n",bbpecnmax);
fflush(DBG);
}
		}
	}
bbHV=bbhvcnum+bbhvclo;
if(bbdebugflag<(DEBUG-1)){
fprintf(DBG,"Out of find-max-pec for native v's.  bbHV=%d\n",bbHV);
fflush(DBG);
}
for(bbe=bbhvclo;bbe<bbHV;bbe++){ /* find max pec for cascade v's */
	for(bbf=0;bbf<30;bbf++){ 
		if(MOL[bba].t[bbb].v[bbe].fcfc[bbf]>bbpeccmax){
			bbpeccmax=MOL[bba].t[bbb].v[bbe].fcfc[bbf];
			}
if(bbdebugflag<(DEBUG-1)){
fprintf(DBG,"bbpeccmax=%11.4e\n",bbpeccmax);
fflush(DBG);
}
		}
	} 
/* find minimum and maximum lower-state v's with PEC's that are above the 
   low-intensity cutoff (JKCUT) */ 
bbHV=bbhvnnum+bbhvnlo;
if(bbdebugflag<DEBUG){
fprintf(DBG,"In cut-off loop:  bbc=%d, bbd=%d, bbHV=%d\n",bbc,bbd,bbHV); 
fprintf(DBG,"\tMOL[%d].t[%d].Nhi=%s, MOL[%d].t[%d].Nlo=%s\n",bba,bbb,\
		MOL[bba].t[bbb].Nhi,bba,bbb,MOL[bba].t[bbb].Nlo);
fflush(DBG);
}
for(bbe=bbhvnlo;bbe<bbHV;bbe++){ /* for native populations */
	bbf=0;
if(bbdebugflag<(DEBUG-1)){
fprintf(DBG,"\tMOL[%d].t[%d].Nhi=%s, MOL[%d].t[%d].Nlo=%s\n",bba,bbb,\
		MOL[bba].t[bbb].Nhi,bba,bbb,MOL[bba].t[bbb].Nlo);
fflush(DBG);
}
	while(bbflg==0){
		bbpecncompb=MOL[bba].t[bbb].v[bbe].fcfn[bbf]/JKCUT;
		if(bbpecncompb<bbpecnmax) bbdum=bbf;
		else bbflg=1;
		bbf++;
		}
	if(bbdum<bbStnlovlo) bbStnlovlo=bbdum;
if(bbdebugflag<(DEBUG-1)){
fprintf(DBG,"bbpecncompb=%11.4e, bbf=%d\n",bbpecncompb,bbf);
fprintf(DBG,"bbStnlovlo=%d\n",bbStnlovlo);
fflush(DBG);
} 
	bbf=29;
	bbflg=0;
if(bbdebugflag<(DEBUG-1)){
fflush(DBG);
}
	while(bbflg==0){
		bbpecncompb=MOL[bba].t[bbb].v[bbe].fcfn[bbf]/JKCUT;
		if(bbpecncompb<bbpecnmax) bbdum=bbf;
		else bbflg=1; 
if(bbdebugflag<(DEBUG-1)){
fprintf(DBG,"bbStnlovhi loop:  bbpecncompb=%11.4e , ",bbpecncompb); 
fprintf(DBG,"bbf=%d , bbdum=%d , bbflg=%d\n",bbf,bbdum,bbflg);
fflush(DBG);
} 
		bbf--;
		}
	if(bbdum>bbStnlovhi) bbStnlovhi=bbdum;
if(bbdebugflag<(DEBUG-1)){
fprintf(DBG,"bbpecncompb=%11.4e, bbf=%d, bbdum=%d\n",bbpecncompb,bbf,bbdum); 
fprintf(DBG,"bbStnlovhi=%d\n",bbStnlovhi);
fflush(DBG);
} 
	}
bbHV=bbhvcnum+bbhvclo;
if(bbdebugflag<(DEBUG-1)){
fprintf(DBG,"out of native population cutoff. bbHV=%d\n",bbHV);
fprintf(DBG,"\tMOL[%d].t[%d].Nhi=%s, MOL[%d].t[%d].Nlo=%s\n",bba,bbb,\
		MOL[bba].t[bbb].Nhi,bba,bbb,MOL[bba].t[bbb].Nlo);
fflush(DBG);
}
for(bbe=bbhvclo;bbe<bbHV;bbe++){ /* for cascade populations */
	bbf=0;
if(bbdebugflag<(DEBUG-1)){

fflush(DBG);
}
	while(bbflg==0){
		bbpecccompb=MOL[bba].t[bbb].v[bbe].fcfc[bbf]/JKCUT;
		if(bbpecccompb<bbpeccmax) bbdum=bbf;
		else bbflg=1;
		bbf++;
		}
if(bbdebugflag<(DEBUG-1)){
fprintf(DBG,"bbpecccompb=%11.4e, bbf=%d, bbdum=%d\n",bbpecccompb,bbf,bbdum);
fflush(DBG);
}
	if(bbdum<bbStclovlo) bbStclovlo=bbdum;
	bbf=29;
	bbflg=0;
	while(bbflg==0){
		bbpecccompb=MOL[bba].t[bbb].v[bbe].fcfc[bbf]/JKCUT;
		if(bbpecccompb<bbpeccmax) bbdum=bbf;
		else bbflg=1;
		bbf--;
		}
if(bbdebugflag<(DEBUG-1)){
fprintf(DBG,"bbpecccompb=%11.4e, bbf=%d, bbdum=%d\n",bbpecccompb,bbf,bbdum);
fflush(DBG);
}
	if(bbdum>bbStclovhi) bbStclovhi=bbdum;
	}
if(bbdebugflag<(DEBUG-1)){
fprintf(DBG,"out of cutoff loop\n\n");
fprintf(DBG,"MOL[%d].t[%d].Nhi=%s, MOL[%d].t[%d].Nlo=%s\n",bba,bbb,\
		MOL[bba].t[bbb].Nhi,bba,bbb,MOL[bba].t[bbb].Nlo);
fflush(DBG);
}

/* determine the lowest necessary destination-state vib level and the
   number of destination-state vibration levels */
if((bbStnlovlo)<(bbStclovlo)) bblovlo=bbStnlovlo;
else bblovlo=bbStclovlo;
if((bbStnlovhi)>(bbStclovhi)) bblovnum=bbStnlovhi+1-bblovlo;
else bblovnum=bbStclovhi+1-bblovlo;
MOL[bba].t[bbb].vnl=bblovnum;
MOL[bba].t[bbb].vllo=bblovlo;
bbinfo.lovnum=bblovnum;
bbinfo.lovlo=bblovlo;
/* allocate memory for the simsets in the transition structure */
MOL[bba].t[bbb].nS=bbhivnum*bblovnum;
MOL[bba].t[bbb].S=(simset*)calloc(MOL[bba].t[bbb].nS,sizeof(simset)); 
MOL[bba].t[bbb].fS=(int*)calloc(MOL[bba].t[bbb].nS,sizeof(int)); 
if(bbdebugflag<DEBUG){
fprintf(DBG,"MOL[%d].t[%d].vnl=%d, MOL[%d].t[%d].vllo=%d ",bba,bbb,\
		MOL[bba].t[bbb].vnl,bba,bbb,MOL[bba].t[bbb].vllo);
fprintf(DBG,"MOL[%d].t[%d].nS=%d\n",bba,bbb,MOL[bba].t[bbb].nS);
fprintf(DBG,"\tMOL[%d].t[%d].Nhi=%s, MOL[%d].t[%d].Nlo=%s\n",bba,bbb,\
		MOL[bba].t[bbb].Nhi,bba,bbb,MOL[bba].t[bbb].Nlo);
fflush(DBG);
}

bbHV=bbhvnnum+bbhvnlo; 
/* get maximum K value present */
for(bbf=bbhvnlo;bbf<bbHV;bbf++){ /* check native populations */ 
if(bbdebugflag<(DEBUG)){
fprintf(DBG,"In native population loop \n");
fprintf(DBG,"\tMOL[%d].t[%d].Nhi=%s, MOL[%d].t[%d].Nlo=%s\n",bba,bbb,\
		MOL[bba].t[bbb].Nhi,bba,bbb,MOL[bba].t[bbb].Nlo);
fflush(DBG);
}
	bbk=MOL[bba].s[bbSthi].r[bbf-bbhvnlo].k;
	if(bbk>bbmaxhiK){ 
		bbmaxhiK=bbk;
		}
	}
if(bbdebugflag<(DEBUG)){
fprintf(DBG,"Native:  bbe=%d, bbf=%d, bbmaxhiK=%d\n",bbe,bbf,bbmaxhiK);
fprintf(DBG,"\tMOL[%d].t[%d].Nhi=%s, MOL[%d].t[%d].Nlo=%s\n",bba,bbb,\
		MOL[bba].t[bbb].Nhi,bba,bbb,MOL[bba].t[bbb].Nlo);
fflush(DBG);
}
bbHV=bbhvcnum+bbhvclo;
for(bbf=bbhvclo;bbf<bbHV;bbf++){ /* check cascade populations */
	bbk=MOL[bba].s[bbSthi].rc[bbf-bbhvclo].k;
	if(bbk>bbmaxhiK){ 
		bbmaxhiK=bbk;
		}
	}
if(bbdebugflag<(DEBUG)){
fprintf(DBG,"Cascade:  bbe=%d, bbf=%d, bbmaxhiK=%d\n",bbe,bbf,bbmaxhiK);
fprintf(DBG,"\tMOL[%d].t[%d].Nhi=%s, MOL[%d].t[%d].Nlo=%s\n",bba,bbb,\
		MOL[bba].t[bbb].Nhi,bba,bbb,MOL[bba].t[bbb].Nlo);
fflush(DBG);
}

/* See if all the needed destination-state v's are calculated.  Check also 
   that enough K/J's are calculated.  Calculate any states that aren't 
   already done.  The following method isn't efficient.  It will calculate 
   more lower-state energy levels than are actually needed. But, this 
   shouldn't be a problem other than doing a few unnecessary calculations. */ 
bbinfo.hiK=bbmaxhiK;
bbmaxhiK++; /* because of the Delta-K = +1 possibility */
if(bbdebugflag<(DEBUG)){
fprintf(DBG,"Before new calculations loop bbmaxhiK=%d\n",bbmaxhiK);
fprintf(DBG,"\tMOL[%d].t[%d].Nhi=%s, MOL[%d].t[%d].Nlo=%s\n",bba,bbb,\
		MOL[bba].t[bbb].Nhi,bba,bbb,MOL[bba].t[bbb].Nlo);
fflush(DBG);
}
/* Calculate any newly needed energy levels. */
/* loop through each lower Omega */
bbLV=bblovlo+bblovnum;
if(bbdebugflag<(DEBUG)){
fprintf(DBG,"bbd=%d, bbLV=%d\n",bbd,bbLV);
fflush(DBG);
}
for(bbf=bblovlo;bbf<bbLV;bbf++){/* check low-state v-levels */
	bbff=bbd*30+bbf;/* Index for the array of rc rotational sets */ 
if(bbdebugflag<(DEBUG)){
fprintf(DBG,"bbf=%d, bbff=%d, bbdum=%d\n",bbf,bbff,bbdum);
fprintf(DBG,"\tMOL[%d].t[%d].Nhi=%s, MOL[%d].t[%d].Nlo=%s\n",bba,bbb,\
		MOL[bba].t[bbb].Nhi,bba,bbb,MOL[bba].t[bbb].Nlo);
fflush(DBG);
}
	if(MOL[bba].s[bbStlo].rc[bbff].kc==0){
/* this v has never been calculated, so call a function to calculate it. */ 
	if(CB[bbCblo].we!=0){
		bbDe=4*pow(CB[bbCblo].Be,3)/(CB[bbCblo].we*CB[bbCblo].we);
		bbBvDv=(CB[bbCblo].Be-CB[bbCblo].ae*(bbf+0.5))/\
			(2*(bbDe+CB[bbCblo].beta*(bbf+0.5)));
		MOL[bba].s[bbStlo].rc[bbff].Kdissoc=\
			floor((sqrt(1+4*bbBvDv)-1)/2);
		}
	else MOL[bba].s[bbStlo].rc[bbff].Kdissoc=RAND_MAX;
		if(bbmaxhiK<MOL[bba].s[bbStlo].rc[bbff].Kdissoc){
			case_b_cascade(&MOL[bba].s[bbStlo].rc[bbff],\
				bbCblo,bbf,0,bbmaxhiK,bbinfo.nlJ);
			}
		else{
			case_b_cascade(&MOL[bba].s[bbStlo].rc[bbff],bbCblo,\
			bbf,0,(int)(MOL[bba].s[bbStlo].rc[bbff].Kdissoc+1),\
			bbinfo.nlJ);
			} 
if(bbdebugflag<(DEBUG)){
fprintf(DBG,"just called case_b_cascade, v never calc'd\n");
fprintf(DBG,"CB[%d].Be=%12.6e; CB[%d].we=%12.6e; CB[%d].ae=%12.6e\n",\
	bbCblo,CB[bbCblo].Be,bbCblo,CB[bbCblo].we,bbCblo,CB[bbCblo].ae);
fprintf(DBG,"bbf=%d, bbDe=%12.6e, bbBvDv=%12.6e, ",bbf,bbDe,bbBvDv);
fprintf(DBG,"MOL[%d].s[%d].rc[%d].Jdissoc=%12.6e\n",bba,bbStlo,bbff,\
		MOL[bba].s[bbStlo].rc[bbff].Kdissoc); 
fflush(DBG);
fprintf(DBG,"j was zero; just called cascade \n");
fprintf(DBG,"MOL[%d].s[%d].rc[%d].kc=%d\n",bba,bbStlo,bbff,\
		MOL[bba].s[bbStlo].rc[bbff].kc);
fprintf(DBG,"\tMOL[%d].t[%d].Nhi=%s, MOL[%d].t[%d].Nlo=%s\n",bba,bbb,\
		MOL[bba].t[bbb].Nhi,bba,bbb,MOL[bba].t[bbb].Nlo);
fflush(DBG);
}
		} 
	else{
/* This v has been calculated.  See if enough J's are calculated */
		if(MOL[bba].s[bbStlo].rc[bbff].kc<bbmaxhiK){
		if(bbmaxhiK<MOL[bba].s[bbStlo].rc[bbff].Kdissoc){
if(bbdebugflag<(DEBUG)){
fprintf(DBG,"j<maxJ; bbmaxhiK < Kdissoc; about to call cascade \n");
fprintf(DBG,"bbmaxhiK=%d; MOL[%d].s[%d].rc[%d].Kdissoc=%12.6e\n",\
	bbmaxhiK,bba,bbStlo,bbff,MOL[bba].s[bbStlo].rc[bbff].Kdissoc);
fprintf(DBG,"MOL[%d].s[%d].rc[%d].kc=%d\n",bba,bbStlo,bbff,\
		MOL[bba].s[bbStlo].rc[bbff].kc);
fprintf(DBG,"\tMOL[%d].t[%d].Nhi=%s, MOL[%d].t[%d].Nlo=%s\n",bba,bbb,\
		MOL[bba].t[bbb].Nhi,bba,bbb,MOL[bba].t[bbb].Nlo);
fflush(DBG);
}
			case_b_cascade(&MOL[bba].s[bbStlo].rc[bbff],bbCblo,\
				bbf,MOL[bba].s[bbStlo].rc[bbff].kc,\
				bbmaxhiK,bbinfo.nlJ);
			}
		else{
if(bbdebugflag<(DEBUG)){
fprintf(DBG,"j<maxJ; about to call cascade \n");
fprintf(DBG,"j<maxJ; bbmaxhiK > Kdissoc; about to call cascade \n");
fprintf(DBG,"bbmaxhiK=%d; MOL[%d].s[%d].rc[%d].Kdissoc=%12.6e\n",\
	bbmaxhiK,bba,bbStlo,bbff,MOL[bba].s[bbStlo].rc[bbff].Kdissoc);
fprintf(DBG,"MOL[%d].s[%d].rc[%d].kc=%d\n",bba,bbStlo,bbff,\
		MOL[bba].s[bbStlo].rc[bbff].kc);
fprintf(DBG,"\tMOL[%d].t[%d].Nhi=%s, MOL[%d].t[%d].Nlo=%s\n",bba,bbb,\
		MOL[bba].t[bbb].Nhi,bba,bbb,MOL[bba].t[bbb].Nlo);
fflush(DBG);
}
			case_b_cascade(&MOL[bba].s[bbStlo].rc[bbff],bbCblo,\
			bbf,MOL[bba].s[bbStlo].rc[bbff].kc,\
		(int)(MOL[bba].s[bbStlo].rc[bbff].Kdissoc+1),bbinfo.nlJ);
			}
if(bbdebugflag<(DEBUG)){
fprintf(DBG,"j was <maxJ; just called cascade \n");
fprintf(DBG,"MOL[%d].s[%d].rc[%d].kc=%d\n",bba,bbStlo,bbff,\
		MOL[bba].s[bbStlo].rc[bbff].kc);
fprintf(DBG,"\tMOL[%d].t[%d].Nhi=%s, MOL[%d].t[%d].Nlo=%s\n",bba,bbb,\
		MOL[bba].t[bbb].Nhi,bba,bbb,MOL[bba].t[bbb].Nlo);
fflush(DBG);
fprintf(DBG,"just called case_b_cascade for cascade only\n");
fflush(DBG);
}
			}
		} 
	} 
/* Reset values of vclo and vcnum for the next transition */ 
	MOL[bba].s[bbStlo].v[0].vclo=bblovlo;
	MOL[bba].s[bbStlo].v[0].vcnum=bblovnum; 
if(bbdebugflag<(DEBUG)){
fprintf(DBG,"MOL[%d].s[%d].v[%d].vclo=%d ",bba,bbStlo,0,\
		MOL[bba].s[bbStlo].v[0].vclo);
fprintf(DBG,"MOL[%d].s[%d].v[%d].vcnum=%d\n",bba,bbStlo,0,\
		MOL[bba].s[bbStlo].v[0].vcnum);
fprintf(DBG,"\tMOL[%d].t[%d].Nhi=%s, MOL[%d].t[%d].Nlo=%s\n",bba,bbb,\
		MOL[bba].t[bbb].Nhi,bba,bbb,MOL[bba].t[bbb].Nlo);
fflush(DBG);
}

/* Check whether larger selection rules (Lambda, S, +/- in Sigma states) are
   violated.  If they are, write a note to that effect to the DBG file. */ 
if(bbdebugflag<(DEBUG)){
fprintf(DBG,"bbHI[0].L=%d, bbLO[0].L=%d\n",bbHI[0].L,bbLO[0].L);
fflush(DBG);
}
if((abs(bbHI[0].L-bbLO[0].L))>1){
	fprintf(PAR,"WARNING!!! Normal selection rule violated for Delta-");
	fprintf(PAR,"Lambda=0,+/-1 for\n\tTransition %s-->%s\n",\
			MOL[bba].t[bbb].Nhi,MOL[bba].t[bbb].Nlo);
	fprintf(PAR,"\tThis is non-fatal; only a warning.\n");
if(bbdebugflag<(DEBUG)){
fprintf(DBG,"Delta-L selection rule violated.\n");
fflush(DBG);
}
	}
if((fabs(bbHI[0].S-bbLO[0].S))>0.01){ /* <<< an overkill against round-off */
	fprintf(PAR,"WARNING!!! Normal selection rule violated for Delta-");
	fprintf(PAR,"S=0 for\n\tTransition %s-->%s\n",\
			MOL[bba].t[bbb].Nhi,MOL[bba].t[bbb].Nlo);
	fprintf(PAR,"\tThis is non-fatal; only a warning.\n");
if(bbdebugflag<(DEBUG)){
fprintf(DBG,"Delta-S selection rule violated\n");
fflush(DBG);
}
	} 
/* loop through vibration levels starting with the highest in case anyone
   wants to include cascade between v-levels one day.  Same for J's and
   K's later.  Certainly *someone* will want to do that... */
bbT=&MOL[bba].t[bbb]; /* an abbreviation */
if(bbdebugflag<DEBUG){
fprintf(DBG,"bba=%d, bbb=%d, MOL[bba].t[bbb].Nlo=%s\n",bba,bbb,bbT[0].Nlo);
fprintf(DBG,"bbT.Nhi=%s, bbT.Nlo=%s\n",bbT[0].Nhi,bbT[0].Nlo);
fprintf(DBG,"High: MOL[%d].s[%d].v[0].vclo=%d ",bba,bbSthi,\
		MOL[bba].s[bbSthi].v[0].vclo);
fprintf(DBG,"High: MOL[%d].s[%d].v[0].vcnum=%d\n",bba,bbSthi,\
		MOL[bba].s[bbSthi].v[0].vcnum);
fprintf(DBG,"High: MOL[%d].s[%d].v[0].vlo=%d ",bba,bbSthi,\
		MOL[bba].s[bbSthi].v[0].vlo);
fprintf(DBG,"High: MOL[%d].s[%d].v[0].vnum=%d\n",bba,bbSthi,\
		MOL[bba].s[bbSthi].v[0].vnum);
fprintf(DBG,"bbhivnum=%d; bblovnum=%d\n",bbhivnum,bblovnum);
fprintf(DBG,"bbHI[0].L=%d, bbLO[0].L=%d\n",bbHI[0].L,bbLO[0].L);
fflush(DBG);
} 
/* call other functions to calculate transitions for various scenarios */ 
if((bbHI[0].L==0)&&(bbLO[0].L==0)){ /* If Sigma to Sigma */
	if(bbHI[0].p!=bbLO[0].p){ /* For the +<->+, -<->- rule */
fprintf(PAR,"\n\nTRANSITION OMITTED!!!!!!!!!  (see below)\n\n");
fprintf(PAR,"Transition called for Sigma(%d)-->Sigma(%d).\n",\
		bbHI[0].p,bbLO[0].p);
fprintf(PAR,"Molecule is %s; Transition %s-->%s\n",MOL[bba].Mol,\
		bbT[0].Nhi,bbT[0].Nlo);
fprintf(PAR,"This program doesn't know how to break that selection rule.\n");
fprintf(PAR,"The program will continue with this transition omitted.\n\n"); 
fflush(PAR);
		}
	else{ /* check g-u rule if homonuclear and call function */
		if((bbHI[0].g!=0)&&(bbHI[0].g==bbLO[0].g)){
fprintf(PAR,"\n\nTRANSITION OMITTED!!!!!!!!!  (see below)\n\n");
fprintf(PAR,"Transition called for two states that are both either ");
fprintf(PAR,"gerade or ungerade.\n");
fprintf(PAR,"Molecule is %s; Transition %s-->%s\n",MOL[bba].Mol,\
		bbT[0].Nhi,bbT[0].Nlo);
fprintf(PAR,"This program doesn't know how to break that selection rule.\n");
fprintf(PAR,"The program will continue with this transition omitted.\n\n"); 
fflush(PAR);
			}
		else bb_SigmaSigma(bbinfo);
		}
	}
else{
	if((bbHI[0].L==0)||(bbLO[0].L==0)) bb_SigmaOther(bbinfo); /* if one
		or the other state is Sigma, but not both */
	else bb_Other(bbinfo); /* neither state is Sigma */
	} 
return; 
} 

/**************** b-b_SigmaSigma *****************/
  
/*This function calculates transitions between Hund's Cases (b) and (b)
   when both the upper and lower states are sigma states. */

void bb_SigmaSigma(BBTinfo bbSS){
/* see parent function and header file for key to variable names */
/* indexes and dummy variables: */
int bbSSe=0,bbSSee=0,bbSSf=0,bbSSff=0;
int bbSSh=0,bbSSk=0,bbSSfrnh=0,bbSSfrcl=0,bbSSfrch=0; 
int bbSSeh=0,bbSSi=0,bbSSehh=0,bbSSel=0,bbSSell=0,bbSSjii=0;
/* variables for info about the two states */
double bbSSJ=0,bbSSSh=0,bbSSSl=0,bbSStrprob=0,bbSSvhpop=0,bbSShipop=0;
/* Intensity multiplier for satellite bands and FCF dummy variable */
double bbSSSATT=0,bbSSfcf=0,bbSSfcfr=0;
/* variables for calculating Holn-London factors */
int bbSSL=0,bbSSB=0,bbSSdL=0; 
/* variables for spin stats if the state is homonuclear */
int bbSSg=0,bbSSisym=0,bbSSpsym=0;
double bbSSIa=0,bbSScascadetemp=0;
/* a few flags for various purposes */
int bbSSflga=0,bbSSflgb=0,bbSSfrnhf=0,bbSSfrchf=0;
/* variables for creating and writing to output files */
char bbSSfile[1000];
FILE *BBSSFILE;
Case_b_stateinfo *bbSSHI,*bbSSLO;

if(bbSSdebugflag<DEBUG){
fprintf(DBG,"bbSS.hvnlo=%d, bbSS.hvnnum=%d\n",\
		bbSS.hvnlo,bbSS.hvnnum);
fprintf(DBG,"bbSS.hvclo=%d, bbSS.hvcnum=%d\n",\
		bbSS.hvclo,bbSS.hvcnum);
fprintf(DBG,"bbSS.lvnlo=%d, bbSS.lvnnum=%d,\n",\
		bbSS.lvnlo,bbSS.lvnnum);
fprintf(DBG,"bbSS.lvclo=%d, bbSS.lvcnum=%d\n",\
		bbSS.lvclo,bbSS.lvcnum); 
fprintf(DBG,"bbSS.hivlo=%d, bbSS.hivnum=%d,\n",\
		bbSS.hivlo,bbSS.hivnum);
fprintf(DBG,"bbSS.lovlo=%d, bbSS.lovnum=%d\n",\
		bbSS.lovlo,bbSS.lovnum); 
fflush(DBG);
} 
bbSSHI=bbSS.Ch;
bbSSLO=bbSS.Cl;
bbSSSh=bbSSHI[0].S;
bbSSSl=bbSSLO[0].S;
bbSSisym=bbSS.Sh[0].Isymm;
bbSSpsym=bbSS.Sh[0].pmsymm;
bbSShipop=bbSS.Sh[0].pop;
bbSStrprob=bbSS.T[0].P[0];
/* check for homonuclear information */
if(bbSSHI[0].g!=0){
	bbSSg=bbSSHI[0].g;
	bbSSIa=bbSSHI[0].I/(bbSSHI[0].I+1);
	}
/* Just to be explicit about it... */
bbSSL=bbSSdL=0; /* Upper-state Lambda is zero and Delta-Lambda is zero. */

/* start loop down through high vib levels */
for(bbSSe=(bbSS.hivlo+bbSS.hivnum-1);bbSSe>=bbSS.hivlo;bbSSe--){ 
	bbSSee=(bbSSe-bbSS.hivlo)*bbSS.lovnum; /* simset posn. in 2D */
	bbSSfrnh=(bbSSe-bbSS.hvnlo); /* high state native vib position */
	bbSSfrch=bbSSe; /* high state cascade vib position */
/* these flags (bbSSfrnhf and bbSSfrchf) tell if this high vib state is 
   populated natively, by cascade, or both.  They will be used later */ 
	if((bbSSe>=bbSS.hvnlo)&&(bbSSe<(bbSS.hvnlo+bbSS.hvnnum))){
		bbSSfrnhf=0;
		bbSSvhpop=bbSS.Vh[0].p[bbSSe-bbSS.hvnlo];
		}
	else bbSSfrnhf=-1; 
	if((bbSSe>=bbSS.hvclo)&&(bbSSe<(bbSS.hvclo+bbSS.hvcnum))){
		bbSSfrchf=0;
		}
	else bbSSfrchf=-1; 
if(bbSSdebugflag<DEBUG){
fprintf(DBG,"bbSSfrnhf=%d, bbSSvhpop=%f, bbSSfrchf=%d, ",\
		bbSSfrnhf,bbSSvhpop,bbSSfrchf);
fprintf(DBG,"bbSSe=%d, bbSSee=%d\nbbSSfrnh=%d, bbSSfrch=%d, ",\
		bbSSe,bbSSee,bbSSfrnh,bbSSfrch);
fflush(DBG);
} 
/* start loop down through low vib levels */
for(bbSSf=(bbSS.lovlo+bbSS.lovnum-1);bbSSf>=bbSS.lovlo;bbSSf--){ 
	bbSSfcf=bbSS.T[0].v[bbSSe].fcfn[bbSSf];
	if(bbSSfcf!=0){bbSSfcfr=bbSS.T[0].v[bbSSe].fcfc[bbSSf]/bbSSfcf;}
	bbSSff=bbSSee+bbSSf-bbSS.lovlo; /* position in last dimension */
	bbSSfrcl=bbSSf;
if(bbSSdebugflag<DEBUG){
fprintf(DBG,"bbSSf=%d, bbSS.lovlo=%d, bbSS.lovnum=%d \n",\
		bbSSf,bbSS.lovlo,bbSS.lovnum);
fprintf(DBG,"bbSS.nT=%d, bbSS.hiK=%d \n",bbSS.nT,bbSS.hiK);
fflush(DBG);
}
	bbSS.T[0].S[bbSSff].n=bbSS.nT*bbSS.hiK;
	bbSS.T[0].S[bbSSff].f=\
		(double*)calloc(bbSS.T[0].S[bbSSff].n,sizeof(double));
	bbSS.T[0].S[bbSSff].ni=\
		(double*)calloc(bbSS.T[0].S[bbSSff].n,sizeof(double));
	bbSS.T[0].S[bbSSff].ci=\
		(double*)calloc(bbSS.T[0].S[bbSSff].n,sizeof(double));

bbSS.Rl=&bbSS.Sl[0].rc[bbSSfrcl]; 
if(bbSSdebugflag<DEBUG){
fprintf(DBG,"bbSSff=%d, bbSSfrcl=%d, bbSS.T[0].S[bbSSff].n=%d\n",\
		bbSSff,bbSSfrcl,bbSS.T[0].S[bbSSff].n);
fflush(DBG);
}
/* start with highest J value (same reason as before), and loop down looking 
   for lower J's to which to transit.  assign intensities.  If both states are
   Omega=0, then assign zero intensity to transitions for J=0<->J=0.  */ 
	if(bbSSfrnhf==0){/* if this v is natively populated.*/
sprintf(bbSSfile,"%s_molecules/%s/%s--%s_v%d--v%d_NAT.dat",\
	PREF,MOL[bbSS.m].Mol,bbSS.T[0].Nhi,bbSS.T[0].Nlo,bbSSe,bbSSf);
		BBSSFILE=fopen(bbSSfile,"w");
		if(BBSSFILE==NULL){
printf("Error opening transition sub-file %s.  Exiting.\n",bbSSfile);
			exit(1);
			}
fprintf(BBSSFILE,"# File created by %s.\n# File contains ",PROGRAM_NAME);
fprintf(BBSSFILE,"P and R");
fprintf(BBSSFILE," branches from the simulation titled %s\n",PREF);
fprintf(BBSSFILE,"# THESE INTENSITIES ARE DUE TO USER-SPECIFIED POPULATIONS");
fprintf(BBSSFILE," ONLY -- NO CASCADE\n");
fprintf(BBSSFILE,"# Transition is Hund's Case (b) to Hund's Case (b).\n");
fprintf(BBSSFILE,"# Molecule %s, \tHigh state: %s, v=%d\n",\
		MOL[bbSS.m].Mol,bbSS.T[0].Nhi,bbSSe);
fprintf(BBSSFILE,"#\t\tLow state: %s, v=%d\n# Columns are: ",\
		bbSS.T[0].Nlo,bbSSf);
fprintf(BBSSFILE,"K_hi K_lo J_hi P(J_hi+1)_v P_i  R(J_hi-1)_v R_i\n#\n");
fflush(BBSSFILE);
		bbSS.Rh=&bbSS.Sh[0].r[bbSSfrnh]; 
		if((bbSS.hiK)<bbSS.Sh[0].r[bbSSfrnh].k){
			bbSSk=bbSS.hiK-1;
			}
		else bbSSk=bbSS.Sh[0].r[bbSSfrnh].k-1; 
		if(bbSS.Rl[0].Kdissoc<(bbSSk-1)){ 
fprintf(PAR,"\nWARNING!! State %s of molecule %s has significant ",\
		bbSS.T[0].Nhi,MOL[bbSS.m].Mol);
fprintf(PAR,"rotational population at level %d\n",bbSSk);
fprintf(PAR,"\tBut lower state %s dissociates at level %.1f\n",\
		bbSS.T[0].Nlo,bbSS.Rl[0].Kdissoc);
fprintf(PAR,"It is likely that the bound upper state is transiting to (an)");
fprintf(PAR," unbound lower state(s).\n");
fprintf(PAR,"Look for continuum emissions trailing off to the low-energy ");
fprintf(PAR,"end of the v(hi)=%d to v(lo)=%d\nband in the real ",bbSSe,bbSSf);
fprintf(PAR,"spectrum (not the simulated one)\nCurrently, this program will");
fprintf(PAR," not simulate bound-->unbound spectra.\n\n");
			bbSSk=(int)(bbSS.Rl[0].Kdissoc);
			}
if(bbSSdebugflag<DEBUG){
fprintf(DBG,"bbSSfrnh=%d, bbSSk=%d, ",bbSSfrnh,bbSSk);
fprintf(DBG,"bbSS.Sh[0].r[%d].k=%d\n",bbSSfrnh,bbSS.Sh[0].r[bbSSfrnh].k);
fflush(DBG);
}
/* To make the output file easy to write, I'm cycling by high-state K first,
   then by high-state J.  From that J/K(hi) pair, the program will look
   for a P transition and an R transition to K+1 then to K-1.  If all the
   necessary states exist, it will check to see if there are any symmetry
   or nuclear spin effects to be concerned about.  Taking these into account,
   it will calculate the transition. */
/* loop first by upper-state K-value */
for(bbSSh=bbSSk;bbSSh>=0;bbSSh--){
	bbSSjii=(bbSSh+1)*bbSS.nT-1; /* position in simset array */
	bbSSeh=bbSSh*bbSS.nhJ; /* position in 2D for high EJ & PJ */ 
/* then loop by upper-state J for this K */
	for(bbSSi=(bbSS.nhJ-1);bbSSi>=0;bbSSi--){ 
		bbSSJ=bbSSh+bbSSi-bbSSSh;
		bbSSehh=bbSSeh+bbSSi;
		bbSSflga=0;
		if(bbSSh<bbSSHI[0].L) bbSSflga=1; /* if this K is less than 
			Lambda for the high state, don't proceed */
		if(bbSSJ<0) bbSSflga=1; /* if J is less than zero, don't */ 
if(bbSSdebugflag<DEBUG){
fprintf(DBG,"bbSSh=%d, bbSSi=%d, bbSSSh=%.1f, bbSSHI[0].L=%d, bbSSflga=%d\n",\
		bbSSh,bbSSi,bbSSSh,bbSSHI[0].L,bbSSflga);
fprintf(DBG,"bbSSjii=%d\n",bbSSjii);
fflush(DBG);
} 
		if(bbSSflga==0){ 
			bbSSflgb=0;
/* Check Delta-J=+1 (P) and then Delta-J=-1 (R) for Delta-K=+1. First, check
   to see if the target K exists: */ 
			if((bbSSh+1)<(bbSSLO[0].L)) bbSSflgb=1;
/* Relative position in the low state E/P array */
			bbSSel=(bbSSh+1)*bbSS.nlJ + bbSS.nlJ-bbSS.nhJ + \
				(int)(bbSSSh-bbSSSl) -1 + bbSSi; 
if(bbSSdebugflag<DEBUG){
fprintf(DBG,"\t(bbSSh+1)*bbSS.nlJ=%d, bbSS.nlJ-bbSS.nhJ=%d\n",\
		(bbSSh+1)*bbSS.nlJ,bbSS.nlJ-bbSS.nhJ);
fprintf(DBG,"\t(int)(bbSSSh-bbSSSl)=%d, bbSSi=%d, bbSSel=%d\n",\
		(int)(bbSSSh-bbSSSl),bbSSi,bbSSel);
fflush(DBG);
} 
			if((bbSSel+1)>=(bbSS.Rl[0].kc*bbSS.nlJ)) bbSSflgb=1;
/* Then check if the target value for J+1 is less than or greater than the 
   allowed J's for K+1. */
			if((bbSSJ+1)>(bbSSh+1+bbSSSl)) bbSSflgb=1;
			if((bbSSJ+1)<(bbSSh+1-bbSSSl)) bbSSflgb=1; 
/* Still OK?  Calculate transition */ 
			if(bbSSflgb==0){
				bbSSB=+1;
				bbSSell=bbSSel+1;
				bbSS.T[0].S[bbSSff].f[bbSSjii]=\
					bbSS.Rh[0].EJ[bbSSehh]-\
					bbSS.Rl[0].EJ[bbSSell]; 
if(bbSSdebugflag<DEBUG){
fprintf(DBG,"bbSS.T[0].S[%d].f[%d]=%13.5e ",bbSSff,bbSSjii,\
		bbSS.T[0].S[bbSSff].f[bbSSjii]);
fprintf(DBG,"bbSS.Rh[0].EJ[%d]=%13.5e ",bbSSehh,bbSS.Rh[0].EJ[bbSSehh]);
fprintf(DBG,"bbSS.Rl[0].EJ[%d]=%13.5e\n",bbSSell,bbSS.Rl[0].EJ[bbSSell]);
fflush(DBG);
} 
				bbSS.T[0].S[bbSSff].ni[bbSSjii]=\
					bbSShipop*bbSStrprob*bbSSvhpop*\
					bbSS.Rh[0].PJ[bbSSehh]*bbSSfcf*\
					HL(bbSSB,bbSSJ,bbSSL,bbSSdL);
if(bbSSdebugflag<DEBUG){
fprintf(DBG,"bbSShipop=%13.5e bbSStrprob=%13.5e  ",bbSShipop,bbSStrprob);
fprintf(DBG,"bbSSvhpop=%13.5e bbSSfcf=%13.6e\n",bbSSvhpop,bbSSfcf); 
fprintf(DBG,"bbSS.T[0].S[%d].ni[%d]=%13.5e ",\
		bbSSff,bbSSjii,bbSS.T[0].S[bbSSff].ni[bbSSjii]); 
fprintf(DBG,"bbSS.Rh[0].PJ[%d]=%13.5e\n",bbSSehh,bbSS.Rh[0].PJ[bbSSehh]);
fprintf(DBG,"HL(%d,%.1f,%d,%d)=%12.6e\n",bbSSB,bbSSJ,bbSSL,bbSSdL,\
		HL(bbSSB,bbSSJ,bbSSL,bbSSdL));
fflush(DBG);
} 
				bbSS.Rl[0].CJ[bbSSell]+=bbSSfcfr*\
					bbSS.T[0].S[bbSSff].ni[bbSSjii];
if(bbSSdebugflag<DEBUG){
fprintf(DBG,"bbSS.Rl[0].CJ[%d]=%13.5e\n",bbSSell,bbSS.Rl[0].CJ[bbSSell]);
fflush(DBG);
}
				bbSSjii--;
				} 
if(bbSSdebugflag<DEBUG){
fprintf(DBG,"%d ",bbSSh);
fflush(DBG); 
if(bbSSflgb!=0) fprintf(DBG," %d %.1f M M ",(bbSSh+1),bbSSJ);
else{
fprintf(DBG,"%d %.1f %18.12e %18.12e (1a)\n",(bbSSh+1),\
	bbSSJ,bbSS.T[0].S[bbSSff].f[bbSSjii+1],\
	bbSS.T[0].S[bbSSff].ni[bbSSjii+1]);
	} 
fflush(DBG);
}
fprintf(BBSSFILE,"%d ",bbSSh);
if(bbSSflgb!=0) fprintf(BBSSFILE," %d %.1f M M ",(bbSSh+1),bbSSJ);
else{
fprintf(BBSSFILE,"%d %.1f %18.12e %18.12e ",(bbSSh+1),\
	bbSSJ,bbSS.T[0].S[bbSSff].f[bbSSjii+1],\
	bbSS.T[0].S[bbSSff].ni[bbSSjii+1]);
	}
fflush(BBSSFILE);
/* Now, Delta-J=-1 (R) for Delta-K=+1.  */ 
			bbSSflgb=0; 
			if((bbSSh+1)<(bbSSLO[0].L)) bbSSflgb=1;
			if((bbSSJ-1)>(bbSSh+1+bbSSSl)) bbSSflgb=1;
			if((bbSSJ-1)<(bbSSh+1-bbSSSl)) bbSSflgb=1;
			if((bbSSJ-1)<0) bbSSflgb=1;
			if((bbSSel-1)<0) bbSSflgb=1;
			if(bbSSflgb==0){
				bbSSB=-1;
				bbSSell=bbSSel-1;
				bbSS.T[0].S[bbSSff].f[bbSSjii]=\
					bbSS.Rh[0].EJ[bbSSehh]-\
					bbSS.Rl[0].EJ[bbSSell]; 
if(bbSSdebugflag<DEBUG){
fprintf(DBG,"bbSS.T[0].S[%d].f[%d]=%13.5e ",bbSSff,bbSSjii,\
		bbSS.T[0].S[bbSSff].f[bbSSjii]);
fprintf(DBG,"bbSS.Rh[0].EJ[%d]=%13.5e ",bbSSehh,bbSS.Rh[0].EJ[bbSSehh]);
fprintf(DBG,"bbSS.Rl[0].EJ[%d]=%13.5e\n",bbSSell,bbSS.Rl[0].EJ[bbSSell]);
fflush(DBG);
} 
				bbSS.T[0].S[bbSSff].ni[bbSSjii]=\
					bbSShipop*bbSStrprob*bbSSvhpop*\
					bbSS.Rh[0].PJ[bbSSehh]*bbSSSATT*\
					HL(bbSSB,bbSSJ,bbSSL,bbSSdL)*bbSSfcf;
if(bbSSdebugflag<DEBUG){
fprintf(DBG,"bbSShipop=%13.5e bbSStrprob=%13.5e bbSSvhpop=%13.5e ",\
		bbSShipop,bbSStrprob,bbSSvhpop);
fprintf(DBG,"bbSSfcf=%13.6e\n",bbSSfcf);
fprintf(DBG,"bbSS.T[0].S[%d].ni[%d]=%13.5e ",\
		bbSSff,bbSSjii,bbSS.T[0].S[bbSSff].ni[bbSSjii]); 
fprintf(DBG,"bbSS.Rh[0].PJ[%d]=%13.5e\n",bbSSehh,bbSS.Rh[0].PJ[bbSSehh]);
fprintf(DBG,"HL(%d,%.1f,%d,%d)=%12.6e\n",bbSSB,bbSSJ,bbSSL,bbSSdL,\
		HL(bbSSB,bbSSJ,bbSSL,bbSSdL));
fflush(DBG);
} 
				bbSS.Rl[0].CJ[bbSSell]+=bbSSfcfr*\
					bbSS.T[0].S[bbSSff].ni[bbSSjii];
if(bbSSdebugflag<DEBUG){
fprintf(DBG,"bbSS.Rl[0].CJ[%d]=%13.5e\n",bbSSell,bbSS.Rl[0].CJ[bbSSell]);
}
				bbSSjii--;
				} 
if(bbSSdebugflag<DEBUG){
if(bbSSflgb!=0) fprintf(DBG,"M M (Nat) (1b)\n");
else{
fprintf(DBG,"%18.12e (Z) %18.12e (1b)\n",bbSS.T[0].S[bbSSff].f[bbSSjii+1],\
	bbSS.T[0].S[bbSSff].ni[bbSSjii+1]);
	} 
}
if(bbSSflgb!=0) fprintf(BBSSFILE,"M M\n");
else{
fprintf(BBSSFILE,"%18.12e %18.12e\n",bbSS.T[0].S[bbSSff].f[bbSSjii+1],\
	bbSS.T[0].S[bbSSff].ni[bbSSjii+1]);
	} 
			bbSSflgb=0;
/* Check Delta-J=+1 (P) and then Delta-J=-1 (R) for Delta-K=-1. First, check
   to see if the target K exists: */ 
			if((bbSSh-1)<(bbSSLO[0].L)) bbSSflgb=1;
/* Relative position in the low state E/P array */
			bbSSel=(bbSSh-1)*bbSS.nlJ + bbSS.nlJ-bbSS.nhJ + \
				(int)(bbSSSh-bbSSSl) + 1 + bbSSi; 
if(bbSSdebugflag<DEBUG){
fprintf(DBG,"\t(bbSSh+1)*bbSS.nlJ=%d, bbSS.nlJ-bbSS.nhJ=%d\n",\
		(bbSSh+1)*bbSS.nlJ,bbSS.nlJ-bbSS.nhJ);
fprintf(DBG,"\t(int)(bbSSSh-bbSSSl)=%d, bbSSi=%d, bbSSel=%d\n",\
		(int)(bbSSSh-bbSSSl),bbSSi,bbSSel);
fflush(DBG);
} 
/* Then check if the target value for J+1 is less than or greater than the 
   allowed J's for K-1. */
			if((bbSSJ+1)>(bbSSh-1+bbSSSl)) bbSSflgb=1;
			if((bbSSJ+1)<(bbSSh-1-bbSSSl)) bbSSflgb=1;
			if((bbSSel+1)>=(bbSS.Rl[0].kc*bbSS.nlJ)) bbSSflgb=1;
/* Still OK?  Calculate transition */
			if(bbSSflgb==0){
				bbSSB=+1;
				bbSSell=bbSSel+1;
				bbSS.T[0].S[bbSSff].f[bbSSjii]=\
					bbSS.Rh[0].EJ[bbSSehh]-\
					bbSS.Rl[0].EJ[bbSSell]; 
				bbSS.T[0].S[bbSSff].ni[bbSSjii]=\
					bbSShipop*bbSStrprob*bbSSvhpop*\
					bbSS.Rh[0].PJ[bbSSehh]*bbSSSATT*\
					HL(bbSSB,bbSSJ,bbSSL,bbSSdL)*bbSSfcf;
				bbSS.Rl[0].CJ[bbSSell]+=bbSSfcfr*\
					bbSS.T[0].S[bbSSff].ni[bbSSjii];
if(bbSSdebugflag<DEBUG){
fprintf(DBG,"bbSSel=%d, bbSSi-1=%d, bbSSell=%d\n",bbSSel,bbSSi-1,bbSSell);
fprintf(DBG,"\tbbSS.T[0].S[%d].f[%d]=%12.6e\n",\
		bbSSff,bbSSjii,bbSS.T[0].S[bbSSff].f[bbSSjii]);
fprintf(DBG,"\tbbSS.Rh[0].EJ[%d]=%12.6e\n",bbSSehh,bbSS.Rh[0].EJ[bbSSehh]);
fprintf(DBG,"\tbbSS.Rl[0].EJ[%d]=%12.6e\n",bbSSell,bbSS.Rl[0].EJ[bbSSell]);
fprintf(DBG,"\tbbSS.T[0].S[%d].ni[%d]=%12.6e\n",\
		bbSSff,bbSSjii,bbSS.T[0].S[bbSSff].ni[bbSSjii]);
fprintf(DBG,"\tbbSShipop=%12.6e, bbSStrprob=%12.6e, bbSSvhpop=%12.6e ",\
		bbSShipop,bbSStrprob,bbSSvhpop);
fprintf(DBG,"bbSSfcf=%13.6e\n",bbSSfcf);
fprintf(DBG,"\tbbSS.Rh[0].PJ[%d]=%12.6e\n",bbSSehh,bbSS.Rh[0].PJ[bbSSehh]);
fprintf(DBG,"\tHL(%d,%.1f,%d,%d)=%12.6e\n",bbSSB,bbSSJ,\
		bbSSL,bbSSdL,HL(bbSSB,bbSSJ,bbSSL,bbSSdL));
fprintf(DBG,"\tbbSS.Rl[0].CJ[%d]=%12.6e\n",bbSSell,bbSS.Rl[0].CJ[bbSSell]);
fprintf(DBG,"\tbbSS.T[0].S[%d].ni[%d]=%12.6e\n",\
		bbSSff,bbSSjii,bbSS.T[0].S[bbSSff].ni[bbSSjii]);
}
				bbSSjii--;
				} 
if(bbSSdebugflag<DEBUG){
fprintf(DBG,"%d ",bbSSh);
if(bbSSflgb!=0) fprintf(DBG," %d %.1f M M  (here)",(bbSSh-1),bbSSJ);
else{
fprintf(DBG,"%d %.1f %18.12e %18.12e (here) ",(bbSSh-1),\
	bbSSJ,bbSS.T[0].S[bbSSff].f[bbSSjii+1],\
	bbSS.T[0].S[bbSSff].ni[bbSSjii+1]);
	}
}
fprintf(BBSSFILE,"%d ",bbSSh);
if(bbSSflgb!=0) fprintf(BBSSFILE," %d %.1f M M ",(bbSSh-1),bbSSJ);
else{
fprintf(BBSSFILE,"%d %.1f %18.12e %18.12e ",(bbSSh-1),\
	bbSSJ,bbSS.T[0].S[bbSSff].f[bbSSjii+1],\
	bbSS.T[0].S[bbSSff].ni[bbSSjii+1]);
	}
/* Now, Delta-J=-1 (R) for Delta-K=+1.  */ 
			bbSSflgb=0; 
			if((bbSSh-1)<(bbSSLO[0].L)) bbSSflgb=1;
			if((bbSSJ-1)>(bbSSh-1+bbSSSl)) bbSSflgb=1;
			if((bbSSJ-1)<(bbSSh-1-bbSSSl)) bbSSflgb=1;
			if((bbSSJ-1)<0) bbSSflgb=1;
			if((bbSSel-1)<0) bbSSflgb=1;
if(bbSSdebugflag<DEBUG){
fprintf(DBG,"bbSSh-1=%d  bbSSLO[0].L=%d\n",bbSSh-1,bbSSLO[0].L);
fprintf(DBG,"bbSSh-1=%d, bbSS.nlJ-1=%d, bbSSel=%d\n",\
		bbSSh-1,bbSS.nlJ-1,bbSSel);
fprintf(DBG,"bbSSJ-1=%.1f  bbSSh-1=%d, bbSSSl=%.1f\n",bbSSJ-1,bbSSh-1,bbSSSl);
fprintf(DBG,"bbSSJ-1=%.1f, bbSSh-1=%d, bbSSSl=%.1f\n",bbSSJ-1,bbSSh-1,bbSSSl);
fprintf(DBG,"bbSSflgb=%d\n",bbSSflgb);
}
			if(bbSSflgb==0){
				bbSSB=-1;
				bbSSell=bbSSel-1;
				bbSS.T[0].S[bbSSff].f[bbSSjii]=\
					bbSS.Rh[0].EJ[bbSSehh]-\
					bbSS.Rl[0].EJ[bbSSell]; 
				bbSS.T[0].S[bbSSff].ni[bbSSjii]=\
					bbSShipop*bbSStrprob*bbSSvhpop*\
					bbSS.Rh[0].PJ[bbSSehh]*bbSSfcf*\
					HL(bbSSB,bbSSJ,bbSSL,bbSSdL);
				bbSS.Rl[0].CJ[bbSSell]+=bbSSfcfr*\
					bbSS.T[0].S[bbSSff].ni[bbSSjii]; 
if(bbSSdebugflag<DEBUG){
fprintf(DBG,"bbSSel=%d, bbSSi-1=%d, bbSSell=%d\n",bbSSel,bbSSi-1,bbSSell);
fprintf(DBG,"\tbbSS.T[0].S[%d].f[%d]=%12.6e\n",\
		bbSSff,bbSSjii,bbSS.T[0].S[bbSSff].f[bbSSjii]);
fprintf(DBG,"\tbbSS.Rh[0].EJ[%d]=%12.6e\n",bbSSehh,bbSS.Rh[0].EJ[bbSSehh]);
fprintf(DBG,"\tbbSS.Rl[0].EJ[%d]=%12.6e\n",bbSSell,bbSS.Rl[0].EJ[bbSSell]);
fprintf(DBG,"\tbbSS.T[0].S[%d].ni[%d]=%12.6e\n",\
		bbSSff,bbSSjii,bbSS.T[0].S[bbSSff].ni[bbSSjii]);
fprintf(DBG,"\tbbSShipop=%12.6e, bbSStrprob=%12.6e, bbSSvhpop=%12.6e ",\
		bbSShipop,bbSStrprob,bbSSvhpop);
fprintf(DBG,"bbSSfcf=%13.6e\n",bbSSfcf);
fprintf(DBG,"\tbbSS.Rh[0].PJ[%d]=%12.6e\n",bbSSehh,bbSS.Rh[0].PJ[bbSSehh]);
fprintf(DBG,"\tHL(%d,%.1f,%d,%d)=%12.6e\n",bbSSB,bbSSJ,\
		bbSSL,bbSSdL,HL(bbSSB,bbSSJ,bbSSL,bbSSdL));
fprintf(DBG,"\tbbSS.Rl[0].CJ[%d]=%12.6e\n",bbSSell,bbSS.Rl[0].CJ[bbSSell]);
fprintf(DBG,"\tbbSS.T[0].S[%d].ni[%d]=%12.6e\n",\
		bbSSff,bbSSjii,bbSS.T[0].S[bbSSff].ni[bbSSjii]);
}
				}
if(bbSSdebugflag<DEBUG){
if(bbSSflgb!=0) fprintf(DBG,"M M\n");
else{
fprintf(DBG,"%18.12e %18.12e\n",bbSS.T[0].S[bbSSff].f[bbSSjii],\
	bbSS.T[0].S[bbSSff].ni[bbSSjii]);
	} 
}
if(bbSSflgb!=0) fprintf(BBSSFILE,"M M\n");
else{
fprintf(BBSSFILE,"%18.12e %18.12e\n",bbSS.T[0].S[bbSSff].f[bbSSjii],\
	bbSS.T[0].S[bbSSff].ni[bbSSjii]);
	} 
			} /* close if bbSSflga is zero. */
		} /* close upper-state J-of-K loop */ 
	} /* close upper-state K loop */ 
	fclose(BBSSFILE);
		} /* close if high-state v is natively populated */ 

	if(bbSSfrchf==0){/* if this v is cascade populated.*/
sprintf(bbSSfile,"%s_molecules/%s/%s--%s_v%d--v%d_CAS.dat",\
	PREF,MOL[bbSS.m].Mol,bbSS.T[0].Nhi,bbSS.T[0].Nlo,bbSSe,bbSSf);
		BBSSFILE=fopen(bbSSfile,"w");
		if(BBSSFILE==NULL){
printf("Error opening transition sub-file %s.  Exiting.\n",bbSSfile);
			exit(1);
			}
fprintf(BBSSFILE,"# File created by %s.\n# File contains ",PROGRAM_NAME);
fprintf(BBSSFILE,"P and R");
fprintf(BBSSFILE," branches from the simulation titled %s\n",PREF);
fprintf(BBSSFILE,"# THESE INTENSITIES ARE DUE TO CASCADE ONLY -- ");
fprintf(BBSSFILE,"NO USER-DEFINED POPULATIONS\n");
fprintf(BBSSFILE,"# Transition is Hund's Case (b) to Hund's Case (b).\n");
fprintf(BBSSFILE,"# Molecule %s, \tHigh state: %s, v=%d\n",\
		MOL[bbSS.m].Mol,bbSS.T[0].Nhi,bbSSe);
fprintf(BBSSFILE,"#\t\tLow state: %s, v=%d\n# Columns are: ",\
		bbSS.T[0].Nlo,bbSSf);
fprintf(BBSSFILE,"K_hi K_lo J_hi P(J_hi+1)_v P_i  R(J_hi-1)_v R_i\n#\n");

		bbSS.Rh=&bbSS.Sh[0].rc[bbSSfrch]; 
		if((bbSS.hiK)<bbSS.Sh[0].rc[bbSSfrch].kc){
			bbSSk=bbSS.hiK-1;
			}
		else bbSSk=bbSS.Sh[0].rc[bbSSfrch].kc-1; 
		if(bbSS.Rl[0].Kdissoc<(bbSSk-1)){ 
fprintf(PAR,"\nWARNING!! State %s of molecule %s has significant ",\
		bbSS.T[0].Nhi,MOL[bbSS.m].Mol);
fprintf(PAR,"CASCADE rotational population at level %d\n",bbSSk);
fprintf(PAR,"\tBut lower state %s dissociates at level %.1f\n",\
		bbSS.T[0].Nlo,bbSS.Rl[0].Kdissoc);
fprintf(PAR,"It is likely that the bound upper state is transiting to (an)");
fprintf(PAR," unbound lower state(s).\n");
fprintf(PAR,"Look for continuum emissions trailing off to the low-energy ");
fprintf(PAR,"end of the v(hi)=%d to v(lo)=%d\nband in the real ",bbSSe,bbSSf);
fprintf(PAR,"spectrum (not the simulated one)\nCurrently, this program will");
fprintf(PAR," not simulate bound-->unbound spectra.\n\n");
			bbSSk=(int)(bbSS.Rl[0].Kdissoc);
			}
if(bbSSdebugflag<DEBUG){
fprintf(DBG,"bbSSfrch=%d, bbSSk=%d, ",bbSSfrch,bbSSk);
fprintf(DBG,"bbSS.Sh[0].rc[%d].kc=%d\n",bbSSfrch,bbSS.Sh[0].rc[bbSSfrch].kc);
fflush(DBG);
}
/* To make the output file easy to write, I'm cycling by high-state K first,
   then by high-state J.  From that J/K(hi) pair, the program will look
   for a P transition and an R transition to K+1 then to K-1.  If all the
   necessary states exist, it will check to see if there are any symmetry
   or nuclear spin effects to be concerned about.  Taking these into account,
   it will calculate the transition. */
/* loop first by upper-state K-value */
for(bbSSh=bbSSk;bbSSh>=0;bbSSh--){
	bbSSjii=(bbSSh+1)*bbSS.nT-1; /* position in simset array */ 
	bbSSeh=bbSSh*bbSS.nhJ; /* position in 2D for high EJ & PJ */ 
/* then loop by upper-state J for this K */
	for(bbSSi=(bbSS.nhJ-1);bbSSi>=0;bbSSi--){ 
		bbSSJ=bbSSh+bbSSi-bbSSSh;
		bbSSehh=bbSSeh+bbSSi;
		bbSSflga=0;
		if(bbSSh<bbSSHI[0].L) bbSSflga=1; /* if this K is less than 
			Lambda for the high state, don't proceed */
		if(bbSSJ<0) bbSSflga=1; /* if J is less than zero, don't */ 
if(bbSSdebugflag<DEBUG){
fprintf(DBG,"bbSSh=%d, bbSSi=%d, bbSSSh=%.1f, bbSSHI[0].L=%d, bbSSflga=%d\n",\
		bbSSh,bbSSi,bbSSSh,bbSSHI[0].L,bbSSflga);
fflush(DBG);
} 
		if(bbSSflga==0){
			bbSSflgb=0;
/* Check Delta-J=+1 (P) and then Delta-J=-1 (R) for Delta-K=+1. First, check
   to see if the target K exists: */ 
			if((bbSSh+1)<(bbSSLO[0].L)) bbSSflgb=1;
/* Relative position in the low state E/P array */
			bbSSel=(bbSSh+1)*bbSS.nlJ + bbSS.nlJ-bbSS.nhJ + \
				(int)(bbSSSh-bbSSSl) -1 + bbSSi; 
if(bbSSdebugflag<DEBUG){
fprintf(DBG,"\t(bbSSh+1)*bbSS.nlJ=%d, bbSS.nlJ-bbSS.nhJ=%d\n",\
		(bbSSh+1)*bbSS.nlJ,bbSS.nlJ-bbSS.nhJ);
fprintf(DBG,"\t(int)(bbSSSh-bbSSSl)=%d, bbSSi=%d, bbSSel=%d\n",\
		(int)(bbSSSh-bbSSSl),bbSSi,bbSSel);
fflush(DBG);
} 
			if((bbSSel+1)>=(bbSS.Rl[0].kc*bbSS.nlJ)) bbSSflgb=1;
/* Then check if the target value for J+1 is less than or greater than the 
   allowed J's for K+1. */
			if((bbSSJ+1)>(bbSSh+1+bbSSSl)) bbSSflgb=1;
			if((bbSSJ+1)<(bbSSh+1-bbSSSl)) bbSSflgb=1; 
/* Still OK?  Calculate transition */
			if(bbSSflgb==0){
				bbSSB=+1;
				bbSSell=bbSSel+1;
				bbSS.T[0].S[bbSSff].f[bbSSjii]=\
					bbSS.Rh[0].EJ[bbSSehh]-\
					bbSS.Rl[0].EJ[bbSSell]; 
				bbSScascadetemp=bbSStrprob*bbSSfcf*\
					bbSS.Rh[0].CJ[bbSSehh]*\
					HL(bbSSB,bbSSJ,bbSSL,bbSSdL);
if(bbSSdebugflag<DEBUG){
fprintf(DBG,"bbSS.T[0].S[%d].f[%d]=%13.5e ",bbSSff,bbSSjii,\
		bbSS.T[0].S[bbSSff].f[bbSSjii]);
fprintf(DBG,"bbSS.Rh[0].EJ[%d]=%13.5e ",bbSSehh,bbSS.Rh[0].EJ[bbSSehh]);
fprintf(DBG,"bbSS.Rl[0].EJ[%d]=%13.5e\n",bbSSell,bbSS.Rl[0].EJ[bbSSell]);
fprintf(DBG,"bbSStrprob=%13.5e, bbSSfcf=%13.6e",bbSStrprob,bbSSfcf);
fprintf(DBG,"bbSScascadetemp=%13.5e\n",bbSScascadetemp);
fprintf(DBG,"bbSS.Rh[0].CJ[%d]=%13.5e, ",bbSSehh,bbSS.Rh[0].CJ[bbSSehh]);
fprintf(DBG,"HL(%d,%.1f,%d,%d)=%13.5e\n",bbSSB,bbSSJ,bbSSL,bbSSdL,\
		HL(bbSSB,bbSSJ,bbSSL,bbSSdL));
fflush(DBG);
} 
				bbSS.T[0].S[bbSSff].ci[bbSSjii]+=\
					bbSScascadetemp;
				bbSS.Rl[0].CJ[bbSSell]+=bbSSfcfr*\
					bbSScascadetemp;
if(bbSSdebugflag<DEBUG){
fprintf(DBG,"bbSS.T[0].S[%d].ci[%d]=%12.6e bbSS.Rl[0].CJ[%d]=%12.6e\n",\
	bbSSff,bbSSjii,bbSS.T[0].S[bbSSff].ci[bbSSjii],\
	bbSSell,bbSS.Rl[0].CJ[bbSSell]);
}
				bbSSjii--;
				} 
if(bbSSdebugflag<DEBUG){
fprintf(DBG,"%d ",bbSSh);
if(bbSSflgb!=0) fprintf(DBG," %d %.1f (1c) M M\n",(bbSSh+1),bbSSJ);
else{
fprintf(DBG,"%d %.1f %18.12e %18.12e (1c)\n",(bbSSh+1),\
	bbSSJ,bbSS.T[0].S[bbSSff].f[bbSSjii+1],\
	bbSS.T[0].S[bbSSff].ci[bbSSjii+1]);
	} 
}
fprintf(BBSSFILE,"%d ",bbSSh);
if(bbSSflgb!=0) fprintf(BBSSFILE," %d %.1f M M ",(bbSSh+1),bbSSJ);
else{
fprintf(BBSSFILE,"%d %.1f %18.12e %18.12e ",(bbSSh+1),\
	bbSSJ,bbSS.T[0].S[bbSSff].f[bbSSjii+1],\
	bbSS.T[0].S[bbSSff].ci[bbSSjii+1]);
	}
/* Now, Delta-J=-1 (R) for Delta-K=+1.  */ 
			bbSSflgb=0; 
			if((bbSSh+1)<(bbSSLO[0].L)) bbSSflgb=1;
			if((bbSSJ-1)>(bbSSh+1+bbSSSl)) bbSSflgb=1;
			if((bbSSJ-1)<(bbSSh+1-bbSSSl)) bbSSflgb=1;
			if((bbSSJ-1)<0) bbSSflgb=1;
			if((bbSSel-1)<0) bbSSflgb=1;
			if(bbSSflgb==0){
				bbSSB=-1;
				bbSSell=bbSSel-1;
				bbSS.T[0].S[bbSSff].f[bbSSjii]=\
					bbSS.Rh[0].EJ[bbSSehh]-\
					bbSS.Rl[0].EJ[bbSSell]; 
				bbSScascadetemp=bbSStrprob*bbSSfcf*\
					bbSS.Rh[0].CJ[bbSSehh]*bbSSSATT*\
					HL(bbSSB,bbSSJ,bbSSL,bbSSdL);
if(bbSSdebugflag<DEBUG){
fprintf(DBG,"bbSS.T[0].S[%d].f[%d]=%13.5e ",bbSSff,bbSSjii,\
		bbSS.T[0].S[bbSSff].f[bbSSjii]);
fprintf(DBG,"bbSS.Rh[0].EJ[%d]=%13.5e ",bbSSehh,bbSS.Rh[0].EJ[bbSSehh]);
fprintf(DBG,"bbSS.Rl[0].EJ[%d]=%13.5e\n",bbSSell,bbSS.Rl[0].EJ[bbSSell]);
fprintf(DBG,"bbSStrprob=%13.5e, bbSSfcf=%13.6e",bbSStrprob,bbSSfcf);
fprintf(DBG,"bbSScascadetemp=%13.5e\n",bbSScascadetemp);
fprintf(DBG,"bbSS.Rh[0].CJ[%d]=%13.5e, ",bbSSehh,bbSS.Rh[0].CJ[bbSSehh]);
fprintf(DBG,"HL(%d,%.1f,%d,%d)=%13.5e\n",bbSSB,bbSSJ,bbSSL,bbSSdL,\
		HL(bbSSB,bbSSJ,bbSSL,bbSSdL));
fflush(DBG);
}
				bbSS.T[0].S[bbSSff].ci[bbSSjii]+=\
					bbSScascadetemp;
				bbSS.Rl[0].CJ[bbSSell]+=bbSSfcfr*\
					bbSScascadetemp;
if(bbSSdebugflag<DEBUG){
fprintf(DBG,"bbSS.T[0].S[%d].ci[%d]=%12.6e bbSS.Rl[0].CJ[%d]=%12.6e\n",\
	bbSSff,bbSSjii,bbSS.T[0].S[bbSSff].ci[bbSSjii],\
	bbSSell,bbSS.Rl[0].CJ[bbSSell]);
}
				bbSSjii--;
				}
if(bbSSdebugflag<DEBUG){
if(bbSSflgb!=0) fprintf(DBG,"M M\n");
else{
fprintf(DBG,"%18.12e %18.12e (1d)\n",bbSS.T[0].S[bbSSff].f[bbSSjii+1],\
	bbSS.T[0].S[bbSSff].ci[bbSSjii+1]);
	} 
}
if(bbSSflgb!=0) fprintf(BBSSFILE,"M M\n");
else{
fprintf(BBSSFILE,"%18.12e %18.12e\n",bbSS.T[0].S[bbSSff].f[bbSSjii+1],\
	bbSS.T[0].S[bbSSff].ci[bbSSjii+1]);
	} 
			bbSSflgb=0;
/* Check Delta-J=+1 (P) and then Delta-J=-1 (R) for Delta-K=-1. First, check
   to see if the target K exists: */ 
			if((bbSSh-1)<(bbSSLO[0].L)) bbSSflgb=1;
/* Relative position in the low state E/P array */
			bbSSel=(bbSSh-1)*bbSS.nlJ + bbSS.nlJ-bbSS.nhJ + \
				(int)(bbSSSh-bbSSSl) +1 + bbSSi; 
if(bbSSdebugflag<DEBUG){
fprintf(DBG,"\t(bbSSh+1)*bbSS.nlJ=%d, bbSS.nlJ-bbSS.nhJ=%d\n",\
		(bbSSh+1)*bbSS.nlJ,bbSS.nlJ-bbSS.nhJ);
fprintf(DBG,"\t(int)(bbSSSh-bbSSSl)=%d, bbSSi=%d, bbSSel=%d\n",\
		(int)(bbSSSh-bbSSSl),bbSSi,bbSSel);
fflush(DBG);
} 
			if((bbSSel+1)>=(bbSS.Rl[0].kc*bbSS.nlJ)) bbSSflgb=1;
/* Then check if the target value for J+1 is less than or greater than the 
   allowed J's for K-1. */
			if((bbSSJ+1)>(bbSSh-1+bbSSSl)) bbSSflgb=1;
			if((bbSSJ+1)<(bbSSh-1-bbSSSl)) bbSSflgb=1;
/* Still OK?  Calculate transition */
			if(bbSSflgb==0){
				bbSSB=+1;
				bbSSell=bbSSel+1;
				bbSS.T[0].S[bbSSff].f[bbSSjii]=\
					bbSS.Rh[0].EJ[bbSSehh]-\
					bbSS.Rl[0].EJ[bbSSell]; 
				bbSScascadetemp=bbSStrprob*bbSSfcf*\
					bbSS.Rh[0].CJ[bbSSehh]*bbSSSATT*\
					HL(bbSSB,bbSSJ,bbSSL,bbSSdL);
if(bbSSdebugflag<DEBUG){
fprintf(DBG,"bbSS.T[0].S[%d].f[%d]=%13.5e ",bbSSff,bbSSjii,\
		bbSS.T[0].S[bbSSff].f[bbSSjii]);
fprintf(DBG,"bbSS.Rh[0].EJ[%d]=%13.5e ",bbSSehh,bbSS.Rh[0].EJ[bbSSehh]);
fprintf(DBG,"bbSS.Rl[0].EJ[%d]=%13.5e\n",bbSSell,bbSS.Rl[0].EJ[bbSSell]);
fprintf(DBG,"bbSStrprob=%13.5e, bbSSfcf=%13.6e ",bbSStrprob,bbSSfcf);
fprintf(DBG,"bbSScascadetemp=%13.5e\n",bbSScascadetemp);
fprintf(DBG,"bbSS.Rh[0].CJ[%d]=%13.5e, ",bbSSehh,bbSS.Rh[0].CJ[bbSSehh]);
fprintf(DBG,"HL(%d,%.1f,%d,%d)=%13.5e\n",bbSSB,bbSSJ,bbSSL,bbSSdL,\
		HL(bbSSB,bbSSJ,bbSSL,bbSSdL));
fflush(DBG);
} 
				bbSS.T[0].S[bbSSff].ci[bbSSjii]+=\
					bbSScascadetemp;
				bbSS.Rl[0].CJ[bbSSell]+=bbSSfcfr*\
					bbSScascadetemp;
if(bbSSdebugflag<DEBUG){
fprintf(DBG,"bbSS.T[0].S[%d].ci[%d]=%12.6e bbSS.Rl[0].CJ[%d]=%12.6e\n",\
	bbSSff,bbSSjii,bbSS.T[0].S[bbSSff].ci[bbSSjii],\
	bbSSell,bbSS.Rl[0].CJ[bbSSell]);
}
				bbSSjii--;
				} 
if(bbSSdebugflag<DEBUG){
fprintf(DBG,"%d ",bbSSh);
if(bbSSflgb!=0) fprintf(DBG," %d %.1f M M ",(bbSSh-1),bbSSJ);
else{
fprintf(DBG,"%d %.1f %18.12e %18.12e (1e)\n",(bbSSh-1),\
	bbSSJ,bbSS.T[0].S[bbSSff].f[bbSSjii+1],\
	bbSS.T[0].S[bbSSff].ci[bbSSjii+1]);
	}
}
fprintf(BBSSFILE,"%d ",bbSSh);
if(bbSSflgb!=0) fprintf(BBSSFILE," %d %.1f M M ",(bbSSh-1),bbSSJ);
else{
fprintf(BBSSFILE,"%d %.1f %18.12e %18.12e ",(bbSSh-1),\
	bbSSJ,bbSS.T[0].S[bbSSff].f[bbSSjii+1],\
	bbSS.T[0].S[bbSSff].ci[bbSSjii+1]);
	}
/* Now, Delta-J=-1 (R) for Delta-K=-1.  */ 
			bbSSflgb=0; 
			if((bbSSh-1)<(bbSSLO[0].L)) bbSSflgb=1;
			if((bbSSJ-1)>(bbSSh-1+bbSSSl)) bbSSflgb=1;
			if((bbSSJ-1)<(bbSSh-1-bbSSSl)) bbSSflgb=1;
			if((bbSSJ-1)<0) bbSSflgb=1;
			if((bbSSel-1)<0) bbSSflgb=1;
			if(bbSSflgb==0){
				bbSSB=-1;
				bbSSell=bbSSel-1;
				bbSS.T[0].S[bbSSff].f[bbSSjii]=\
					bbSS.Rh[0].EJ[bbSSehh]-\
					bbSS.Rl[0].EJ[bbSSell]; 
				bbSScascadetemp=bbSStrprob*bbSSfcf*\
					bbSS.Rh[0].CJ[bbSSehh]*\
					HL(bbSSB,bbSSJ,bbSSL,bbSSdL);
if(bbSSdebugflag<DEBUG){
fprintf(DBG,"bbSS.T[0].S[%d].f[%d]=%13.5e ",bbSSff,bbSSjii,\
		bbSS.T[0].S[bbSSff].f[bbSSjii]);
fprintf(DBG,"bbSS.Rh[0].EJ[%d]=%13.5e ",bbSSehh,bbSS.Rh[0].EJ[bbSSehh]);
fprintf(DBG,"bbSS.Rl[0].EJ[%d]=%13.5e\n",bbSSell,bbSS.Rl[0].EJ[bbSSell]);
fprintf(DBG,"bbSStrprob=%13.5e, bbSSfcf=%13.6e ",bbSStrprob,bbSSfcf);
fprintf(DBG,"bbSScascadetemp=%13.5e\n",bbSScascadetemp);
fprintf(DBG,"bbSS.Rh[0].CJ[%d]=%13.5e, ",bbSSehh,bbSS.Rh[0].CJ[bbSSehh]);
fprintf(DBG,"HL(%d,%.1f,%d,%d)=%13.5e\n",bbSSB,bbSSJ,bbSSL,bbSSdL,\
		HL(bbSSB,bbSSJ,bbSSL,bbSSdL));
fflush(DBG);
}
				bbSS.T[0].S[bbSSff].ci[bbSSjii]+=\
					bbSScascadetemp;
				bbSS.Rl[0].CJ[bbSSell]+=bbSSfcfr*\
					bbSScascadetemp; 
if(bbSSdebugflag<DEBUG){
fprintf(DBG,"bbSS.T[0].S[%d].ci[%d]=%12.6e bbSS.Rl[0].CJ[%d]=%12.6e\n",\
	bbSSff,bbSSjii,bbSS.T[0].S[bbSSff].ci[bbSSjii],\
	bbSSell,bbSS.Rl[0].CJ[bbSSell]);
}
				}
if(bbSSdebugflag<DEBUG){
if(bbSSflgb!=0) fprintf(DBG,"M M\n");
else{
fprintf(DBG,"%18.12e %18.12e (1f)\n",bbSS.T[0].S[bbSSff].f[bbSSjii],\
	bbSS.T[0].S[bbSSff].ci[bbSSjii]);
	} 
}
if(bbSSflgb!=0) fprintf(BBSSFILE,"M M\n");
else{
fprintf(BBSSFILE,"%18.12e %18.12e\n",bbSS.T[0].S[bbSSff].f[bbSSjii],\
	bbSS.T[0].S[bbSSff].ci[bbSSjii]);
	} 
			} /* close if bbSSflga is zero. */
		} /* close upper-state J-of-K loop */ 
	} /* close upper-state K loop */ 
	fclose(BBSSFILE);
		} /* close if high-state v is cascade populated */ 
	} /* close low-state v */
	} /* close high-state v */
/* check to see if this is all to do... */
return; 
}

/**************** b-b_SigmaOther *****************/
  
/* This function calculates transitions between Hund's Cases (b) and (b)
   when both the upper and lower states are sigma states. */

void bb_SigmaOther(BBTinfo bbSO){
/* see parent function and header file for key to variable names */
/* indexes and dummy variables: */
int bbSOe=0,bbSOee=0,bbSOf=0,bbSOff=0;
int bbSOh=0,bbSOk=0,bbSOfrnh=0,bbSOfrcl=0,bbSOfrch=0; 
int bbSOeh=0,bbSOehh=0,bbSOi=0,bbSOel=0,bbSOell=0,bbSOjii=0;
/* variables for info about the two states */
double bbSOJ=0,bbSOSh=0,bbSOSl=0,bbSOtrprob=0,bbSOvhpop=0,bbSOhipop=0;
/* Intensity multiplier for satellite bands and FCF dummy variable*/
double bbSOSATT=0,bbSOfcf=0,bbSOfcfr=0;
/* variables for calculating Holn-London factors */
int bbSOL=0,bbSOB=0,bbSOdL=0;
/* variables for spin stats if the state is homonuclear */
int bbSOg=0,bbSOisym=0,bbSOpsym=0;
double bbSOIa=0,bbSOcascadetemp=0;
/* a few flags for various purposes */
int bbSOflga=0,bbSOflgb=0,bbSOfrnhf=0,bbSOfrchf=0;
/* variables for creating and writing to output files */
char bbSOfile[1000];
FILE *BBSOFILE;
Case_b_stateinfo *bbSOHI,*bbSOLO;

if(bbSOdebugflag<DEBUG){
fprintf(DBG,"bbSO.hvnlo=%d, bbSO.hvnnum=%d\n",\
		bbSO.hvnlo,bbSO.hvnnum);
fprintf(DBG,"bbSO.hvclo=%d, bbSO.hvcnum=%d\n",\
		bbSO.hvclo,bbSO.hvcnum);
fprintf(DBG,"bbSO.lvnlo=%d, bbSO.lvnnum=%d,\n",\
		bbSO.lvnlo,bbSO.lvnnum);
fprintf(DBG,"bbSO.lvclo=%d, bbSO.lvcnum=%d\n",\
		bbSO.lvclo,bbSO.lvcnum); 
fprintf(DBG,"bbSO.hivlo=%d, bbSO.hivnum=%d,\n",\
		bbSO.hivlo,bbSO.hivnum);
fprintf(DBG,"bbSO.lovlo=%d, bbSO.lovnum=%d\n",\
		bbSO.lovlo,bbSO.lovnum); 
fflush(DBG);
} 
bbSOHI=bbSO.Ch; 
bbSOLO=bbSO.Cl;
bbSOSh=bbSOHI[0].S;
bbSOSl=bbSOLO[0].S;
bbSOhipop=bbSO.Sh[0].pop;
bbSOtrprob=bbSO.T[0].P[0]; 
if(bbSOdebugflag<DEBUG){
fprintf(DBG,"bbSOSh=%.1f, bbSOSl=%.1f,\n",bbSOSh,bbSOSl);
fprintf(DBG,"bbSOhipop=%12.6e, bbSOtrprob=%12.6e\n",bbSOhipop,bbSOtrprob); 
fflush(DBG);
} 

if(bbSOHI[0].L==0){
	bbSOisym=bbSO.Sh[0].Isymm;
	bbSOpsym=bbSO.Sh[0].pmsymm;
	if(bbSOHI[0].g!=0){
		bbSOg=bbSOHI[0].g;
		bbSOIa=bbSOHI[0].I/(bbSOHI[0].I+1);
		}
	}
if(bbSOLO[0].L==0){
	bbSOisym=bbSO.Sl[0].Isymm;
	bbSOpsym=bbSO.Sl[0].pmsymm;
	if(bbSOLO[0].g!=0){
		bbSOg=bbSOLO[0].g;
		bbSOIa=bbSOLO[0].I/(bbSOLO[0].I+1);
		}
	} 
if((bbSOHI[0].L-bbSOLO[0].L)==+1) bbSOdL=+1;
if((bbSOHI[0].L-bbSOLO[0].L)==-1) bbSOdL=-1;
bbSOL=bbSOHI[0].L;

/* start loop down through high vib levels */
for(bbSOe=(bbSO.hivlo+bbSO.hivnum-1);bbSOe>=bbSO.hivlo;bbSOe--){ 
	bbSOee=(bbSOe-bbSO.hivlo)*bbSO.lovnum; /* simset posn. in 2D */
	bbSOfrnh=(bbSOe-bbSO.hvnlo); /* high state native vib position */
	bbSOfrch=bbSOe; /* high state cascade vib position */ 
/* these flags (bbSOfrnhf and bbSOfrchf) tell if this high vib state is 
   populated natively, by cascade, or both.  They will be used later */ 
	if((bbSOe>=bbSO.hvnlo)&&(bbSOe<(bbSO.hvnlo+bbSO.hvnnum))){
		bbSOfrnhf=0; 
		bbSOvhpop=bbSO.Vh[0].p[bbSOe-bbSO.hvnlo]; 
		}
	else bbSOfrnhf=-1; 
	if((bbSOe>=bbSO.hvclo)&&(bbSOe<(bbSO.hvclo+bbSO.hvcnum))){
		bbSOfrchf=0;
		}
	else bbSOfrchf=-1; 
if(bbSOdebugflag<DEBUG){
fprintf(DBG,"bbSOfrnhf=%d, bbSOvhpop=%12.6e, bbSOfrchf=%d, ",\
		bbSOfrnhf,bbSOvhpop,bbSOfrchf);
fprintf(DBG,"bbSOe=%d, bbSOee=%d\nbbSOfrnh=%d, bbSOfrch=%d,\n",\
		bbSOe,bbSOee,bbSOfrnh,bbSOfrch);
fprintf(DBG,"bbSO.Vh[0].p[%d-%d]=%12.6e\n",bbSOe,bbSO.hvnlo,\
	       	bbSO.Vh[0].p[bbSOe-bbSO.hvnlo]); 
fflush(DBG);
} 
/* start loop down through low vib levels */
for(bbSOf=(bbSO.lovlo+bbSO.lovnum-1);bbSOf>=bbSO.lovlo;bbSOf--){ 
	bbSOfcf=bbSO.T[0].v[bbSOe].fcfn[bbSOf];
	if(bbSOfcf!=0){bbSOfcfr=bbSO.T[0].v[bbSOe].fcfc[bbSOf]/bbSOfcf;}
	bbSOff=bbSOee+bbSOf-bbSO.lovlo; /* position in last dimension */
	bbSOfrcl=bbSOf;

if(bbSOdebugflag<DEBUG){
fprintf(DBG,"bbSOf=%d, bbSO.lovlo=%d, bbSO.lovnum=%d, bbSOfcf=%12.6e \n",\
		bbSOf,bbSO.lovlo,bbSO.lovnum,bbSOfcf);
fflush(DBG);
}
bbSOff=bbSOee+bbSOf-bbSO.lovlo; /* simset posn. in final dimension */
	bbSO.T[0].S[bbSOff].n=bbSO.nT*bbSO.hiK;
	bbSO.T[0].S[bbSOff].f=\
		(double*)calloc(bbSO.T[0].S[bbSOff].n,sizeof(double));
	bbSO.T[0].S[bbSOff].ni=\
		(double*)calloc(bbSO.T[0].S[bbSOff].n,sizeof(double));
	bbSO.T[0].S[bbSOff].ci=\
		(double*)calloc(bbSO.T[0].S[bbSOff].n,sizeof(double));

bbSO.Rl=&bbSO.Sl[0].rc[bbSOfrcl]; 
if(bbSOdebugflag<DEBUG){
fprintf(DBG,"bbSOff=%d, bbSOfrcl=%d, bbSO.T[0].S[bbSOff].n=%d\n",\
		bbSOff,bbSOfrcl,bbSO.T[0].S[bbSOff].n);
fprintf(DBG,"bbSO.nT=%d, bbSO.hiK=%d\n",bbSO.nT,bbSO.hiK);
fflush(DBG);
}
/* start with highest J value (same reason as before), and loop down looking 
   for lower J's to which to transit.  Assign intensities.  If both states are 
   Omega=0, then assign zero intensity to transitions for J=0<->J=0.  */ 
	if(bbSOfrnhf==0){/* if this v is natively populated.*/
sprintf(bbSOfile,"%s_molecules/%s/%s--%s_v%d--v%d_NAT.dat",\
	PREF,MOL[bbSO.m].Mol,bbSO.T[0].Nhi,bbSO.T[0].Nlo,bbSOe,bbSOf);
		BBSOFILE=fopen(bbSOfile,"w");
		if(BBSOFILE==NULL){
printf("Error opening transition sub-file %s.  Exiting.\n",bbSOfile);
			exit(1);
			}
fprintf(BBSOFILE,"# File created by %s.\n# File contains ",PROGRAM_NAME);
fprintf(BBSOFILE,"P, Q and R");
fprintf(BBSOFILE," branches from the simulation titled %s\n",PREF);
fprintf(BBSOFILE,"# THESE INTENSITIES ARE DUE TO USER-SPECIFIED POPULATIONS");
fprintf(BBSOFILE," ONLY -- NO CASCADE\n");
fprintf(BBSOFILE,"# Transition is Hund's Case (b) to Hund's Case (b).\n");
fprintf(BBSOFILE,"# Molecule %s, \tHigh state: %s, v=%d\n",\
		MOL[bbSO.m].Mol,bbSO.T[0].Nhi,bbSOe);
fprintf(BBSOFILE,"#\t\tLow state: %s, v=%d\n# Columns are: ",\
		bbSO.T[0].Nlo,bbSOf);
fprintf(BBSOFILE,"K_hi K_lo J_hi P(J_hi+1)_v P_i ");
fprintf(BBSOFILE,"Q(J_hi+0)_v Q_i R(J_hi-1)_v R_i\n#\n");

		bbSO.Rh=&bbSO.Sh[0].r[bbSOfrnh]; 
		if((bbSO.hiK)<bbSO.Sh[0].r[bbSOfrnh].k){
			bbSOk=bbSO.hiK-1;
			}
		else bbSOk=bbSO.Sh[0].r[bbSOfrnh].k-1; 
		if(bbSO.Rl[0].Kdissoc<(bbSOk-1)){ 
fprintf(PAR,"\nWARNING!! State %s of molecule %s has significant ",\
		bbSO.T[0].Nhi,MOL[bbSO.m].Mol);
fprintf(PAR,"rotational population at level %d\n",bbSOk);
fprintf(PAR,"\tBut lower state %s dissociates at level %.1f\n",\
		bbSO.T[0].Nlo,bbSO.Rl[0].Kdissoc);
fprintf(PAR,"It is likely that the bound upper state is transiting to (an)");
fprintf(PAR," unbound lower state(s).\n");
fprintf(PAR,"Look for continuum emissions trailing off to the low-energy ");
fprintf(PAR,"end of the v(hi)=%d to v(lo)=%d\nband in the real ",bbSOe,bbSOf);
fprintf(PAR,"spectrum (not the simulated one)\nCurrently, this program will");
fprintf(PAR," not simulate bound-->unbound spectra.\n\n");
			bbSOk=(int)(bbSO.Rl[0].Kdissoc);
			}
if(bbSOdebugflag<DEBUG){
fprintf(DBG,"bbSOfrnh=%d, bbSOk=%d, ",bbSOfrnh,bbSOk);
fprintf(DBG,"bbSO.Sh[0].r[%d].k=%d\n",bbSOfrnh,bbSO.Sh[0].r[bbSOfrnh].k);
fflush(DBG);
}
/* To make the output file easy to write, I'm cycling by high-state K first,
   then by high-state J.  From that J/K(hi) pair, the program will look
   for a P transition and an R transition to K+1 then to K-1.  If all the
   necessary states exist, it will check to see if there are any symmetry
   or nuclear spin effects to be concerned about.  Taking these into account,
   it will calculate the transition. */
/* loop first by upper-state K-value */
for(bbSOh=bbSOk;bbSOh>=0;bbSOh--){
	bbSOjii=(bbSOh+1)*bbSO.nT-1; /* position in simset array */
	bbSOeh=bbSOh*bbSO.nhJ; /* position in 2D for high EJ & PJ */ 
/* then loop by upper-state J for this K */
	for(bbSOi=(bbSO.nhJ-1);bbSOi>=0;bbSOi--){ 
		bbSOJ=bbSOh+bbSOi-bbSOSh;
		bbSOehh=bbSOeh+bbSOi;
		bbSOflga=0;
		if(bbSOh<bbSOHI[0].L) bbSOflga=1; /* if this K is less than 
			Lambda for the high state, don't proceed */
		if(bbSOJ<0) bbSOflga=1; /* if J is less than zero, don't */ 
if(bbSOdebugflag<DEBUG){
fprintf(DBG,"bbSOh=%d, bbSOi=%d, bbSOSh=%.1f, bbSOHI[0].L=%d, bbSOflga=%d\n",\
		bbSOh,bbSOi,bbSOSh,bbSOHI[0].L,bbSOflga);
fflush(DBG);
} 
		if(bbSOflga==0){ 
/* Check Delta-J=+1 (P), then Delta-J=0 (Q), and then Delta-J=-1 (R) for 
Delta-K=+1. First, check to see if the target K exists: */ 
			bbSOflgb=0;
			if((bbSOh+1)<(bbSOLO[0].L)) bbSOflgb=1;
/* Relative position in the low state E/P array */
			bbSOel=(bbSOh+1)*bbSO.nlJ + bbSO.nlJ-bbSO.nhJ + \
				(int)(bbSOSh-bbSOSl) -1 + bbSOi; 
if(bbSOdebugflag<DEBUG){
fprintf(DBG,"\t(bbSOh+1)*bbSO.nlJ=%d, bbSO.nlJ-bbSO.nhJ=%d\n",\
		(bbSOh+1)*bbSO.nlJ,bbSO.nlJ-bbSO.nhJ);
fprintf(DBG,"\t(int)(bbSOSh-bbSOSl)=%d, bbSOi=%d, bbSOel=%d\n",\
		(int)(bbSOSh-bbSOSl),bbSOi,bbSOel);
fflush(DBG);
} 
			if((bbSOel+1)>=(bbSO.Rl[0].kc*bbSO.nlJ)) bbSOflgb=1;
/* Then check if the target value for J+1 is less than or greater than the 
   allowed J's for K+1. */
			if((bbSOJ+1)>(bbSOh+1+bbSOSl)) bbSOflgb=1;
			if((bbSOJ+1)<(bbSOh+1-bbSOSl)) bbSOflgb=1;
/* Still OK?  Calculate transition */
			if(bbSOflgb==0){
				bbSOB=+1;
				bbSOell=bbSOel+1;
				bbSO.T[0].S[bbSOff].f[bbSOjii]=\
					bbSO.Rh[0].EJ[bbSOehh]-\
					bbSO.Rl[0].EJ[bbSOell]; 
				bbSO.T[0].S[bbSOff].ni[bbSOjii]=\
					bbSOhipop*bbSOtrprob*bbSOvhpop*\
					bbSO.Rh[0].PJ[bbSOehh]*bbSOfcf*\
					HL(bbSOB,bbSOJ,bbSOL,bbSOdL);
if(bbSOdebugflag<DEBUG){
fprintf(DBG,"bbSO.Rh[0].EJ[%d]=%12.6e ",bbSOehh,bbSO.Rh[0].EJ[bbSOehh]);
fprintf(DBG,"bbSO.Rl[0].EJ[%d]=%12.6e ",bbSOell,bbSO.Rl[0].EJ[bbSOell]); 
fprintf(DBG,"bbSO.T[0].S[%d].f[%d]=%12.6e ",bbSOff,bbSOjii,\
		bbSO.T[0].S[bbSOff].f[bbSOjii]);
fprintf(DBG,"bbSOhipop%12.6e, bbSOtrprob%12.6e, bbSOvhpop=%12.6e\n",\
		bbSOhipop,bbSOtrprob,bbSOvhpop);
fprintf(DBG,"bbSOfcf=%13.6e ",bbSOfcf);
fprintf(DBG,"bbSO.Rh[0].PJ[%d]=%12.6e ",bbSOehh,bbSO.Rh[0].PJ[bbSOehh]);
fprintf(DBG,"HL(%d,%.1f,%d,%d)=%12.6e ",bbSOB,bbSOJ,bbSOL,bbSOdL,\
	HL(bbSOB,bbSOJ,bbSOL,bbSOdL));
fprintf(DBG,"bbSO.T[0].S[%d].ni[%d]=%12.6e\n",bbSOff,bbSOjii,\
	bbSO.T[0].S[bbSOff].ni[bbSOjii]);
}
/* Check for nuclear effects and modify intensity if necessary.  No need
   to do this if the high state is Sigma.  */
		if((bbSOg!=0)&&(bbSOLO[0].L==0)){ /* low state is Sigma */
			if((((bbSOh+1)%2)==0)&&(bbSOisym==-1)){
				bbSO.T[0].S[bbSOff].ni[bbSOjii]*=bbSOIa;
if(bbSOdebugflag<DEBUG){
fprintf(DBG,"Even low state is Ia intensity\n");
fflush(DBG);
}
				}
			if((((bbSOh+1)%2)==1)&&(bbSOisym==+1)){
				bbSO.T[0].S[bbSOff].ni[bbSOjii]*=bbSOIa;
if(bbSOdebugflag<DEBUG){
fprintf(DBG,"Even low state is Ia intensity\n");
fflush(DBG);
}
				}
			} 
/* set low-state population for cascade */
				bbSO.Rl[0].CJ[bbSOell]+=bbSOfcfr*\
					bbSO.T[0].S[bbSOff].ni[bbSOjii];
if(bbSOdebugflag<DEBUG){
fprintf(DBG,"bbSO.Rl[0].CJ[%d]=%12.6e, bbSO.T[0].S[%d].ni[%d]=%12.6e\n",\
	bbSOell,bbSO.Rl[0].CJ[bbSOell],bbSOff,\
	bbSOjii,bbSO.T[0].S[bbSOff].ni[bbSOjii]); 
}
				bbSOjii--;
				} 
if(bbSOdebugflag<DEBUG){
fprintf(DBG,"bbSO, K+1, J+1 native\n");
fprintf(DBG,"%d ",bbSOh);
if(bbSOflgb!=0) fprintf(BBSOFILE," %d %.1f M M ",(bbSOh+1),bbSOJ);
else{
fprintf(DBG,"%d %.1f %18.12e %18.12e ",(bbSOh+1),\
	bbSOJ,bbSO.T[0].S[bbSOff].f[bbSOjii+1],\
	bbSO.T[0].S[bbSOff].ni[bbSOjii+1]);
	} 
}
fprintf(BBSOFILE,"%d ",bbSOh);
if(bbSOflgb!=0) fprintf(BBSOFILE," %d %.1f M M ",(bbSOh+1),bbSOJ);
else{
fprintf(BBSOFILE,"%d %.1f %18.12e %18.12e ",(bbSOh+1),\
	bbSOJ,bbSO.T[0].S[bbSOff].f[bbSOjii+1],\
	bbSO.T[0].S[bbSOff].ni[bbSOjii+1]);
	} 
/* Now Delta-J=0 (Q) (for Delta-K=+1)  */ 
			bbSOflgb=0;
			if((bbSOh+1)<(bbSOLO[0].L)) bbSOflgb=1;
/* Then check if the target value for J+1 is less than or greater than the 
   allowed J's for K+1. */
			if((bbSOJ+0)>(bbSOh+1+bbSOSl)) bbSOflgb=1;
			if((bbSOJ+0)<(bbSOh+1-bbSOSl)) bbSOflgb=1;
			if((bbSOJ+0)<0) bbSOflgb=1;
/* Still OK?  Calculate transition */
			if(bbSOflgb==0){
				bbSOB=+0;
				bbSOell=bbSOel+0;
				bbSO.T[0].S[bbSOff].f[bbSOjii]=\
					bbSO.Rh[0].EJ[bbSOehh]-\
					bbSO.Rl[0].EJ[bbSOell]; 
				bbSO.T[0].S[bbSOff].ni[bbSOjii]=\
					bbSOhipop*bbSOtrprob*bbSOvhpop*\
					bbSO.Rh[0].PJ[bbSOehh]*bbSOSATT*\
					HL(bbSOB,bbSOJ,bbSOL,bbSOdL)*bbSOfcf;
if(bbSOdebugflag<DEBUG){
fprintf(DBG,"bbSO.Rh[0].EJ[%d]=%12.6e ",bbSOehh,bbSO.Rh[0].EJ[bbSOehh]);
fprintf(DBG,"bbSO.Rl[0].EJ[%d]=%12.6e ",bbSOell,bbSO.Rl[0].EJ[bbSOell]); 
fprintf(DBG,"bbSO.T[0].S[%d].f[%d]=%12.6e ",bbSOff,bbSOjii,\
		bbSO.T[0].S[bbSOff].f[bbSOjii]);
fprintf(DBG,"bbSOhipop %12.6e, bbSOtrprob %12.6e, bbSOvhpop= %12.6e\n",\
		bbSOhipop,bbSOtrprob,bbSOvhpop);
fprintf(DBG,"bbSOfcf=%13.6e ",bbSOfcf);
fprintf(DBG,"bbSO.Rh[0].PJ[%d]=%12.6e ",bbSOehh,bbSO.Rh[0].PJ[bbSOehh]);
fprintf(DBG,"HL(%d,%.1f,%d,%d)=%12.6e ",bbSOB,bbSOJ,bbSOL,bbSOdL,\
	HL(bbSOB,bbSOJ,bbSOL,bbSOdL));
fprintf(DBG,"bbSO.T[0].S[%d].ni[%d]=%12.6e\n",bbSOff,bbSOjii,\
	bbSO.T[0].S[bbSOff].ni[bbSOjii]);
}
/* Check for nuclear effects and modify intensity if necessary. */
		if((bbSOg!=0)&&(bbSOLO[0].L==0)){ /* low state is Sigma */
			if((((bbSOh+1)%2)==0)&&(bbSOisym==-1)){
				bbSO.T[0].S[bbSOff].ni[bbSOjii]*=bbSOIa;
				}
			if((((bbSOh+1)%2)==1)&&(bbSOisym==+1)){
				bbSO.T[0].S[bbSOff].ni[bbSOjii]*=bbSOIa;
				}
			} 
/* set low-state population for cascade */
				bbSO.Rl[0].CJ[bbSOell]+=bbSOfcfr*\
					bbSO.T[0].S[bbSOff].ni[bbSOjii];
if(bbSOdebugflag<DEBUG){
fprintf(DBG,"bbSO.Rl[0].CJ[%d]=%12.6e, bbSO.T[0].S[%d].ni[%d]=%12.6e\n",\
	bbSOell,bbSO.Rl[0].CJ[bbSOell],bbSOff,\
	bbSOjii,bbSO.T[0].S[bbSOff].ni[bbSOjii]); 
}
				bbSOjii--;
				} 
if(bbSOdebugflag<DEBUG){
fprintf(DBG,"bbSO, K+1, J+0 native\n");
fprintf(DBG,"%d ",bbSOh);
if(bbSOflgb!=0) fprintf(DBG," %d %.1f M M ",(bbSOh+1),bbSOJ);
else{
fprintf(DBG,"%d %.1f %18.12e %18.12e ",(bbSOh+1),\
	bbSOJ,bbSO.T[0].S[bbSOff].f[bbSOjii+1],\
	bbSO.T[0].S[bbSOff].ni[bbSOjii+1]);
	} 
}
if(bbSOflgb!=0) fprintf(BBSOFILE,"  M M ");
else{
fprintf(BBSOFILE," %18.12e %18.12e ",bbSO.T[0].S[bbSOff].f[bbSOjii+1],\
	bbSO.T[0].S[bbSOff].ni[bbSOjii+1]);
	}
/* Now, Delta-J=-1 (R) for Delta-K=+1.  */ 
			bbSOflgb=0; 
			if((bbSOh+1)<(bbSOLO[0].L)) bbSOflgb=1;
			if((bbSOJ-1)>(bbSOh+1+bbSOSl)) bbSOflgb=1;
			if((bbSOJ-1)<(bbSOh+1-bbSOSl)) bbSOflgb=1;
			if((bbSOJ-1)<0) bbSOflgb=1;
			if(bbSOflgb==0){
				bbSOB=-1;
				bbSOell=bbSOel-1;
				bbSO.T[0].S[bbSOff].f[bbSOjii]=\
					bbSO.Rh[0].EJ[bbSOehh]-\
					bbSO.Rl[0].EJ[bbSOell]; 
				bbSO.T[0].S[bbSOff].ni[bbSOjii]=\
					bbSOhipop*bbSOtrprob*bbSOvhpop*\
					bbSO.Rh[0].PJ[bbSOehh]*bbSOSATT*\
					HL(bbSOB,bbSOJ,bbSOL,bbSOdL)*bbSOfcf;
if(bbSOdebugflag<DEBUG){
fprintf(DBG,"bbSO.Rh[0].EJ[%d]=%12.6e ",bbSOehh,bbSO.Rh[0].EJ[bbSOehh]);
fprintf(DBG,"bbSO.Rl[0].EJ[%d]=%12.6e ",bbSOell,bbSO.Rl[0].EJ[bbSOell]); 
fprintf(DBG,"bbSO.T[0].S[%d].f[%d]=%12.6e ",bbSOff,bbSOjii,\
		bbSO.T[0].S[bbSOff].f[bbSOjii]);
fprintf(DBG,"bbSOhipop%12.6e, bbSOtrprob%12.6e, bbSOvhpop=%12.6e\n",\
		bbSOhipop,bbSOtrprob,bbSOvhpop);
fprintf(DBG,"bbSOfcf=%12.6e ",bbSOfcf);
fprintf(DBG,"bbSO.Rh[0].PJ[%d]=%12.6e ",bbSOehh,bbSO.Rh[0].PJ[bbSOehh]);
fprintf(DBG,"HL(%d,%.1f,%d,%d)=%12.6e ",bbSOB,bbSOJ,bbSOL,bbSOdL,\
	HL(bbSOB,bbSOJ,bbSOL,bbSOdL));
fprintf(DBG,"bbSO.T[0].S[%d].ni[%d]=%12.6e\n",bbSOff,bbSOjii,\
	bbSO.T[0].S[bbSOff].ni[bbSOjii]);
}
		if((bbSOg!=0)&&(bbSOLO[0].L==0)){ /* low state is Sigma */
			if((((bbSOh+1)%2)==0)&&(bbSOisym==-1)){
				bbSO.T[0].S[bbSOff].ni[bbSOjii]*=bbSOIa;
				}
			if((((bbSOh+1)%2)==1)&&(bbSOisym==+1)){
				bbSO.T[0].S[bbSOff].ni[bbSOjii]*=bbSOIa;
				}
			} 
				bbSO.Rl[0].CJ[bbSOell]+=bbSOfcfr*\
					bbSO.T[0].S[bbSOff].ni[bbSOjii];
if(bbSOdebugflag<DEBUG){
fprintf(DBG,"bbSO.Rl[0].CJ[%d]=%12.6e, bbSO.T[0].S[%d].ni[%d]=%12.6e\n",\
	bbSOell,bbSO.Rl[0].CJ[bbSOell],bbSOff,\
	bbSOjii,bbSO.T[0].S[bbSOff].ni[bbSOjii]); 
}
				bbSOjii--;
				}
if(bbSOdebugflag<DEBUG){
fprintf(DBG,"bbSO, K+1, J-1 native\n");
if(bbSOflgb!=0) fprintf(DBG,"M M\n");
else{
fprintf(DBG,"%18.12e %18.12e\n",bbSO.T[0].S[bbSOff].f[bbSOjii+1],\
	bbSO.T[0].S[bbSOff].ni[bbSOjii+1]);
	} 
}
if(bbSOflgb!=0) fprintf(BBSOFILE,"M M\n");
else{
fprintf(BBSOFILE,"%18.12e %18.12e\n",bbSO.T[0].S[bbSOff].f[bbSOjii+1],\
	bbSO.T[0].S[bbSOff].ni[bbSOjii+1]);
	} 

/* Check Delta-J=+1 (P), then Delta-J=0 (Q), and then Delta-J=-1 (R) for 
Delta-K=+0. First, check to see if the target K exists: */ 
			bbSOflgb=0;
			if((bbSOh+0)<(bbSOLO[0].L)) bbSOflgb=1;
/* Relative position in the low state E/P array */
			bbSOel=(bbSOh+0)*bbSO.nlJ + bbSO.nlJ-bbSO.nhJ + \
				(int)(bbSOSh-bbSOSl) +0 + bbSOi; 
if(bbSOdebugflag<DEBUG){
fprintf(DBG,"\t(bbSOh+0)*bbSO.nlJ=%d, bbSO.nlJ-bbSO.nhJ=%d\n",\
		(bbSOh+0)*bbSO.nlJ,bbSO.nlJ-bbSO.nhJ);
fprintf(DBG,"\t(int)(bbSOSh-bbSOSl)=%d, bbSOi=%d, bbSOel=%d\n",\
		(int)(bbSOSh-bbSOSl),bbSOi,bbSOel);
fflush(DBG);
} 
			if((bbSOel+1)>=(bbSO.Rl[0].kc*bbSO.nlJ)) bbSOflgb=1;
/* Then check if the target value for J+1 is less than or greater than the 
   allowed J's for K+0. */
			if((bbSOJ+1)>(bbSOh+0+bbSOSl)) bbSOflgb=1;
			if((bbSOJ+1)<(bbSOh+0-bbSOSl)) bbSOflgb=1;
/* Still OK?  Calculate transition */
			if(bbSOflgb==0){
				bbSOB=+1;
				bbSOell=bbSOel+1;
				bbSO.T[0].S[bbSOff].f[bbSOjii]=\
					bbSO.Rh[0].EJ[bbSOehh]-\
					bbSO.Rl[0].EJ[bbSOell]; 
				bbSO.T[0].S[bbSOff].ni[bbSOjii]=\
					bbSOhipop*bbSOtrprob*bbSOvhpop*\
					bbSO.Rh[0].PJ[bbSOehh]*bbSOSATT*\
					HL(bbSOB,bbSOJ,bbSOL,bbSOdL)*bbSOfcf;
if(bbSOdebugflag<DEBUG){
fprintf(DBG,"bbSO.Rh[0].EJ[%d]=%12.6e ",bbSOehh,bbSO.Rh[0].EJ[bbSOehh]);
fprintf(DBG,"bbSO.Rl[0].EJ[%d]=%12.6e ",bbSOell,bbSO.Rl[0].EJ[bbSOell]); 
fprintf(DBG,"bbSO.T[0].S[%d].f[%d]=%12.6e ",bbSOff,bbSOjii,\
		bbSO.T[0].S[bbSOff].f[bbSOjii]);
fprintf(DBG,"bbSOhipop=%12.6e, bbSOtrprob=%12.6e, bbSOvhpop=%12.6e\n",\
		bbSOhipop,bbSOtrprob,bbSOvhpop);
fprintf(DBG,"bbSOfcf=%12.6e ",bbSOfcf);
fprintf(DBG,"bbSO.Rh[0].PJ[%d]=%12.6e ",bbSOehh,bbSO.Rh[0].PJ[bbSOehh]);
fprintf(DBG,"HL(%d,%.1f,%d,%d)=%12.6e ",bbSOB,bbSOJ,bbSOL,bbSOdL,\
	HL(bbSOB,bbSOJ,bbSOL,bbSOdL));
fprintf(DBG,"bbSO.T[0].S[%d].ni[%d]=%12.6e\n",bbSOff,bbSOjii,\
	bbSO.T[0].S[bbSOff].ni[bbSOjii]);
}
/* Check for nuclear effects and modify intensity if necessary. */
		if((bbSOg!=0)&&(bbSOLO[0].L==0)){ /* low state is Sigma */
			if((((bbSOh+0)%2)==0)&&(bbSOisym==-1)){
				bbSO.T[0].S[bbSOff].ni[bbSOjii]*=bbSOIa;
				}
			if((((bbSOh+0)%2)==1)&&(bbSOisym==+1)){
				bbSO.T[0].S[bbSOff].ni[bbSOjii]*=bbSOIa;
				}
			} 
/* set low-state population for cascade */
				bbSO.Rl[0].CJ[bbSOell]+=bbSOfcfr*\
					bbSO.T[0].S[bbSOff].ni[bbSOjii];
if(bbSOdebugflag<DEBUG){
fprintf(DBG,"bbSO.Rl[0].CJ[%d]=%12.6e, bbSO.T[0].S[%d].ni[%d]=%12.6e\n",\
	bbSOell,bbSO.Rl[0].CJ[bbSOell],bbSOff,\
	bbSOjii,bbSO.T[0].S[bbSOff].ni[bbSOjii]); 
}
				bbSOjii--;
				} 
if(bbSOdebugflag<DEBUG){
fprintf(DBG,"bbSO, K+0, J+1 native\n");
fprintf(DBG,"%d ",bbSOh);
if(bbSOflgb!=0) fprintf(DBG," %d %.1f M M ",(bbSOh+0),bbSOJ);
else{
fprintf(DBG,"%d %.1f %18.12e %18.12e ",(bbSOh+0),\
	bbSOJ,bbSO.T[0].S[bbSOff].f[bbSOjii+1],\
	bbSO.T[0].S[bbSOff].ni[bbSOjii+1]);
	} 
}
fprintf(BBSOFILE,"%d ",bbSOh);
if(bbSOflgb!=0) fprintf(BBSOFILE," %d %.1f M M ",(bbSOh+0),bbSOJ);
else{
fprintf(BBSOFILE,"%d %.1f %18.12e %18.12e ",(bbSOh+0),\
	bbSOJ,bbSO.T[0].S[bbSOff].f[bbSOjii+1],\
	bbSO.T[0].S[bbSOff].ni[bbSOjii+1]);
	}

/* Now Delta-J=0 (Q) (for Delta-K=+0)  */ 
			bbSOflgb=0;
			if((bbSOh+0)<(bbSOLO[0].L)) bbSOflgb=1;
/* Then check if the target value for J+1 is less than or greater than the 
   allowed J's for K+0. */
			if((bbSOJ+0)>(bbSOh+0+bbSOSl)) bbSOflgb=1;
			if((bbSOJ+0)<(bbSOh+0-bbSOSl)) bbSOflgb=1;
/* Still OK?  Calculate transition */
			if(bbSOflgb==0){
				bbSOB=0;
				bbSOell=bbSOel+0;
				bbSO.T[0].S[bbSOff].f[bbSOjii]=\
					bbSO.Rh[0].EJ[bbSOehh]-\
					bbSO.Rl[0].EJ[bbSOell]; 
				bbSO.T[0].S[bbSOff].ni[bbSOjii]=\
					bbSOhipop*bbSOtrprob*bbSOvhpop*\
					bbSO.Rh[0].PJ[bbSOehh]*bbSOfcf*\
					HL(bbSOB,bbSOJ,bbSOL,bbSOdL);
if(bbSOdebugflag<DEBUG){
fprintf(DBG,"bbSO.Rh[0].EJ[%d]=%12.6e ",bbSOehh,bbSO.Rh[0].EJ[bbSOehh]);
fprintf(DBG,"bbSO.Rl[0].EJ[%d]=%12.6e ",bbSOell,bbSO.Rl[0].EJ[bbSOell]); 
fprintf(DBG,"bbSO.T[0].S[%d].f[%d]=%12.6e ",bbSOff,bbSOjii,\
		bbSO.T[0].S[bbSOff].f[bbSOjii]);
fprintf(DBG,"bbSOhipop%12.6e, bbSOtrprob%12.6e, bbSOvhpop=%12.6e\n",\
		bbSOhipop,bbSOtrprob,bbSOvhpop);
fprintf(DBG,"bbSOfcf=%12.6e ",bbSOfcf);
fprintf(DBG,"bbSO.Rh[0].PJ[%d]=%12.6e ",bbSOehh,bbSO.Rh[0].PJ[bbSOehh]);
fprintf(DBG,"HL(%d,%.1f,%d,%d)=%12.6e ",bbSOB,bbSOJ,bbSOL,bbSOdL,\
	HL(bbSOB,bbSOJ,bbSOL,bbSOdL));
fprintf(DBG,"bbSO.T[0].S[%d].ni[%d]=%12.6e\n",bbSOff,bbSOjii,\
	bbSO.T[0].S[bbSOff].ni[bbSOjii]);
}
/* Check for nuclear effects and modify intensity if necessary. */
		if((bbSOg!=0)&&(bbSOLO[0].L==0)){ /* low state is Sigma */
			if((((bbSOh+0)%2)==0)&&(bbSOisym==-1)){
				bbSO.T[0].S[bbSOff].ni[bbSOjii]*=bbSOIa;
				}
			if((((bbSOh+0)%2)==1)&&(bbSOisym==+1)){
				bbSO.T[0].S[bbSOff].ni[bbSOjii]*=bbSOIa;
				}
			} 
/* set low-state population for cascade */
				bbSO.Rl[0].CJ[bbSOell]+=bbSOfcfr*\
					bbSO.T[0].S[bbSOff].ni[bbSOjii];
if(bbSOdebugflag<DEBUG){
fprintf(DBG,"bbSO.Rl[0].CJ[%d]=%12.6e, bbSO.T[0].S[%d].ni[%d]=%12.6e\n",\
	bbSOell,bbSO.Rl[0].CJ[bbSOell],bbSOff,\
	bbSOjii,bbSO.T[0].S[bbSOff].ni[bbSOjii]); 
}
				bbSOjii--;
				} 
if(bbSOdebugflag<DEBUG){
fprintf(DBG,"bbSO, K+0, J+0 native\n");
fprintf(DBG,"%d ",bbSOh);
if(bbSOflgb!=0) fprintf(DBG," M M ");
else{
fprintf(DBG," %18.12e %18.12e ",bbSO.T[0].S[bbSOff].f[bbSOjii+1],\
	bbSO.T[0].S[bbSOff].ni[bbSOjii+1]);
	} 
}
if(bbSOflgb!=0) fprintf(BBSOFILE," M M ");
else{
fprintf(BBSOFILE,"%18.12e %18.12e ",bbSO.T[0].S[bbSOff].f[bbSOjii+1],\
	bbSO.T[0].S[bbSOff].ni[bbSOjii+1]);
	} 
/* Now, Delta-J=-1 (R) for Delta-K=+0.  */ 
			bbSOflgb=0; 
			if((bbSOh+0)<(bbSOLO[0].L)) bbSOflgb=1;
			if((bbSOJ-1)>(bbSOh+0+bbSOSl)) bbSOflgb=1;
			if((bbSOJ-1)<(bbSOh+0-bbSOSl)) bbSOflgb=1;
			if((bbSOJ-1)<0) bbSOflgb=1;
			if(bbSOflgb==0){
				bbSOB=-1;
				bbSOell=bbSOel-1;
				bbSO.T[0].S[bbSOff].f[bbSOjii]=\
					bbSO.Rh[0].EJ[bbSOehh]-\
					bbSO.Rl[0].EJ[bbSOell]; 
				bbSO.T[0].S[bbSOff].ni[bbSOjii]=\
					bbSOhipop*bbSOtrprob*bbSOvhpop*\
					bbSO.Rh[0].PJ[bbSOehh]*bbSOSATT*\
					HL(bbSOB,bbSOJ,bbSOL,bbSOdL)*bbSOfcf;
if(bbSOdebugflag<DEBUG){
fprintf(DBG,"bbSO.Rh[0].EJ[%d]=%12.6e ",bbSOehh,bbSO.Rh[0].EJ[bbSOehh]);
fprintf(DBG,"bbSO.Rl[0].EJ[%d]=%12.6e ",bbSOell,bbSO.Rl[0].EJ[bbSOell]); 
fprintf(DBG,"bbSO.T[0].S[%d].f[%d]=%12.6e ",bbSOff,bbSOjii,\
		bbSO.T[0].S[bbSOff].f[bbSOjii]);
fprintf(DBG,"bbSOhipop=%12.6e, bbSOtrprob=%12.6e, bbSOvhpop=%12.6e\n",\
		bbSOhipop,bbSOtrprob,bbSOvhpop);
fprintf(DBG,"bbSOfcf=%12.6e ",bbSOfcf);
fprintf(DBG,"bbSO.Rh[0].PJ[%d]=%12.6e ",bbSOehh,bbSO.Rh[0].PJ[bbSOehh]);
fprintf(DBG,"HL(%d,%.1f,%d,%d)=%12.6e ",bbSOB,bbSOJ,bbSOL,bbSOdL,\
	HL(bbSOB,bbSOJ,bbSOL,bbSOdL));
fprintf(DBG,"bbSO.T[0].S[%d].ni[%d]=%12.6e\n",bbSOff,bbSOjii,\
	bbSO.T[0].S[bbSOff].ni[bbSOjii]);
}
		if((bbSOg!=0)&&(bbSOLO[0].L==0)){ /* low state is Sigma */
			if((((bbSOh+0)%2)==0)&&(bbSOisym==-1)){
				bbSO.T[0].S[bbSOff].ni[bbSOjii]*=bbSOIa;
				}
			if((((bbSOh+0)%2)==1)&&(bbSOisym==+1)){
				bbSO.T[0].S[bbSOff].ni[bbSOjii]*=bbSOIa;
				}
			} 
				bbSO.Rl[0].CJ[bbSOell]+=bbSOfcfr*\
					bbSO.T[0].S[bbSOff].ni[bbSOjii];
if(bbSOdebugflag<DEBUG){
fprintf(DBG,"bbSO.Rl[0].CJ[%d]=%12.6e, bbSO.T[0].S[%d].ni[%d]=%12.6e\n",\
	bbSOell,bbSO.Rl[0].CJ[bbSOell],bbSOff,\
	bbSOjii,bbSO.T[0].S[bbSOff].ni[bbSOjii]); 
}
				bbSOjii--;
				}
if(bbSOdebugflag<DEBUG){
fprintf(DBG,"bbSO, K+0, J-1 native\n");
if(bbSOflgb!=0) fprintf(DBG,"M M\n");
else{
fprintf(DBG,"%18.12e %18.12e\n",bbSO.T[0].S[bbSOff].f[bbSOjii+1],\
	bbSO.T[0].S[bbSOff].ni[bbSOjii+1]);
	} 
}
if(bbSOflgb!=0) fprintf(BBSOFILE,"M M\n");
else{
fprintf(BBSOFILE,"%18.12e %18.12e\n",bbSO.T[0].S[bbSOff].f[bbSOjii+1],\
	bbSO.T[0].S[bbSOff].ni[bbSOjii+1]);
	} 
/* Check Delta-J=+1 (P), then Delta-J=0 and Delta-J=-1 (R) for Delta-K=-1. 
   First, check to see if the target K exists: */ 
			bbSOflgb=0;
			if((bbSOh-1)<(bbSOLO[0].L)) bbSOflgb=1;
/* Relative position in the low state E/P array */
			bbSOel=(bbSOh-1)*bbSO.nlJ + bbSO.nlJ-bbSO.nhJ + \
				(int)(bbSOSh-bbSOSl) + 1 + bbSOi; 
if(bbSOdebugflag<DEBUG){
fprintf(DBG,"\t(bbSOh-1)*bbSO.nlJ=%d, bbSO.nlJ-bbSO.nhJ=%d\n",\
		(bbSOh-1)*bbSO.nlJ,bbSO.nlJ-bbSO.nhJ);
fprintf(DBG,"\t(int)(bbSOSh-bbSOSl)=%d, bbSOi=%d, bbSOel=%d\n",\
		(int)(bbSOSh-bbSOSl),bbSOi,bbSOel);
fflush(DBG);
} 
			if((bbSOel+1)>=(bbSO.Rl[0].kc*bbSO.nlJ)) bbSOflgb=1;
/* Then check if the target value for J+1 is less than or greater than the 
   allowed J's for K-1. */
			if((bbSOJ+1)>(bbSOh-1+bbSOSl)) bbSOflgb=1;
			if((bbSOJ+1)<(bbSOh-1-bbSOSl)) bbSOflgb=1;
/* Still OK?  Calculate transition */
			if(bbSOflgb==0){
				bbSOB=+1;
				bbSOell=bbSOel+1;
				bbSO.T[0].S[bbSOff].f[bbSOjii]=\
					bbSO.Rh[0].EJ[bbSOehh]-\
					bbSO.Rl[0].EJ[bbSOell]; 
				bbSO.T[0].S[bbSOff].ni[bbSOjii]=\
					bbSOhipop*bbSOtrprob*bbSOvhpop*\
					bbSO.Rh[0].PJ[bbSOehh]*bbSOSATT*\
					HL(bbSOB,bbSOJ,bbSOL,bbSOdL)*bbSOfcf;
if(bbSOdebugflag<DEBUG){
fprintf(DBG,"bbSO.Rh[0].EJ[%d]=%12.6e ",bbSOehh,bbSO.Rh[0].EJ[bbSOehh]);
fprintf(DBG,"bbSO.Rl[0].EJ[%d]=%12.6e ",bbSOell,bbSO.Rl[0].EJ[bbSOell]); 
fprintf(DBG,"bbSO.T[0].S[%d].f[%d]=%12.6e ",bbSOff,bbSOjii,\
		bbSO.T[0].S[bbSOff].f[bbSOjii]);
fprintf(DBG,"bbSOhipop%12.6e, bbSOtrprob%12.6e, bbSOvhpop=%12.6e\n",\
		bbSOhipop,bbSOtrprob,bbSOvhpop);
fprintf(DBG,"bbSOfcf=%12.6e ",bbSOfcf);
fprintf(DBG,"bbSO.Rh[0].PJ[%d]=%12.6e ",bbSOehh,bbSO.Rh[0].PJ[bbSOehh]);
fprintf(DBG,"HL(%d,%.1f,%d,%d)=%12.6e ",bbSOB,bbSOJ,bbSOL,bbSOdL,\
	HL(bbSOB,bbSOJ,bbSOL,bbSOdL));
fprintf(DBG,"bbSO.T[0].S[%d].ni[%d]=%12.6e\n",bbSOff,bbSOjii,\
	bbSO.T[0].S[bbSOff].ni[bbSOjii]);
}
/* Check for nuclear effects and modify intensity if necessary. */
		if((bbSOg!=0)&&(bbSOLO[0].L==0)){ /* low state is Sigma */
			if((((bbSOh-1)%2)==0)&&(bbSOisym==-1)){
				bbSO.T[0].S[bbSOff].ni[bbSOjii]*=bbSOIa;
				}
			if((((bbSOh-1)%2)==1)&&(bbSOisym==+1)){
				bbSO.T[0].S[bbSOff].ni[bbSOjii]*=bbSOIa;
				}
			} 
				bbSO.Rl[0].CJ[bbSOell]+=bbSOfcfr*\
					bbSO.T[0].S[bbSOff].ni[bbSOjii];
if(bbSOdebugflag<DEBUG){
fprintf(DBG,"bbSO.Rl[0].CJ[%d]=%12.6e, bbSO.T[0].S[%d].ni[%d]=%12.6e\n",\
	bbSOell,bbSO.Rl[0].CJ[bbSOell],bbSOff,\
	bbSOjii,bbSO.T[0].S[bbSOff].ni[bbSOjii]); 
}
				bbSOjii--;
				} 
if(bbSOdebugflag<DEBUG){
fprintf(DBG,"bbSO, K-1, J+1 native\n");
fprintf(DBG,"%d ",bbSOh);
if(bbSOflgb!=0) fprintf(DBG," %d %.1f M M ",(bbSOh-1),bbSOJ);
else{
fprintf(DBG,"%d %.1f %18.12e %18.12e ",(bbSOh-1),\
	bbSOJ,bbSO.T[0].S[bbSOff].f[bbSOjii+1],\
	bbSO.T[0].S[bbSOff].ni[bbSOjii+1]);
	}
}
fprintf(BBSOFILE,"%d ",bbSOh);
if(bbSOflgb!=0) fprintf(BBSOFILE," %d %.1f M M ",(bbSOh-1),bbSOJ);
else{
fprintf(BBSOFILE,"%d %.1f %18.12e %18.12e ",(bbSOh-1),\
	bbSOJ,bbSO.T[0].S[bbSOff].f[bbSOjii+1],\
	bbSO.T[0].S[bbSOff].ni[bbSOjii+1]);
	} 
/* Now, Delta-J=+0 (Q) for Delta-K=-1. */ 
			bbSOflgb=0;
			if((bbSOh-1)<(bbSOLO[0].L)) bbSOflgb=1;
/* Then check if the target value for J+1 is less than or greater than the 
   allowed J's for K-1. */
			if((bbSOJ+0)>(bbSOh-1+bbSOSl)) bbSOflgb=1;
			if((bbSOJ+0)<(bbSOh-1-bbSOSl)) bbSOflgb=1;
/* Still OK?  Calculate transition */
			if(bbSOflgb==0){
				bbSOB=+0;
				bbSOell=bbSOel+0;
				bbSO.T[0].S[bbSOff].f[bbSOjii]=\
					bbSO.Rh[0].EJ[bbSOehh]-\
					bbSO.Rl[0].EJ[bbSOell]; 
				bbSO.T[0].S[bbSOff].ni[bbSOjii]=\
					bbSOhipop*bbSOtrprob*bbSOvhpop*\
					bbSO.Rh[0].PJ[bbSOehh]*bbSOSATT*\
					HL(bbSOB,bbSOJ,bbSOL,bbSOdL)*bbSOfcf;
if(bbSOdebugflag<DEBUG){
fprintf(DBG,"bbSO.Rh[0].EJ[%d]=%12.6e ",bbSOehh,bbSO.Rh[0].EJ[bbSOehh]);
fprintf(DBG,"bbSO.Rl[0].EJ[%d]=%12.6e ",bbSOell,bbSO.Rl[0].EJ[bbSOell]); 
fprintf(DBG,"bbSO.T[0].S[%d].f[%d]=%12.6e ",bbSOff,bbSOjii,\
		bbSO.T[0].S[bbSOff].f[bbSOjii]);
fprintf(DBG,"bbSOhipop=%12.6e, bbSOtrprob=%12.6e, bbSOvhpop=%12.6e\n",\
		bbSOhipop,bbSOtrprob,bbSOvhpop);
fprintf(DBG,"bbSOfcf=%12.6e ",bbSOfcf);
fprintf(DBG,"bbSO.Rh[0].PJ[%d]=%12.6e ",bbSOehh,bbSO.Rh[0].PJ[bbSOehh]);
fprintf(DBG,"HL(%d,%.1f,%d,%d)=%12.6e ",bbSOB,bbSOJ,bbSOL,bbSOdL,\
	HL(bbSOB,bbSOJ,bbSOL,bbSOdL));
fprintf(DBG,"bbSO.T[0].S[%d].ni[%d]=%12.6e\n",bbSOff,bbSOjii,\
	bbSO.T[0].S[bbSOff].ni[bbSOjii]);
}
/* Check for nuclear effects and modify intensity if necessary. */
		if((bbSOg!=0)&&(bbSOLO[0].L==0)){ /* low state is Sigma */
			if((((bbSOh-1)%2)==0)&&(bbSOisym==-1)){
				bbSO.T[0].S[bbSOff].ni[bbSOjii]*=bbSOIa;
				}
			if((((bbSOh-1)%2)==1)&&(bbSOisym==+1)){
				bbSO.T[0].S[bbSOff].ni[bbSOjii]*=bbSOIa;
				}
			} 
				bbSO.Rl[0].CJ[bbSOell]+=bbSOfcfr*\
					bbSO.T[0].S[bbSOff].ni[bbSOjii];
if(bbSOdebugflag<DEBUG){
fprintf(DBG,"bbSO.Rl[0].CJ[%d]=%12.6e, bbSO.T[0].S[%d].ni[%d]=%12.6e\n",\
	bbSOell,bbSO.Rl[0].CJ[bbSOell],bbSOff,\
	bbSOjii,bbSO.T[0].S[bbSOff].ni[bbSOjii]); 
}
				bbSOjii--;
				} 
if(bbSOdebugflag<DEBUG){
fprintf(DBG,"bbSO, K-1, J+0 native\n");
fprintf(DBG,"%d ",bbSOh);
if(bbSOflgb!=0) fprintf(DBG," %d %.1f M M ",(bbSOh-1),bbSOJ);
else{
fprintf(DBG,"%d %.1f %18.12e %18.12e ",(bbSOh-1),\
	bbSOJ,bbSO.T[0].S[bbSOff].f[bbSOjii+1],\
	bbSO.T[0].S[bbSOff].ni[bbSOjii+1]);
	}
}
if(bbSOflgb!=0) fprintf(BBSOFILE,"M M ");
else{
fprintf(BBSOFILE," %18.12e %18.12e ",bbSO.T[0].S[bbSOff].f[bbSOjii+1],\
	bbSO.T[0].S[bbSOff].ni[bbSOjii+1]);
	} 
/* Now, Delta-J=-1 (R) for Delta-K=-1.  */ 
			bbSOflgb=0; 
			if((bbSOh-1)<(bbSOLO[0].L)) bbSOflgb=1;
			if((bbSOJ-1)>(bbSOh-1+bbSOSl)) bbSOflgb=1;
			if((bbSOJ-1)<(bbSOh-1-bbSOSl)) bbSOflgb=1;
			if((bbSOJ-1)<0) bbSOflgb=1;
			if(bbSOflgb==0){
				bbSOB=-1;
				bbSOell=bbSOel-1;
				bbSO.T[0].S[bbSOff].f[bbSOjii]=\
					bbSO.Rh[0].EJ[bbSOehh]-\
					bbSO.Rl[0].EJ[bbSOell]; 
				bbSO.T[0].S[bbSOff].ni[bbSOjii]=\
					bbSOhipop*bbSOtrprob*bbSOvhpop*\
					bbSO.Rh[0].PJ[bbSOehh]*bbSOfcf*\
				       HL(bbSOB,bbSOJ,bbSOL,bbSOdL);	
if(bbSOdebugflag<DEBUG){
fprintf(DBG,"bbSO.Rh[0].EJ[%d]=%12.6e ",bbSOehh,bbSO.Rh[0].EJ[bbSOehh]);
fprintf(DBG,"bbSO.Rl[0].EJ[%d]=%12.6e ",bbSOell,bbSO.Rl[0].EJ[bbSOell]); 
fprintf(DBG,"bbSO.T[0].S[%d].f[%d]=%12.6e ",bbSOff,bbSOjii,\
		bbSO.T[0].S[bbSOff].f[bbSOjii]);
fprintf(DBG,"bbSOhipop=%12.6e, bbSOtrprob=%12.6e, bbSOvhpop=%12.6e\n",\
		bbSOhipop,bbSOtrprob,bbSOvhpop);
fprintf(DBG,"bbSOfcf=%12.6e ",bbSOfcf);
fprintf(DBG,"bbSO.Rh[0].PJ[%d]=%12.6e ",bbSOehh,bbSO.Rh[0].PJ[bbSOehh]);
fprintf(DBG,"HL(%d,%.1f,%d,%d)=%12.6e ",bbSOB,bbSOJ,bbSOL,bbSOdL,\
	HL(bbSOB,bbSOJ,bbSOL,bbSOdL));
fprintf(DBG,"bbSO.T[0].S[%d].ni[%d]=%12.6e\n",bbSOff,bbSOjii,\
	bbSO.T[0].S[bbSOff].ni[bbSOjii]);
}
		if((bbSOg!=0)&&(bbSOLO[0].L==0)){ /* low state is Sigma */
			if((((bbSOh-1)%2)==0)&&(bbSOisym==-1)){
				bbSO.T[0].S[bbSOff].ni[bbSOjii]*=bbSOIa;
				}
			if((((bbSOh-1)%2)==1)&&(bbSOisym==+1)){
				bbSO.T[0].S[bbSOff].ni[bbSOjii]*=bbSOIa;
				}
			} 
				bbSO.Rl[0].CJ[bbSOell]+=bbSOfcfr*\
					bbSO.T[0].S[bbSOff].ni[bbSOjii]; 
if(bbSOdebugflag<DEBUG){
fprintf(DBG,"bbSO.Rl[0].CJ[%d]=%12.6e, bbSO.T[0].S[%d].ni[%d]=%12.6e\n",\
	bbSOell,bbSO.Rl[0].CJ[bbSOell],bbSOff,\
	bbSOjii,bbSO.T[0].S[bbSOff].ni[bbSOjii]); 
}
				}
if(bbSOdebugflag<DEBUG){
fprintf(DBG,"bbSO, K-1, J-1 native\n");
if(bbSOflgb!=0) fprintf(DBG,"M M\n");
else{
fprintf(DBG,"%18.12e %18.12e\n",bbSO.T[0].S[bbSOff].f[bbSOjii],\
	bbSO.T[0].S[bbSOff].ni[bbSOjii]);
	} 
}
if(bbSOflgb!=0) fprintf(BBSOFILE,"M M\n");
else{
fprintf(BBSOFILE,"%18.12e %18.12e\n",bbSO.T[0].S[bbSOff].f[bbSOjii],\
	bbSO.T[0].S[bbSOff].ni[bbSOjii]);
	} 
			} /* close if bbSOflga is zero. */
		} /* close upper-state J-of-K loop */ 
	} /* close upper-state K loop */ 
	fclose(BBSOFILE);
		} /* close if high-state v is natively populated */ 

	if(bbSOfrchf==0){/* if this v is cascade populated.*/
sprintf(bbSOfile,"%s_molecules/%s/%s--%s_v%d--v%d_CAS.dat",\
	PREF,MOL[bbSO.m].Mol,bbSO.T[0].Nhi,bbSO.T[0].Nlo,bbSOe,bbSOf);
		BBSOFILE=fopen(bbSOfile,"w");
		if(BBSOFILE==NULL){
printf("Error opening transition sub-file %s.  Exiting.\n",bbSOfile);
			exit(1);
			}
fprintf(BBSOFILE,"# File created by %s.\n# File contains ",PROGRAM_NAME);
fprintf(BBSOFILE,"P and R");
fprintf(BBSOFILE," branches from the simulation titled %s\n",PREF);
fprintf(BBSOFILE,"# THESE INTENSITIES ARE DUE TO CASCADE ONLY -- ");
fprintf(BBSOFILE,"NO USER-DEFINED POPULATIONS\n");
fprintf(BBSOFILE,"# Transition is Hund's Case (b) to Hund's Case (b).\n");
fprintf(BBSOFILE,"# Molecule %s, \tHigh state: %s, v=%d\n",\
		MOL[bbSO.m].Mol,bbSO.T[0].Nhi,bbSOe);
fprintf(BBSOFILE,"#\t\tLow state: %s, v=%d\n# Columns are: ",\
		bbSO.T[0].Nlo,bbSOf);
fprintf(BBSOFILE,"K_hi K_lo J_hi P(J_hi+1)_v P_i  R(J_hi-1)_v R_i\n#\n");

		bbSO.Rh=&bbSO.Sh[0].rc[bbSOfrch]; 
		if((bbSO.hiK)<bbSO.Sh[0].rc[bbSOfrch].kc){
			bbSOk=bbSO.hiK-1;
			}
		else bbSOk=bbSO.Sh[0].rc[bbSOfrch].kc-1; 
		if(bbSO.Rl[0].Kdissoc<(bbSOk-1)){ 
fprintf(PAR,"\nWARNING!! State %s of molecule %s has significant ",\
		bbSO.T[0].Nhi,MOL[bbSO.m].Mol);
fprintf(PAR,"CASCADE rotational population at level %d\n",bbSOk);
fprintf(PAR,"\tBut lower state %s dissociates at level %.1f\n",\
		bbSO.T[0].Nlo,bbSO.Rl[0].Kdissoc);
fprintf(PAR,"It is likely that the bound upper state is transiting to (an)");
fprintf(PAR," unbound lower state(s).\n");
fprintf(PAR,"Look for continuum emissions trailing off to the low-energy ");
fprintf(PAR,"end of the v(hi)=%d to v(lo)=%d\nband in the real ",bbSOe,bbSOf);
fprintf(PAR,"spectrum (not the simulated one)\nCurrently, this program will");
fprintf(PAR," not simulate bound-->unbound spectra.\n\n");
			bbSOk=(int)(bbSO.Rl[0].Kdissoc);
			}
if(bbSOdebugflag<DEBUG){
fprintf(DBG,"bbSOfrch=%d, bbSOk=%d, ",bbSOfrch,bbSOk);
fprintf(DBG,"bbSO.Sh[0].rc[%d].kc=%d\n",bbSOfrch,bbSO.Sh[0].rc[bbSOfrch].kc);
fflush(DBG);
}
/* To make the output file easy to write, I'm cycling by high-state K first,
   then by high-state J.  From that J/K(hi) pair, the program will look
   for a P transition and an R transition to K+1 then to K-1.  If all the
   necessary states exist, it will check to see if there are any symmetry
   or nuclear spin effects to be concerned about.  Taking these into account,
   it will calculate the transition. */
/* loop first by upper-state K-value */
for(bbSOh=bbSOk;bbSOh>=0;bbSOh--){ 
	bbSOjii=(bbSOh+1)*bbSO.nT-1; /* position in simset array */
	bbSOeh=bbSOh*bbSO.nhJ; /* position in 2D for high EJ & PJ */ 
/* then loop by upper-state J for this K */
	for(bbSOi=(bbSO.nhJ-1);bbSOi>=0;bbSOi--){ 
		bbSOJ=bbSOh+bbSOi-bbSOSh;
		bbSOehh=bbSOeh+bbSOi;
		bbSOflga=0;
		if(bbSOh<bbSOHI[0].L) bbSOflga=1; /* if this K is less than 
			Lambda for the high state, don't proceed */
		if(bbSOJ<0) bbSOflga=1; /* if J is less than zero, don't */ 
if(bbSOdebugflag<DEBUG){
fprintf(DBG,"bbSOh=%d, bbSOi=%d, bbSOSh=%.1f, bbSOHI[0].L=%d, bbSOflga=%d\n",\
		bbSOh,bbSOi,bbSOSh,bbSOHI[0].L,bbSOflga);
fflush(DBG);
} 
		if(bbSOflga==0){ 
/* Check Delta-J=+1 (P), then Delta-J=0 (Q) and Delta-J=-1 (R) for Delta-K=+1.
   First, check to see if the target K exists: */ 
			bbSOflgb=0;
			if((bbSOh+1)<(bbSOLO[0].L)) bbSOflgb=1;
/* Relative position in the low state E/P array */
			bbSOel=(bbSOh+1)*bbSO.nlJ + bbSO.nlJ-bbSO.nhJ + \
				(int)(bbSOSh-bbSOSl) -1 + bbSOi; 
if(bbSOdebugflag<DEBUG){
fprintf(DBG,"\t(bbSOh+1)*bbSO.nlJ=%d, bbSO.nlJ-bbSO.nhJ=%d\n",\
		(bbSOh+1)*bbSO.nlJ,bbSO.nlJ-bbSO.nhJ);
fprintf(DBG,"\t(int)(bbSOSh-bbSOSl)=%d, bbSOi=%d, bbSOel=%d\n",\
		(int)(bbSOSh-bbSOSl),bbSOi,bbSOel);
fflush(DBG);
} 
			if((bbSOel+1)>=(bbSO.Rl[0].kc*bbSO.nlJ)) bbSOflgb=1;
/* Then check if the target value for J+1 is less than or greater than the 
   allowed J's for K+1. */
			if((bbSOJ+1)>(bbSOh+1+bbSOSl)) bbSOflgb=1;
			if((bbSOJ+1)<(bbSOh+1-bbSOSl)) bbSOflgb=1;
/* Still OK?  Calculate transition */
			if(bbSOflgb==0){
				bbSOB=+1;
				bbSOell=bbSOel+1;
				bbSO.T[0].S[bbSOff].f[bbSOjii]=\
					bbSO.Rh[0].EJ[bbSOehh]-\
					bbSO.Rl[0].EJ[bbSOell]; 
				bbSOcascadetemp=bbSOtrprob*bbSOfcf*\
					bbSO.Rh[0].CJ[bbSOehh]*\
					HL(bbSOB,bbSOJ,bbSOL,bbSOdL);
/* Check for nuclear effects and modify intensity if necessary. */ 
		if((bbSOg!=0)&&(bbSOLO[0].L==0)){ /* low state is Sigma */
			if((((bbSOh+1)%2)==0)&&(bbSOisym==-1)){
				bbSOcascadetemp*=bbSOIa;
				}
			if((((bbSOh+1)%2)==1)&&(bbSOisym==+1)){
				bbSOcascadetemp*=bbSOIa;
				}
			} 
				bbSO.T[0].S[bbSOff].ci[bbSOjii]+=\
					bbSOcascadetemp;
				bbSO.Rl[0].CJ[bbSOell]+=bbSOfcfr*\
					bbSOcascadetemp;
				bbSOjii--;
				} 
if(bbSOdebugflag<DEBUG){
fprintf(DBG,"%d ",bbSOh);
if(bbSOflgb!=0) fprintf(DBG," %d %.1f (1c) M M\n",(bbSOh+1),bbSOJ);
else{
fprintf(DBG,"%d %.1f %18.12e %18.12e (1c)\n",(bbSOh+1),\
	bbSOJ,bbSO.T[0].S[bbSOff].f[bbSOjii+1],\
	bbSO.T[0].S[bbSOff].ci[bbSOjii+1]);
	} 
}
fprintf(BBSOFILE,"%d ",bbSOh);
if(bbSOflgb!=0) fprintf(BBSOFILE," %d %.1f M M ",(bbSOh+1),bbSOJ);
else{
fprintf(BBSOFILE,"%d %.1f %18.12e %18.12e ",(bbSOh+1),\
	bbSOJ,bbSO.T[0].S[bbSOff].f[bbSOjii+1],\
	bbSO.T[0].S[bbSOff].ci[bbSOjii+1]);
	} 
/* Now, Delta-J=0 (Q), Delta-K=+1 */ 
			bbSOflgb=0;
			if((bbSOh+1)<(bbSOLO[0].L)) bbSOflgb=1;
/* Then check if the target value for J+1 is less than or greater than the 
   allowed J's for K+1. */
			if((bbSOJ+0)>(bbSOh+1+bbSOSl)) bbSOflgb=1;
			if((bbSOJ+0)<(bbSOh+1-bbSOSl)) bbSOflgb=1;
			if((bbSOJ+0)<0) bbSOflgb=1;
/* Still OK?  Calculate transition */
			if(bbSOflgb==0){
				bbSOB=+0;
				bbSOell=bbSOel+0;
				bbSO.T[0].S[bbSOff].f[bbSOjii]=\
					bbSO.Rh[0].EJ[bbSOehh]-\
					bbSO.Rl[0].EJ[bbSOell]; 
				bbSOcascadetemp=bbSOtrprob*bbSOfcf*\
					bbSO.Rh[0].CJ[bbSOehh]*bbSOSATT*\
					HL(bbSOB,bbSOJ,bbSOL,bbSOdL);
/* Check for nuclear effects and modify intensity if necessary. */ 
		if((bbSOg!=0)&&(bbSOLO[0].L==0)){ /* low state is Sigma */
			if((((bbSOh+1)%2)==0)&&(bbSOisym==-1)){
				bbSOcascadetemp*=bbSOIa;
				}
			if((((bbSOh+1)%2)==1)&&(bbSOisym==+1)){
				bbSOcascadetemp*=bbSOIa;
				}
			} 
				bbSO.T[0].S[bbSOff].ci[bbSOjii]+=\
					bbSOcascadetemp;
				bbSO.Rl[0].CJ[bbSOell]+=bbSOfcfr*\
					bbSOcascadetemp;
				bbSOjii--;
				} 
if(bbSOdebugflag<DEBUG){
if(bbSOflgb!=0) fprintf(DBG," M M\n");
else{
fprintf(DBG,"%18.12e %18.12e (1c)\n",bbSO.T[0].S[bbSOff].f[bbSOjii+1],\
	bbSO.T[0].S[bbSOff].ci[bbSOjii+1]);
	} 
}
if(bbSOflgb!=0) fprintf(BBSOFILE," M M ");
else{
fprintf(BBSOFILE," %18.12e %18.12e ",bbSO.T[0].S[bbSOff].f[bbSOjii+1],\
	bbSO.T[0].S[bbSOff].ci[bbSOjii+1]);
	} 
/* Now, Delta-J=-1 (R) for Delta-K=+1.  */ 
			bbSOflgb=0; 
			if((bbSOh+1)<(bbSOLO[0].L)) bbSOflgb=1;
			if((bbSOJ-1)>(bbSOh+1+bbSOSl)) bbSOflgb=1;
			if((bbSOJ-1)<(bbSOh+1-bbSOSl)) bbSOflgb=1;
			if((bbSOJ-1)<0) bbSOflgb=1;
			if(bbSOflgb==0){
				bbSOB=-1;
				bbSOell=bbSOel-1;
				bbSO.T[0].S[bbSOff].f[bbSOjii]=\
					bbSO.Rh[0].EJ[bbSOehh]-\
					bbSO.Rl[0].EJ[bbSOell]; 
				bbSOcascadetemp=bbSOtrprob*bbSOfcf*\
					bbSO.Rh[0].CJ[bbSOehh]*bbSOSATT*\
				       HL(bbSOB,bbSOJ,bbSOL,bbSOdL);	
		if((bbSOg!=0)&&(bbSOLO[0].L==0)){ /* low state is Sigma */
			if((((bbSOh+1)%2)==0)&&(bbSOisym==-1)){
				bbSOcascadetemp*=bbSOIa;
				}
			if((((bbSOh+1)%2)==1)&&(bbSOisym==+1)){
				bbSOcascadetemp*=bbSOIa;
				}
			} 
				bbSO.T[0].S[bbSOff].ci[bbSOjii]+=\
					bbSOcascadetemp;
				bbSO.Rl[0].CJ[bbSOell]+=bbSOfcfr*\
					bbSOcascadetemp;
				bbSOjii--;
				}
if(bbSOdebugflag<DEBUG){
if(bbSOflgb!=0) fprintf(DBG,"M M\n");
else{
fprintf(DBG,"%18.12e %18.12e (1d)\n",bbSO.T[0].S[bbSOff].f[bbSOjii+1],\
	bbSO.T[0].S[bbSOff].ci[bbSOjii+1]);
	} 
}
if(bbSOflgb!=0) fprintf(BBSOFILE,"M M\n");
else{
fprintf(BBSOFILE,"%18.12e %18.12e\n",bbSO.T[0].S[bbSOff].f[bbSOjii+1],\
	bbSO.T[0].S[bbSOff].ci[bbSOjii+1]);
	} 

/* Now, again for Delta-K=+0.  Delta-J=+1 (P), then (Q) and (R) */ 
			bbSOflgb=0;
			if((bbSOh+0)<(bbSOLO[0].L)) bbSOflgb=1;
/* Relative position in the low state E/P array */
			bbSOel=(bbSOh+0)*bbSO.nlJ + bbSO.nlJ-bbSO.nhJ + \
				(int)(bbSOSh-bbSOSl) + 0 + bbSOi; 
if(bbSOdebugflag<DEBUG){
fprintf(DBG,"\t(bbSOh+0)*bbSO.nlJ=%d, bbSO.nlJ-bbSO.nhJ=%d\n",\
		(bbSOh+0)*bbSO.nlJ,bbSO.nlJ-bbSO.nhJ);
fprintf(DBG,"\t(int)(bbSOSh-bbSOSl)=%d, bbSOi=%d, bbSOel=%d\n",\
		(int)(bbSOSh-bbSOSl),bbSOi,bbSOel);
fflush(DBG);
} 
			if((bbSOel+1)>=(bbSO.Rl[0].kc*bbSO.nlJ)) bbSOflgb=1;
/* Then check if the target value for J+1 is less than or greater than the 
   allowed J's for K+0. */
			if((bbSOJ+1)>(bbSOh+0+bbSOSl)) bbSOflgb=1;
			if((bbSOJ+1)<(bbSOh+0-bbSOSl)) bbSOflgb=1;
/* Still OK?  Calculate transition */
			if(bbSOflgb==0){
				bbSOB=+1;
				bbSOell=bbSOel+1;
				bbSO.T[0].S[bbSOff].f[bbSOjii]=\
					bbSO.Rh[0].EJ[bbSOehh]-\
					bbSO.Rl[0].EJ[bbSOell]; 
				bbSOcascadetemp=bbSOtrprob*bbSOfcf*\
					bbSO.Rh[0].CJ[bbSOehh]*bbSOSATT*\
					HL(bbSOB,bbSOJ,bbSOL,bbSOdL);
/* Check for nuclear effects and modify intensity if necessary. */ 
		if((bbSOg!=0)&&(bbSOLO[0].L==0)){ /* low state is Sigma */
			if((((bbSOh+0)%2)==0)&&(bbSOisym==-1)){
				bbSOcascadetemp*=bbSOIa;
				}
			if((((bbSOh+0)%2)==1)&&(bbSOisym==+1)){
				bbSOcascadetemp*=bbSOIa;
				}
			} 
				bbSO.T[0].S[bbSOff].ci[bbSOjii]+=\
					bbSOcascadetemp;
				bbSO.Rl[0].CJ[bbSOell]+=bbSOfcfr*\
					bbSOcascadetemp;
				bbSOjii--;
				} 
if(bbSOdebugflag<DEBUG){
fprintf(DBG,"%d ",bbSOh);
if(bbSOflgb!=0) fprintf(DBG," %d %.1f (1c) M M\n",(bbSOh+0),bbSOJ);
else{
fprintf(DBG,"%d %.1f %18.12e %18.12e (1c)\n",(bbSOh+0),\
	bbSOJ,bbSO.T[0].S[bbSOff].f[bbSOjii+1],\
	bbSO.T[0].S[bbSOff].ci[bbSOjii+1]);
	} 
}
fprintf(BBSOFILE,"%d ",bbSOh);
if(bbSOflgb!=0) fprintf(BBSOFILE," %d %.1f M M ",(bbSOh+0),bbSOJ);
else{
fprintf(BBSOFILE,"%d %.1f %18.12e %18.12e ",(bbSOh+0),\
	bbSOJ,bbSO.T[0].S[bbSOff].f[bbSOjii+1],\
	bbSO.T[0].S[bbSOff].ci[bbSOjii+1]);
	} 
/* Now, Delta-J=0 (Q), Delta-K=+0 */ 
			bbSOflgb=0;
			if((bbSOh+0)<(bbSOLO[0].L)) bbSOflgb=1;
/* Then check if the target value for J+1 is less than or greater than the 
   allowed J's for K+0. */
			if((bbSOJ+0)>(bbSOh+0+bbSOSl)) bbSOflgb=1;
			if((bbSOJ+0)<(bbSOh+0-bbSOSl)) bbSOflgb=1;
/* Still OK?  Calculate transition */
			if(bbSOflgb==0){
				bbSOB=+0;
				bbSOell=bbSOel+0;
				bbSO.T[0].S[bbSOff].f[bbSOjii]=\
					bbSO.Rh[0].EJ[bbSOehh]-\
					bbSO.Rl[0].EJ[bbSOell]; 
				bbSOcascadetemp=bbSOtrprob*bbSOfcf*\
					bbSO.Rh[0].CJ[bbSOehh]*\
					HL(bbSOB,bbSOJ,bbSOL,bbSOdL);
/* Check for nuclear effects and modify intensity if necessary. */ 
		if((bbSOg!=0)&&(bbSOLO[0].L==0)){ /* low state is Sigma */
			if((((bbSOh+0)%2)==0)&&(bbSOisym==-1)){
				bbSOcascadetemp*=bbSOIa;
				}
			if((((bbSOh+0)%2)==1)&&(bbSOisym==+1)){
				bbSOcascadetemp*=bbSOIa;
				}
			} 
				bbSO.T[0].S[bbSOff].ci[bbSOjii]+=\
					bbSOcascadetemp;
				bbSO.Rl[0].CJ[bbSOell]+=bbSOfcfr*\
					bbSOcascadetemp;
				bbSOjii--;
				} 
if(bbSOdebugflag<DEBUG){
if(bbSOflgb!=0) fprintf(DBG," M M\n");
else{
fprintf(DBG,"%18.12e %18.12e (1c)\n",bbSO.T[0].S[bbSOff].f[bbSOjii+1],\
	bbSO.T[0].S[bbSOff].ci[bbSOjii+1]);
	} 
}
if(bbSOflgb!=0) fprintf(BBSOFILE," M M ");
else{
fprintf(BBSOFILE," %18.12e %18.12e ",bbSO.T[0].S[bbSOff].f[bbSOjii+1],\
	bbSO.T[0].S[bbSOff].ci[bbSOjii+1]);
	} 
/* Now, Delta-J=-1 (R) for Delta-K=+0.  */ 
			bbSOflgb=0; 
			if((bbSOh+0)<(bbSOLO[0].L)) bbSOflgb=1;
			if((bbSOJ-1)>(bbSOh+0+bbSOSl)) bbSOflgb=1;
			if((bbSOJ-1)<(bbSOh+0-bbSOSl)) bbSOflgb=1;
			if((bbSOJ-1)<0) bbSOflgb=1;
			if(bbSOflgb==0){
				bbSOB=-1;
				bbSOell=bbSOel-1;
				bbSO.T[0].S[bbSOff].f[bbSOjii]=\
					bbSO.Rh[0].EJ[bbSOehh]-\
					bbSO.Rl[0].EJ[bbSOell]; 
				bbSOcascadetemp=bbSOtrprob*bbSOfcf*\
					bbSO.Rh[0].CJ[bbSOehh]*bbSOSATT*\
				       HL(bbSOB,bbSOJ,bbSOL,bbSOdL);	
		if((bbSOg!=0)&&(bbSOLO[0].L==0)){ /* low state is Sigma */
			if((((bbSOh+0)%2)==0)&&(bbSOisym==-1)){
				bbSOcascadetemp*=bbSOIa;
				}
			if((((bbSOh+0)%2)==1)&&(bbSOisym==+1)){
				bbSOcascadetemp*=bbSOIa;
				}
			} 
				bbSO.T[0].S[bbSOff].ci[bbSOjii]+=\
					bbSOcascadetemp;
				bbSO.Rl[0].CJ[bbSOell]+=bbSOfcfr*\
					bbSOcascadetemp;
				bbSOjii--;
				}
if(bbSOdebugflag<DEBUG){
if(bbSOflgb!=0) fprintf(DBG,"M M\n");
else{
fprintf(DBG,"%18.12e %18.12e (1d)\n",bbSO.T[0].S[bbSOff].f[bbSOjii+1],\
	bbSO.T[0].S[bbSOff].ci[bbSOjii+1]);
	} 
}
if(bbSOflgb!=0) fprintf(BBSOFILE,"M M\n");
else{
fprintf(BBSOFILE,"%18.12e %18.12e\n",bbSO.T[0].S[bbSOff].f[bbSOjii+1],\
	bbSO.T[0].S[bbSOff].ci[bbSOjii+1]);
	} 

/* Check Delta-J=+1 (P), then Delta-J=0 (Q) and Delta-J=-1 (R) for Delta-K=-1.
   First, check to see if the target K exists: */ 
			bbSOflgb=0;
			if((bbSOh-1)<(bbSOLO[0].L)) bbSOflgb=1;
/* Relative position in the low state E/P array */
			bbSOel=(bbSOh-1)*bbSO.nlJ + bbSO.nlJ-bbSO.nhJ + \
				(int)(bbSOSh-bbSOSl) + 1 + bbSOi; 
if(bbSOdebugflag<DEBUG){
fprintf(DBG,"\t(bbSOh-1)*bbSO.nlJ=%d, bbSO.nlJ-bbSO.nhJ=%d\n",\
		(bbSOh-1)*bbSO.nlJ,bbSO.nlJ-bbSO.nhJ);
fprintf(DBG,"\t(int)(bbSOSh-bbSOSl)=%d, bbSOi=%d, bbSOel=%d\n",\
		(int)(bbSOSh-bbSOSl),bbSOi,bbSOel);
fflush(DBG);
} 
			if((bbSOel+1)>=(bbSO.Rl[0].kc*bbSO.nlJ)) bbSOflgb=1;
/* Then check if the target value for J+1 is less than or greater than the 
   allowed J's for K+1. */
			if((bbSOJ+1)>(bbSOh-1+bbSOSl)) bbSOflgb=1;
			if((bbSOJ+1)<(bbSOh-1-bbSOSl)) bbSOflgb=1;
/* Still OK?  Calculate transition */
			if(bbSOflgb==0){
				bbSOB=+1;
				bbSOell=bbSOel+1;
				bbSO.T[0].S[bbSOff].f[bbSOjii]=\
					bbSO.Rh[0].EJ[bbSOehh]-\
					bbSO.Rl[0].EJ[bbSOell]; 
				bbSOcascadetemp=bbSOtrprob*bbSOfcf*\
					bbSO.Rh[0].CJ[bbSOehh]*bbSOSATT*\
					HL(bbSOB,bbSOJ,bbSOL,bbSOdL);
/* Check for nuclear effects and modify intensity if necessary. */
		if((bbSOg!=0)&&(bbSOLO[0].L==0)){ /* low state is Sigma */
			if((((bbSOh-1)%2)==0)&&(bbSOisym==-1)){
				bbSOcascadetemp*=bbSOIa;
				}
			if((((bbSOh-1)%2)==1)&&(bbSOisym==+1)){
				bbSOcascadetemp*=bbSOIa;
				}
			} 
				bbSO.T[0].S[bbSOff].ci[bbSOjii]+=\
					bbSOcascadetemp;
				bbSO.Rl[0].CJ[bbSOell]+=bbSOfcfr*\
					bbSOcascadetemp;
				bbSOjii--;
				} 
if(bbSOdebugflag<DEBUG){
fprintf(DBG,"%d ",bbSOh);
if(bbSOflgb!=0) fprintf(DBG," %d %.1f M M ",(bbSOh-1),bbSOJ);
else{
fprintf(DBG,"%d %.1f %18.12e %18.12e (1e)\n",(bbSOh-1),\
	bbSOJ,bbSO.T[0].S[bbSOff].f[bbSOjii+1],\
	bbSO.T[0].S[bbSOff].ci[bbSOjii+1]);
	}
}
fprintf(BBSOFILE,"%d ",bbSOh);
if(bbSOflgb!=0) fprintf(BBSOFILE," %d %.1f M M ",(bbSOh-1),bbSOJ);
else{
fprintf(BBSOFILE,"%d %.1f %18.12e %18.12e ",(bbSOh-1),\
	bbSOJ,bbSO.T[0].S[bbSOff].f[bbSOjii+1],\
	bbSO.T[0].S[bbSOff].ci[bbSOjii+1]);
	} 
/* Now, Delta-J=+0 (Q),  Delta-K=-1. */ 
			bbSOflgb=0;
			if((bbSOh-1)<(bbSOLO[0].L)) bbSOflgb=1;
/* Then check if the target value for J+1 is less than or greater than the 
   allowed J's for K+1. */
			if((bbSOJ+0)>(bbSOh-1+bbSOSl)) bbSOflgb=1;
			if((bbSOJ+0)<(bbSOh-1-bbSOSl)) bbSOflgb=1;
/* Still OK?  Calculate transition */
			if(bbSOflgb==0){
				bbSOB=+0;
				bbSOell=bbSOel+0;
				bbSO.T[0].S[bbSOff].f[bbSOjii]=\
					bbSO.Rh[0].EJ[bbSOehh]-\
					bbSO.Rl[0].EJ[bbSOell]; 
				bbSOcascadetemp=bbSOtrprob*bbSOfcf*\
					bbSO.Rh[0].CJ[bbSOehh]*bbSOSATT*\
					HL(bbSOB,bbSOJ,bbSOL,bbSOdL);
/* Check for nuclear effects and modify intensity if necessary. */
		if((bbSOg!=0)&&(bbSOLO[0].L==0)){ /* low state is Sigma */
			if((((bbSOh-1)%2)==0)&&(bbSOisym==-1)){
				bbSOcascadetemp*=bbSOIa;
				}
			if((((bbSOh-1)%2)==1)&&(bbSOisym==+1)){
				bbSOcascadetemp*=bbSOIa;
				}
			} 
				bbSO.T[0].S[bbSOff].ci[bbSOjii]+=\
					bbSOcascadetemp;
				bbSO.Rl[0].CJ[bbSOell]+=bbSOfcfr*\
					bbSOcascadetemp;
				bbSOjii--;
				} 
if(bbSOdebugflag<DEBUG){
if(bbSOflgb!=0) fprintf(DBG," M M ");
else{
fprintf(DBG,"%18.12e %18.12e \n",bbSO.T[0].S[bbSOff].f[bbSOjii+1],\
	bbSO.T[0].S[bbSOff].ci[bbSOjii+1]);
	}
}
if(bbSOflgb!=0) fprintf(BBSOFILE," M M ");
else{
fprintf(BBSOFILE,"%18.12e %18.12e ",bbSO.T[0].S[bbSOff].f[bbSOjii+1],\
	bbSO.T[0].S[bbSOff].ci[bbSOjii+1]);
	} 
/* Now, Delta-J=-1 (R) for Delta-K=-1.  */ 
			bbSOflgb=0; 
			if((bbSOh-1)<(bbSOLO[0].L)) bbSOflgb=1;
			if((bbSOJ-1)>(bbSOh-1+bbSOSl)) bbSOflgb=1;
			if((bbSOJ-1)<(bbSOh-1-bbSOSl)) bbSOflgb=1;
			if((bbSOJ-1)<0) bbSOflgb=1;
			if(bbSOflgb==0){
				bbSOB=-1;
				bbSOell=bbSOel-1;
				bbSO.T[0].S[bbSOff].f[bbSOjii]=\
					bbSO.Rh[0].EJ[bbSOehh]-\
					bbSO.Rl[0].EJ[bbSOell]; 
				bbSOcascadetemp=bbSOtrprob*bbSOfcf*\
					bbSO.Rh[0].CJ[bbSOehh]*\
				       HL(bbSOB,bbSOJ,bbSOL,bbSOdL);	
/* If the high state is Sigma, the nuclear spin-stats (if homonuclear) 
   were taken care of in the rotational distribution function, so no
   need to bother with them here. */
		if((bbSOg!=0)&&(bbSOLO[0].L==0)){ /* low state is Sigma */
			if((((bbSOh-1)%2)==0)&&(bbSOisym==-1)){
				bbSOcascadetemp*=bbSOIa;
				}
			if((((bbSOh-1)%2)==1)&&(bbSOisym==+1)){
				bbSOcascadetemp*=bbSOIa;
				}
			} 
				bbSO.T[0].S[bbSOff].ci[bbSOjii]+=\
					bbSOcascadetemp;
				bbSO.Rl[0].CJ[bbSOell]+=bbSOfcfr*\
					bbSOcascadetemp; 
				}
if(bbSOdebugflag<DEBUG){
if(bbSOflgb!=0) fprintf(DBG,"M M\n");
else{
fprintf(DBG,"%18.12e %18.12e (1f)\n",bbSO.T[0].S[bbSOff].f[bbSOjii],\
	bbSO.T[0].S[bbSOff].ci[bbSOjii]);
	} 
}
if(bbSOflgb!=0) fprintf(BBSOFILE,"M M\n");
else{
fprintf(BBSOFILE,"%18.12e %18.12e\n",bbSO.T[0].S[bbSOff].f[bbSOjii],\
	bbSO.T[0].S[bbSOff].ci[bbSOjii]);
	} 
			} /* close if bbSOflga is zero. */
		} /* close upper-state J-of-K loop */ 
	} /* close upper-state K loop */ 
	fclose(BBSOFILE);
		} /* close if high-state v is cascade populated */ 

	} /* close low-state v */
	} /* close high-state v */
return; 
}

/****************  b-b_Other *****************/
  
/* This function calculates transitions between Hund's Cases (b) and (b)
   when both the upper and lower states are sigma states. */

void bb_Other(BBTinfo bbOO){
/* see parent function and header file for key to variable names */
/* indexes and dummy variables: */
int bbOOe=0,bbOOee=0,bbOOf=0,bbOOff=0,bbOOi=0,bbOOjii=0;
int bbOOeh=0,bbOOehh=0,bbOOel=0,bbOOell=0;
int bbOOh=0,bbOOk=0,bbOOfrnh=0,bbOOfrcl=0,bbOOfrch=0; 
/* variables for info about the two states */
double bbOOJ=0,bbOOSh=0,bbOOSl=0,bbOOtrprob=0,bbOOvhpop=0,bbOOhipop;
/* Intensity multiplier for satellite bands, temporary variable and fcf var */
double bbOOSATT=0,bbOOcascadetemp=0,bbOOfcf=0,bbOOfcfr=0;
/* variables for calculating Holn-London factors */
int bbOOL=0,bbOOB=0,bbOOdL=0;
/* a few flags for various purposes */
int bbOOflga=0,bbOOflgb=0,bbOOfrnhf=0,bbOOfrchf=0,bbOOLLflg=1;
/* variables for creating and writing to output files */
char bbOOfile[1000];
FILE *BBOOFILE;
Case_b_stateinfo *bbOOHI,*bbOOLO;

if(bbOOdebugflag<DEBUG){
fprintf(DBG,"bbOO.hvnlo=%d, bbOO.hvnnum=%d\n",\
		bbOO.hvnlo,bbOO.hvnnum);
fprintf(DBG,"bbOO.hvclo=%d, bbOO.hvcnum=%d\n",\
		bbOO.hvclo,bbOO.hvcnum);
fprintf(DBG,"bbOO.lvnlo=%d, bbOO.lvnnum=%d,\n",\
		bbOO.lvnlo,bbOO.lvnnum);
fprintf(DBG,"bbOO.lvclo=%d, bbOO.lvcnum=%d\n",\
		bbOO.lvclo,bbOO.lvcnum); 
fprintf(DBG,"bbOO.hivlo=%d, bbOO.hivnum=%d,\n",\
		bbOO.hivlo,bbOO.hivnum);
fprintf(DBG,"bbOO.lovlo=%d, bbOO.lovnum=%d\n",\
		bbOO.lovlo,bbOO.lovnum); 
fflush(DBG);
} 
bbOOHI=bbOO.Ch;
bbOOLO=bbOO.Cl;
bbOOSh=bbOOHI[0].S;
bbOOSl=bbOOLO[0].S;
bbOOhipop=bbOO.Sh[0].pop;
bbOOtrprob=bbOO.T[0].P[0]; 
if(bbOOHI[0].L==bbOOLO[0].L){
	bbOOLLflg=0;
	bbOOdL=0;
	} 
if((bbOOHI[0].L-bbOOLO[0].L)==+1) bbOOdL=+1;
if((bbOOHI[0].L-bbOOLO[0].L)==-1) bbOOdL=-1;
bbOOL=bbOOHI[0].L;

/* start loop down through high vib levels */
for(bbOOe=(bbOO.hivlo+bbOO.hivnum-1);bbOOe>=bbOO.hivlo;bbOOe--){ 
	bbOOee=(bbOOe-bbOO.hivlo)*bbOO.lovnum; /* simset posn. in 2D */
	bbOOfrnh=(bbOOe-bbOO.hvnlo); /* high state native vib position */
	bbOOfrch=bbOOe; /* high state cascade vib position */
/* these flags (bbOOfrnhf and bbOOfrchf) tell if this high vib state is 
   populated natively, by cascade, or both.  They will be used later */ 
	if((bbOOe>=bbOO.hvnlo)&&(bbOOe<(bbOO.hvnlo+bbOO.hvnnum))){
		bbOOfrnhf=0;
		bbOOvhpop=bbOO.Vh[0].p[bbOOe-bbOO.hvnlo];
		}
	else bbOOfrnhf=-1; 
	if((bbOOe>=bbOO.hvclo)&&(bbOOe<(bbOO.hvclo+bbOO.hvcnum))){
		bbOOfrchf=0;
		}
	else bbOOfrchf=-1; 
if(bbOOdebugflag<DEBUG){
fprintf(DBG,"bbOOfrnhf=%d, bbOOvhpop=%f, bbOOfrchf=%d, ",\
		bbOOfrnhf,bbOOvhpop,bbOOfrchf);
fprintf(DBG,"bbOOe=%d, bbOOee=%d\nbbOOfrnh=%d, bbOOfrch=%d, ",\
		bbOOe,bbOOee,bbOOfrnh,bbOOfrch);
fflush(DBG);
} 
/* start loop down through low vib levels */
for(bbOOf=(bbOO.lovlo+bbOO.lovnum-1);bbOOf>=bbOO.lovlo;bbOOf--){ 
	bbOOfcf=bbOO.T[0].v[bbOOe].fcfn[bbOOf];
	if(bbOOfcf!=0){bbOOfcfr=bbOO.T[0].v[bbOOe].fcfc[bbOOf]/bbOOfcf;}
	bbOOff=bbOOee+bbOOf-bbOO.lovlo; /* position in last dimension */
	bbOOfrcl=bbOOf;
if(bbOOdebugflag<DEBUG){
fprintf(DBG,"bbOOf=%d, bbOO.lovlo=%d, bbOO.lovnum=%d \n",\
		bbOOf,bbOO.lovlo,bbOO.lovnum);
fflush(DBG);
}
bbOOff=bbOOee+bbOOf-bbOO.lovlo; /* simset posn. in final dimension */
	bbOO.T[0].S[bbOOff].n=bbOO.nT*bbOO.hiK;
	bbOO.T[0].S[bbOOff].f=\
		(double*)calloc(bbOO.T[0].S[bbOOff].n,sizeof(double));
	bbOO.T[0].S[bbOOff].ni=\
		(double*)calloc(bbOO.T[0].S[bbOOff].n,sizeof(double));
	bbOO.T[0].S[bbOOff].ci=\
		(double*)calloc(bbOO.T[0].S[bbOOff].n,sizeof(double)); 
bbOO.Rl=&bbOO.Sl[0].rc[bbOOfrcl]; 
if(bbOOdebugflag<DEBUG){
fprintf(DBG,"bbOOff=%d, bbOOfrcl=%d, bbOO.T[0].S[bbOOff].n=%d\n",\
		bbOOff,bbOOfrcl,bbOO.T[0].S[bbOOff].n);
fflush(DBG);
}
/* start with highest J value (same reason as before), and loop down looking 
   for lower J's to which to transit.  assign intensities.  If both states are
   Omega=0, then assign zero intensity to transitions for J=0<->J=0.  */ 
	if(bbOOfrnhf==0){/* if this v is natively populated.*/
sprintf(bbOOfile,"%s_molecules/%s/%s--%s_v%d--v%d_NAT.dat",\
	PREF,MOL[bbOO.m].Mol,bbOO.T[0].Nhi,bbOO.T[0].Nlo,bbOOe,bbOOf);
		BBOOFILE=fopen(bbOOfile,"w");
		if(BBOOFILE==NULL){
printf("Error opening transition sub-file %s.  Exiting.\n",bbOOfile);
			exit(1);
			}
fprintf(BBOOFILE,"# File created by %s.\n# File contains ",PROGRAM_NAME);
fprintf(BBOOFILE,"P, Q and R");
fprintf(BBOOFILE," branches from the simulation titled %s\n",PREF);
fprintf(BBOOFILE,"# THESE INTENSITIES ARE DUE TO USER-SPECIFIED POPULATIONS");
fprintf(BBOOFILE," ONLY -- NO CASCADE\n");
fprintf(BBOOFILE,"# Transition is Hund's Case (b) to Hund's Case (b).\n");
fprintf(BBOOFILE,"# Molecule %s, \tHigh state: %s, v=%d\n",\
		MOL[bbOO.m].Mol,bbOO.T[0].Nhi,bbOOe);
fprintf(BBOOFILE,"#\t\tLow state: %s, v=%d\n# Columns are: ",\
		bbOO.T[0].Nlo,bbOOf);
fprintf(BBOOFILE,"K_hi K_lo J_hi P(J_hi+1)_v P_i ");
fprintf(BBOOFILE,"Q(J_hi+0)_v Q_i R(J_hi-1)_v R_i\n#\n");

		bbOO.Rh=&bbOO.Sh[0].r[bbOOfrnh]; 
		if((bbOO.hiK)<bbOO.Sh[0].r[bbOOfrnh].k){
			bbOOk=bbOO.hiK-1;
			}
		else bbOOk=bbOO.Sh[0].r[bbOOfrnh].k-1; 
		if(bbOO.Rl[0].Kdissoc<(bbOOk-1)){ 
fprintf(PAR,"\nWARNING!! State %s of molecule %s has significant ",\
		bbOO.T[0].Nhi,MOL[bbOO.m].Mol);
fprintf(PAR,"rotational population at level %d\n",bbOOk);
fprintf(PAR,"\tBut lower state %s dissociates at level %.1f\n",\
		bbOO.T[0].Nlo,bbOO.Rl[0].Kdissoc);
fprintf(PAR,"It is likely that the bound upper state is transiting to (an)");
fprintf(PAR," unbound lower state(s).\n");
fprintf(PAR,"Look for continuum emissions trailing off to the low-energy ");
fprintf(PAR,"end of the v(hi)=%d to v(lo)=%d\nband in the real ",bbOOe,bbOOf);
fprintf(PAR,"spectrum (not the simulated one)\nCurrently, this program will");
fprintf(PAR," not simulate bound-->unbound spectra.\n\n");
			bbOOk=(int)(bbOO.Rl[0].Kdissoc);
			}
if(bbOOdebugflag<DEBUG){
fprintf(DBG,"bbOOfrnh=%d, bbOOk=%d, ",bbOOfrnh,bbOOk);
fprintf(DBG,"bbOO.Sh[0].r[%d].k=%d\n",bbOOfrnh,bbOO.Sh[0].r[bbOOfrnh].k);
fflush(DBG);
}
/* To make the output file easy to write, I'm cycling by high-state K first,
   then by high-state J.  From that J/K(hi) pair, the program will look
   for a P transition and an R transition to K+1 then to K-1.  If all the
   necessary states exist, it will check to see if there are any symmetry
   or nuclear spin effects to be concerned about.  Taking these into account,
   it will calculate the transition. */
/* loop first by upper-state K-value */
for(bbOOh=bbOOk;bbOOh>=0;bbOOh--){
	bbOOjii=(bbOOh+1)*bbOO.nT-1; /* position in simset array */
	bbOOeh=bbOOh*bbOO.nhJ; /* position in 2D for high EJ & PJ */ 
/* then loop by upper-state J for this K */
	for(bbOOi=(bbOO.nhJ-1);bbOOi>=0;bbOOi--){ 
		bbOOJ=bbOOh+bbOOi-bbOOSh;
		bbOOehh=bbOOeh+bbOOi;
		bbOOflga=0;
		if(bbOOh<bbOOHI[0].L) bbOOflga=1; /* if this K is less than 
			Lambda for the high state, don't proceed */
		if(bbOOJ<0) bbOOflga=1; /* if J is less than zero, don't */ 
if(bbOOdebugflag<DEBUG){
fprintf(DBG,"bbOOh=%d, bbOOi=%d, bbOOSh=%.1f, bbOOHI[0].L=%d, bbOOflga=%d\n",\
		bbOOh,bbOOi,bbOOSh,bbOOHI[0].L,bbOOflga);
fflush(DBG);
} 
		if(bbOOflga==0){ 
/* Check Delta-J=+1 (P), then Delta-J=0 (Q), and then Delta-J=-1 (R) for 
Delta-K=+1. First, check to see if the target K exists: */ 
			bbOOflgb=0;
			if((bbOOh+1)<(bbOOLO[0].L)) bbOOflgb=1;
/* Relative position in the low state E/P array */
			bbOOel=(bbOOh+1)*bbOO.nlJ + bbOO.nlJ-bbOO.nhJ + \
				(int)(bbOOSh-bbOOSl) -1 + bbOOi; 
if(bbOOdebugflag<DEBUG){
fprintf(DBG,"\t(bbOOh+1)*bbOO.nlJ=%d, bbOO.nlJ-bbOO.nhJ=%d\n",\
		(bbOOh+1)*bbOO.nlJ,bbOO.nlJ-bbOO.nhJ);
fprintf(DBG,"\t(int)(bbOOSh-bbOOSl)=%d, bbOOi=%d, bbOOel=%d\n",\
		(int)(bbOOSh-bbOOSl),bbOOi,bbOOel);
fflush(DBG);
} 
			if((bbOOel+1)>=(bbOO.Rl[0].kc*bbOO.nlJ)) bbOOflgb=1;
/* Then check if the target value for J+1 is less than or greater than the 
   allowed J's for K+1. */
			if((bbOOJ+1)>(bbOOh+1+bbOOSl)) bbOOflgb=1;
			if((bbOOJ+1)<(bbOOh+1-bbOOSl)) bbOOflgb=1;
/* Still OK?  Calculate transition */
			if(bbOOflgb==0){
				bbOOB=+1;
				bbOOell=bbOOel+1;
				bbOO.T[0].S[bbOOff].f[bbOOjii]=\
					bbOO.Rh[0].EJ[bbOOehh]-\
					bbOO.Rl[0].EJ[bbOOell]; 
				bbOO.T[0].S[bbOOff].ni[bbOOjii]=\
					bbOOhipop*bbOOtrprob*bbOOvhpop*\
					bbOO.Rh[0].PJ[bbOOehh]*bbOOfcf*\
					HL(bbOOB,bbOOJ,bbOOL,bbOOdL); 
/* set low-state population for cascade */
				bbOO.Rl[0].CJ[bbOOell]+=bbOOfcfr*\
					bbOO.T[0].S[bbOOff].ni[bbOOjii];
				bbOOjii--;
				} 
if(bbOOdebugflag<DEBUG){
fprintf(DBG,"%d ",bbOOh);
if(bbOOflgb!=0) fprintf(DBG," %d %.1f M M ",(bbOOh+1),bbOOJ);
else{
fprintf(DBG,"%d %.1f %18.12e %18.12e ",(bbOOh+1),\
	bbOOJ,bbOO.T[0].S[bbOOff].f[bbOOjii+1],\
	bbOO.T[0].S[bbOOff].ni[bbOOjii+1]);
	} 
}
fprintf(BBOOFILE,"%d ",bbOOh);
if(bbOOflgb!=0) fprintf(BBOOFILE," %d %.1f M M ",(bbOOh+1),bbOOJ);
else{
fprintf(BBOOFILE,"%d %.1f %18.12e %18.12e ",(bbOOh+1),\
	bbOOJ,bbOO.T[0].S[bbOOff].f[bbOOjii+1],\
	bbOO.T[0].S[bbOOff].ni[bbOOjii+1]); }

/* Now Delta-J=0 (Q) (for Delta-K=+1)  */ 
			bbOOflgb=0;
			if((bbOOh+1)<(bbOOLO[0].L)) bbOOflgb=1;
/* Then check if the target value for J+1 is less than or greater than the 
   allowed J's for K+1. */
			if((bbOOJ+0)>(bbOOh+1+bbOOSl)) bbOOflgb=1;
			if((bbOOJ+0)<(bbOOh+1-bbOOSl)) bbOOflgb=1;
/* Still OK?  Calculate transition */
			if(bbOOflgb==0){
				bbOOB=0;
				bbOOell=bbOOel+0;
				bbOO.T[0].S[bbOOff].f[bbOOjii]=\
					bbOO.Rh[0].EJ[bbOOehh]-\
					bbOO.Rl[0].EJ[bbOOell]; 
				bbOO.T[0].S[bbOOff].ni[bbOOjii]=\
					bbOOhipop*bbOOtrprob*bbOOvhpop*\
					bbOO.Rh[0].PJ[bbOOehh]*bbOOSATT*\
					HL(bbOOB,bbOOJ,bbOOL,bbOOdL)*bbOOfcf;
/* set low-state population for cascade */
				bbOO.Rl[0].CJ[bbOOell]+=bbOOfcfr*\
					bbOO.T[0].S[bbOOff].ni[bbOOjii];
				bbOOjii--;
				} 
if(bbOOdebugflag<DEBUG){
fprintf(DBG,"%d ",bbOOh);
if(bbOOflgb!=0) fprintf(DBG," %d %.1f M M ",(bbOOh+1),bbOOJ);
else{
fprintf(DBG,"%d %.1f %18.12e %18.12e ",(bbOOh+1),\
	bbOOJ,bbOO.T[0].S[bbOOff].f[bbOOjii+1],\
	bbOO.T[0].S[bbOOff].ni[bbOOjii+1]);
	} 
}
if(bbOOflgb!=0) fprintf(BBOOFILE," M M ");
else{
fprintf(BBOOFILE,"%18.12e %18.12e ",bbOO.T[0].S[bbOOff].f[bbOOjii+1],\
	bbOO.T[0].S[bbOOff].ni[bbOOjii+1]);
	} 
/* Now, Delta-J=-1 (R) for Delta-K=+1.  */ 
			bbOOflgb=0; 
			if((bbOOh+1)<(bbOOLO[0].L)) bbOOflgb=1;
			if((bbOOJ-1)>(bbOOh+1+bbOOSl)) bbOOflgb=1;
			if((bbOOJ-1)<(bbOOh+1-bbOOSl)) bbOOflgb=1;
			if((bbOOJ-1)<0) bbOOflgb=1;
			if(bbOOflgb==0){
				bbOOB=-1;
				bbOOell=bbOOel-1;
				bbOO.T[0].S[bbOOff].f[bbOOjii]=\
					bbOO.Rh[0].EJ[bbOOehh]-\
					bbOO.Rl[0].EJ[bbOOell]; 
				bbOO.T[0].S[bbOOff].ni[bbOOjii]=\
					bbOOhipop*bbOOtrprob*bbOOvhpop*\
					bbOO.Rh[0].PJ[bbOOehh]*bbOOSATT*\
					HL(bbOOB,bbOOJ,bbOOL,bbOOdL)*bbOOfcf;
				bbOO.Rl[0].CJ[bbOOell]+=bbOOfcfr*\
					bbOO.T[0].S[bbOOff].ni[bbOOjii];
				bbOOjii--;
				}
if(bbOOdebugflag<DEBUG){
if(bbOOflgb!=0) fprintf(DBG,"M M\n");
else{
fprintf(DBG,"%18.12e %18.12e\n",bbOO.T[0].S[bbOOff].f[bbOOjii+1],\
	bbOO.T[0].S[bbOOff].ni[bbOOjii+1]);
	} 
}
if(bbOOflgb!=0) fprintf(BBOOFILE,"M M\n");
else{
fprintf(BBOOFILE,"%18.12e %18.12e\n",bbOO.T[0].S[bbOOff].f[bbOOjii+1],\
	bbOO.T[0].S[bbOOff].ni[bbOOjii+1]);
	} 
/* Check Delta-J=+1 (P), then Delta-J=0 (Q), and then Delta-J=-1 (R) for 
Delta-K=+0. First, check to see if the target K exists: */ 
			bbOOflgb=0;
			if((bbOOh+0)<(bbOOLO[0].L)) bbOOflgb=1;
/* Relative position in the low state E/P array */
			bbOOel=(bbOOh+0)*bbOO.nlJ + bbOO.nlJ-bbOO.nhJ + \
				(int)(bbOOSh-bbOOSl) + 0 + bbOOi; 
if(bbOOdebugflag<DEBUG){
fprintf(DBG,"\t(bbOOh+0)*bbOO.nlJ=%d, bbOO.nlJ-bbOO.nhJ=%d\n",\
		(bbOOh+0)*bbOO.nlJ,bbOO.nlJ-bbOO.nhJ);
fprintf(DBG,"\t(int)(bbOOSh-bbOOSl)=%d, bbOOi=%d, bbOOel=%d\n",\
		(int)(bbOOSh-bbOOSl),bbOOi,bbOOel);
fflush(DBG);
} 
			if((bbOOel+1)>=(bbOO.Rl[0].kc*bbOO.nlJ)) bbOOflgb=1;
/* Then check if the target value for J+1 is less than or greater than the 
   allowed J's for K+0. */
			if((bbOOJ+1)>(bbOOh+0+bbOOSl)) bbOOflgb=1;
			if((bbOOJ+1)<(bbOOh+0-bbOOSl)) bbOOflgb=1;
/* Still OK?  Calculate transition */
			if(bbOOflgb==0){
				bbOOB=+1;
				bbOOell=bbOOel+1;
				bbOO.T[0].S[bbOOff].f[bbOOjii]=\
					bbOO.Rh[0].EJ[bbOOehh]-\
					bbOO.Rl[0].EJ[bbOOell]; 
				bbOO.T[0].S[bbOOff].ni[bbOOjii]=\
					bbOOhipop*bbOOtrprob*bbOOvhpop*\
					bbOO.Rh[0].PJ[bbOOehh]*bbOOSATT*\
					HL(bbOOB,bbOOJ,bbOOL,bbOOdL)*bbOOfcf;
/* set low-state population for cascade */
				bbOO.Rl[0].CJ[bbOOell]+=bbOOfcfr*\
					bbOO.T[0].S[bbOOff].ni[bbOOjii];
				bbOOjii--;
				} 
if(bbOOdebugflag<DEBUG){
fprintf(DBG,"%d ",bbOOh);
if(bbOOflgb!=0) fprintf(DBG," %d %.1f M M ",(bbOOh+0),bbOOJ);
else{
fprintf(DBG,"%d %.1f %18.12e %18.12e ",(bbOOh+0),\
	bbOOJ,bbOO.T[0].S[bbOOff].f[bbOOjii+1],\
	bbOO.T[0].S[bbOOff].ni[bbOOjii+1]);
	} 
}
fprintf(BBOOFILE,"%d ",bbOOh);
if(bbOOflgb!=0) fprintf(BBOOFILE," %d %.1f M M ",(bbOOh+0),bbOOJ);
else{
fprintf(BBOOFILE,"%d %.1f %18.12e %18.12e ",(bbOOh+0),\
	bbOOJ,bbOO.T[0].S[bbOOff].f[bbOOjii+1],\
	bbOO.T[0].S[bbOOff].ni[bbOOjii+1]);
	} 
/* Now Delta-J=0 (Q) (for Delta-K=+0)  */ 
			bbOOflgb=0;
			if((bbOOh+0)<(bbOOLO[0].L)) bbOOflgb=1;
/* Then check if the target value for J+1 is less than or greater than the 
   allowed J's for K+0. */
			if((bbOOJ+0)>(bbOOh+0+bbOOSl)) bbOOflgb=1;
			if((bbOOJ+0)<(bbOOh+0-bbOOSl)) bbOOflgb=1;
/* Still OK?  Calculate transition */
			if(bbOOflgb==0){
				bbOOB=0;
				bbOOell=bbOOel+0;
				bbOO.T[0].S[bbOOff].f[bbOOjii]=\
					bbOO.Rh[0].EJ[bbOOehh]-\
					bbOO.Rl[0].EJ[bbOOell]; 
				bbOO.T[0].S[bbOOff].ni[bbOOjii]=\
					bbOOhipop*bbOOtrprob*bbOOvhpop*\
					bbOO.Rh[0].PJ[bbOOehh]*bbOOfcf*\
					HL(bbOOB,bbOOJ,bbOOL,bbOOdL);
/* set low-state population for cascade */
				bbOO.Rl[0].CJ[bbOOell]+=bbOOfcfr*\
					bbOO.T[0].S[bbOOff].ni[bbOOjii];
				bbOOjii--;
				} 
if(bbOOdebugflag<DEBUG){
fprintf(DBG,"%d ",bbOOh);
if(bbOOflgb!=0) fprintf(DBG," M M ");
else{
fprintf(DBG," %18.12e %18.12e ",bbOO.T[0].S[bbOOff].f[bbOOjii+1],\
	bbOO.T[0].S[bbOOff].ni[bbOOjii+1]);
	} 
}
if(bbOOflgb!=0) fprintf(BBOOFILE," M M ");
else{
fprintf(BBOOFILE,"%18.12e %18.12e ",bbOO.T[0].S[bbOOff].f[bbOOjii+1],\
	bbOO.T[0].S[bbOOff].ni[bbOOjii+1]);
	} 
/* Now, Delta-J=-1 (R) for Delta-K=+0.  */ 
			bbOOflgb=0; 
			if((bbOOh+0)<(bbOOLO[0].L)) bbOOflgb=1;
			if((bbOOJ-1)>(bbOOh+0+bbOOSl)) bbOOflgb=1;
			if((bbOOJ-1)<(bbOOh+0-bbOOSl)) bbOOflgb=1;
			if((bbOOJ-1)<0) bbOOflgb=1;
			if(bbOOflgb==0){
				bbOOB=-1;
				bbOOell=bbOOel-1;
				bbOO.T[0].S[bbOOff].f[bbOOjii]=\
					bbOO.Rh[0].EJ[bbOOehh]-\
					bbOO.Rl[0].EJ[bbOOell]; 
				bbOO.T[0].S[bbOOff].ni[bbOOjii]=\
					bbOOhipop*bbOOtrprob*bbOOvhpop*\
					bbOO.Rh[0].PJ[bbOOehh]*bbOOSATT*\
					HL(bbOOB,bbOOJ,bbOOL,bbOOdL)*bbOOfcf;
				bbOO.Rl[0].CJ[bbOOell]+=bbOOfcfr*\
					bbOO.T[0].S[bbOOff].ni[bbOOjii];
				bbOOjii--;
				}
if(bbOOdebugflag<DEBUG){
if(bbOOflgb!=0) fprintf(DBG,"M M\n");
else{
fprintf(DBG,"%18.12e %18.12e\n",bbOO.T[0].S[bbOOff].f[bbOOjii+1],\
	bbOO.T[0].S[bbOOff].ni[bbOOjii+1]);
	} 
}
if(bbOOflgb!=0) fprintf(BBOOFILE,"M M\n");
else{
fprintf(BBOOFILE,"%18.12e %18.12e\n",bbOO.T[0].S[bbOOff].f[bbOOjii+1],\
	bbOO.T[0].S[bbOOff].ni[bbOOjii+1]);
	} 
/* Check Delta-J=+1 (P), then Delta-J=0 and Delta-J=-1 (R) for Delta-K=-1. 
   First, check to see if the target K exists: */ 
			bbOOflgb=0;
			if((bbOOh-1)<(bbOOLO[0].L)) bbOOflgb=1;
/* Relative position in the low state E/P array */
			bbOOel=(bbOOh-1)*bbOO.nlJ + bbOO.nlJ-bbOO.nhJ + \
				(int)(bbOOSh-bbOOSl) + 1 + bbOOi; 
if(bbOOdebugflag<DEBUG){
fprintf(DBG,"\t(bbOOh-1)*bbOO.nlJ=%d, bbOO.nlJ-bbOO.nhJ=%d\n",\
		(bbOOh-1)*bbOO.nlJ,bbOO.nlJ-bbOO.nhJ);
fprintf(DBG,"\t(int)(bbOOSh-bbOOSl)=%d, bbOOi=%d, bbOOel=%d\n",\
		(int)(bbOOSh-bbOOSl),bbOOi,bbOOel);
fflush(DBG);
} 
			if((bbOOel+1)>=(bbOO.Rl[0].kc*bbOO.nlJ)) bbOOflgb=1;
/* Then check if the target value for J+1 is less than or greater than the 
   allowed J's for K-1. */
			if((bbOOJ+1)>(bbOOh-1+bbOOSl)) bbOOflgb=1;
			if((bbOOJ+1)<(bbOOh-1-bbOOSl)) bbOOflgb=1;
/* Still OK?  Calculate transition */
			if(bbOOflgb==0){
				bbOOB=+1;
				bbOOell=bbOOel+1;
				bbOO.T[0].S[bbOOff].f[bbOOjii]=\
					bbOO.Rh[0].EJ[bbOOehh]-\
					bbOO.Rl[0].EJ[bbOOell]; 
				bbOO.T[0].S[bbOOff].ni[bbOOjii]=\
					bbOOhipop*bbOOtrprob*bbOOvhpop*\
					bbOO.Rh[0].PJ[bbOOehh]*bbOOSATT*\
					HL(bbOOB,bbOOJ,bbOOL,bbOOdL)*bbOOfcf;
				bbOO.Rl[0].CJ[bbOOell]+=bbOOfcfr*\
					bbOO.T[0].S[bbOOff].ni[bbOOjii];
				bbOOjii--;
				} 
if(bbOOdebugflag<DEBUG){
fprintf(DBG,"%d ",bbOOh);
if(bbOOflgb!=0) fprintf(DBG," %d %.1f M M ",(bbOOh-1),bbOOJ);
else{
fprintf(DBG,"%d %.1f %18.12e %18.12e ",(bbOOh-1),\
	bbOOJ,bbOO.T[0].S[bbOOff].f[bbOOjii+1],\
	bbOO.T[0].S[bbOOff].ni[bbOOjii+1]);
	}
}
fprintf(BBOOFILE,"%d ",bbOOh);
if(bbOOflgb!=0) fprintf(BBOOFILE," %d %.1f M M ",(bbOOh-1),bbOOJ);
else{
fprintf(BBOOFILE,"%d %.1f %18.12e %18.12e ",(bbOOh-1),\
	bbOOJ,bbOO.T[0].S[bbOOff].f[bbOOjii+1],\
	bbOO.T[0].S[bbOOff].ni[bbOOjii+1]);
	} 
/* Now, Delta-J=+0 (Q) for Delta-K=-1. */ 
			bbOOflgb=0;
			if((bbOOh-1)<(bbOOLO[0].L)) bbOOflgb=1;
/* Then check if the target value for J+1 is less than or greater than the 
   allowed J's for K-1. */
			if((bbOOJ+0)>(bbOOh-1+bbOOSl)) bbOOflgb=1;
			if((bbOOJ+0)<(bbOOh-1-bbOOSl)) bbOOflgb=1;
/* Still OK?  Calculate transition */
			if(bbOOflgb==0){
				bbOOB=0;
				bbOOell=bbOOel+0;
				bbOO.T[0].S[bbOOff].f[bbOOjii]=\
					bbOO.Rh[0].EJ[bbOOehh]-\
					bbOO.Rl[0].EJ[bbOOell]; 
				bbOO.T[0].S[bbOOff].ni[bbOOjii]=\
					bbOOhipop*bbOOtrprob*bbOOvhpop*\
					bbOO.Rh[0].PJ[bbOOehh]*bbOOSATT*\
					HL(bbOOB,bbOOJ,bbOOL,bbOOdL)*bbOOfcf;
				bbOO.Rl[0].CJ[bbOOell]+=bbOOfcfr*\
					bbOO.T[0].S[bbOOff].ni[bbOOjii];
				bbOOjii--;
				} 
if(bbOOdebugflag<DEBUG){
fprintf(DBG,"%d ",bbOOh);
if(bbOOflgb!=0) fprintf(DBG," %d %.1f M M ",(bbOOh-1),bbOOJ);
else{
fprintf(DBG,"%d %.1f %18.12e %18.12e ",(bbOOh-1),\
	bbOOJ,bbOO.T[0].S[bbOOff].f[bbOOjii+1],\
	bbOO.T[0].S[bbOOff].ni[bbOOjii+1]);
	}
}
if(bbOOflgb!=0) fprintf(BBOOFILE," M M ");
else{
fprintf(BBOOFILE," %18.12e %18.12e ",bbOO.T[0].S[bbOOff].f[bbOOjii+1],\
	bbOO.T[0].S[bbOOff].ni[bbOOjii+1]);
	} 
/* Now, Delta-J=-1 (R) for Delta-K=-1.  */ 
			bbOOflgb=0; 
			if((bbOOh-1)<(bbOOLO[0].L)) bbOOflgb=1;
			if((bbOOJ-1)>(bbOOh-1+bbOOSl)) bbOOflgb=1;
			if((bbOOJ-1)<(bbOOh-1-bbOOSl)) bbOOflgb=1;
			if((bbOOJ-1)<0) bbOOflgb=1;
			if(bbOOflgb==0){
				bbOOB=-1;
				bbOOell=bbOOel-1;
				bbOO.T[0].S[bbOOff].f[bbOOjii]=\
					bbOO.Rh[0].EJ[bbOOehh]-\
					bbOO.Rl[0].EJ[bbOOell]; 
				bbOO.T[0].S[bbOOff].ni[bbOOjii]=\
					bbOOhipop*bbOOtrprob*bbOOvhpop*\
					bbOO.Rh[0].PJ[bbOOehh]*bbOOfcf*\
				       HL(bbOOB,bbOOJ,bbOOL,bbOOdL);	
				bbOO.Rl[0].CJ[bbOOell]+=bbOOfcfr*\
					bbOO.T[0].S[bbOOff].ni[bbOOjii]; 
				}
if(bbOOdebugflag<DEBUG){
if(bbOOflgb!=0) fprintf(DBG,"M M\n");
else{
fprintf(DBG,"%18.12e %18.12e\n",bbOO.T[0].S[bbOOff].f[bbOOjii],\
	bbOO.T[0].S[bbOOff].ni[bbOOjii]);
	} 
}
if(bbOOflgb!=0) fprintf(BBOOFILE,"M M\n");
else{
fprintf(BBOOFILE,"%18.12e %18.12e\n",bbOO.T[0].S[bbOOff].f[bbOOjii],\
	bbOO.T[0].S[bbOOff].ni[bbOOjii]);
	} 
			} /* close if bbOOflga is zero. */
		} /* close upper-state J-of-K loop */ 
	} /* close upper-state K loop */ 
	fclose(BBOOFILE);
		} /* close if high-state v is natively populated */ 

	if(bbOOfrchf==0){/* if this v is cascade populated.*/
sprintf(bbOOfile,"%s_molecules/%s/%s--%s_v%d--v%d_CAS.dat",\
	PREF,MOL[bbOO.m].Mol,bbOO.T[0].Nhi,bbOO.T[0].Nlo,bbOOe,bbOOf);
		BBOOFILE=fopen(bbOOfile,"w");
		if(BBOOFILE==NULL){
printf("Error opening transition sub-file %s.  Exiting.\n",bbOOfile);
			exit(1);
			}
fprintf(BBOOFILE,"# File created by %s.\n# File contains ",PROGRAM_NAME);
fprintf(BBOOFILE,"P and R");
fprintf(BBOOFILE," branches from the simulation titled %s\n",PREF);
fprintf(BBOOFILE,"# THESE INTENSITIES ARE DUE TO CASCADE ONLY -- ");
fprintf(BBOOFILE,"NO USER-DEFINED POPULATIONS\n");
fprintf(BBOOFILE,"# Transition is Hund's Case (b) to Hund's Case (b).\n");
fprintf(BBOOFILE,"# Molecule %s, \tHigh state: %s, v=%d\n",\
		MOL[bbOO.m].Mol,bbOO.T[0].Nhi,bbOOe);
fprintf(BBOOFILE,"#\t\tLow state: %s, v=%d\n# Columns are: ",\
		bbOO.T[0].Nlo,bbOOf);
fprintf(BBOOFILE,"K_hi K_lo J_hi P(J_hi+1)_v P_i  R(J_hi-1)_v R_i\n#\n");

		bbOO.Rh=&bbOO.Sh[0].rc[bbOOfrch]; 
		if((bbOO.hiK)<bbOO.Sh[0].rc[bbOOfrch].kc){
			bbOOk=bbOO.hiK-1;
			}
		else bbOOk=bbOO.Sh[0].rc[bbOOfrch].kc-1; 
		if(bbOO.Rl[0].Kdissoc<(bbOOk-1)){ 
fprintf(PAR,"\nWARNING!! State %s of molecule %s has significant ",\
		bbOO.T[0].Nhi,MOL[bbOO.m].Mol);
fprintf(PAR,"CASCADE rotational population at level %d\n",bbOOk);
fprintf(PAR,"\tBut lower state %s dissociates at level %.1f\n",\
		bbOO.T[0].Nlo,bbOO.Rl[0].Kdissoc);
fprintf(PAR,"It is likely that the bound upper state is transiting to (an)");
fprintf(PAR," unbound lower state(s).\n");
fprintf(PAR,"Look for continuum emissions trailing off to the low-energy ");
fprintf(PAR,"end of the v(hi)=%d to v(lo)=%d\nband in the real ",bbOOe,bbOOf);
fprintf(PAR,"spectrum (not the simulated one)\nCurrently, this program will");
fprintf(PAR," not simulate bound-->unbound spectra.\n\n");
			bbOOk=(int)(bbOO.Rl[0].Kdissoc);
			}
if(bbOOdebugflag<DEBUG){
fprintf(DBG,"bbOOfrch=%d, bbOOk=%d, ",bbOOfrch,bbOOk);
fprintf(DBG,"bbOO.Sh[0].rc[%d].kc=%d\n",bbOOfrch,bbOO.Sh[0].rc[bbOOfrch].kc);
fflush(DBG);
}
/* To make the output file easy to write, I'm cycling by high-state K first,
   then by high-state J.  From that J/K(hi) pair, the program will look
   for a P transition and an R transition to K+1 then to K-1.  If all the
   necessary states exist, it will check to see if there are any symmetry
   or nuclear spin effects to be concerned about.  Taking these into account,
   it will calculate the transition. */
/* loop first by upper-state K-value */
for(bbOOh=(bbOOk-1);bbOOh>=0;bbOOh--){ 
	bbOOjii=(bbOOh+1)*bbOO.nT-1; /* position in simset array */
	bbOOeh=bbOOh*bbOO.nhJ; /* position in 2D for high EJ & PJ */ 
/* then loop by upper-state J for this K */
	for(bbOOi=(bbOO.nhJ-1);bbOOi>=0;bbOOi--){ 
		bbOOJ=bbOOh+bbOOi-bbOOSh;
		bbOOehh=bbOOeh+bbOOi;
		bbOOflga=0;
		if(bbOOh<bbOOHI[0].L) bbOOflga=1; /* if this K is less than 
			Lambda for the high state, don't proceed */
		if(bbOOJ<0) bbOOflga=1; /* if J is less than zero, don't */ 
if(bbOOdebugflag<DEBUG){
fprintf(DBG,"bbOOh=%d, bbOOi=%d, bbOOSh=%.1f, bbOOHI[0].L=%d, bbOOflga=%d\n",\
		bbOOh,bbOOi,bbOOSh,bbOOHI[0].L,bbOOflga);
fflush(DBG);
} 
		if(bbOOflga==0){ 
/* Check Delta-J=+1 (P), then Delta-J=0 (Q) and Delta-J=-1 (R) for Delta-K=+1.
   First, check to see if the target K exists: */ 
			bbOOflgb=0;
			if((bbOOh+1)<(bbOOLO[0].L)) bbOOflgb=1;
/* Relative position in the low state E/P array */
			bbOOel=(bbOOh+1)*bbOO.nlJ + bbOO.nlJ-bbOO.nhJ + \
				(int)(bbOOSh-bbOOSl) -1 + bbOOi; 
if(bbOOdebugflag<DEBUG){
fprintf(DBG,"\t(bbOOh+1)*bbOO.nlJ=%d, bbOO.nlJ-bbOO.nhJ=%d\n",\
		(bbOOh+1)*bbOO.nlJ,bbOO.nlJ-bbOO.nhJ);
fprintf(DBG,"\t(int)(bbOOSh-bbOOSl)=%d, bbOOi=%d, bbOOel=%d\n",\
		(int)(bbOOSh-bbOOSl),bbOOi,bbOOel);
fflush(DBG);
} 
			if((bbOOel+1)>=(bbOO.Rl[0].kc*bbOO.nlJ)) bbOOflgb=1;
/* Then check if the target value for J+1 is less than or greater than the 
   allowed J's for K+1. */
			if((bbOOJ+1)>(bbOOh+1+bbOOSl)) bbOOflgb=1;
			if((bbOOJ+1)<(bbOOh+1-bbOOSl)) bbOOflgb=1;
/* Still OK?  Calculate transition */
			if(bbOOflgb==0){
				bbOOB=+1;
				bbOOell=bbOOel+1;
				bbOO.T[0].S[bbOOff].f[bbOOjii]=\
					bbOO.Rh[0].EJ[bbOOehh]-\
					bbOO.Rl[0].EJ[bbOOell]; 
				bbOOcascadetemp=bbOOtrprob*bbOOfcf*\
					bbOO.Rh[0].CJ[bbOOehh]*\
					HL(bbOOB,bbOOJ,bbOOL,bbOOdL);
				bbOO.T[0].S[bbOOff].ci[bbOOjii]+=\
					bbOOcascadetemp;
				bbOO.Rl[0].CJ[bbOOell]+=bbOOfcfr*\
					bbOOcascadetemp;
				bbOOjii--;
				} 
if(bbOOdebugflag<DEBUG){
fprintf(DBG,"%d ",bbOOh);
if(bbOOflgb!=0) fprintf(DBG," %d %.1f (1c) M M\n",(bbOOh+1),bbOOJ);
else{
fprintf(DBG,"%d %.1f %18.12e %18.12e (1c)\n",(bbOOh+1),\
	bbOOJ,bbOO.T[0].S[bbOOff].f[bbOOjii+1],\
	bbOO.T[0].S[bbOOff].ci[bbOOjii+1]);
	} 
}
fprintf(BBOOFILE,"%d ",bbOOh);
if(bbOOflgb!=0) fprintf(BBOOFILE," %d %.1f M M ",(bbOOh+1),bbOOJ);
else{
fprintf(BBOOFILE,"%d %.1f %18.12e %18.12e ",(bbOOh+1),\
	bbOOJ,bbOO.T[0].S[bbOOff].f[bbOOjii+1],\
	bbOO.T[0].S[bbOOff].ci[bbOOjii+1]);
	} 
/* Now, Delta-J=0 (Q), Delta-K=+1 */ 
			bbOOflgb=0;
			if((bbOOh+1)<(bbOOLO[0].L)) bbOOflgb=1;
/* Then check if the target value for J+1 is less than or greater than the 
   allowed J's for K+1. */
			if((bbOOJ+0)>(bbOOh+1+bbOOSl)) bbOOflgb=1;
			if((bbOOJ+0)<(bbOOh+1-bbOOSl)) bbOOflgb=1;
/* Still OK?  Calculate transition */
			if(bbOOflgb==0){
				bbOOB=0;
				bbOOell=bbOOel+0;
				bbOO.T[0].S[bbOOff].f[bbOOjii]=\
					bbOO.Rh[0].EJ[bbOOehh]-\
					bbOO.Rl[0].EJ[bbOOell]; 
				bbOOcascadetemp=bbOOtrprob*bbOOfcf*\
					bbOO.Rh[0].CJ[bbOOehh]*bbOOSATT*\
					HL(bbOOB,bbOOJ,bbOOL,bbOOdL);
				bbOO.T[0].S[bbOOff].ci[bbOOjii]+=\
					bbOOcascadetemp;
				bbOO.Rl[0].CJ[bbOOell]+=bbOOfcfr*\
					bbOOcascadetemp;
				bbOOjii--;
				} 
if(bbOOdebugflag<DEBUG){
if(bbOOflgb!=0) fprintf(DBG," M M\n");
else{
fprintf(DBG,"%18.12e %18.12e (1c)\n",bbOO.T[0].S[bbOOff].f[bbOOjii+1],\
	bbOO.T[0].S[bbOOff].ci[bbOOjii+1]);
	} 
}
if(bbOOflgb!=0) fprintf(BBOOFILE," M M ");
else{
fprintf(BBOOFILE," %18.12e %18.12e ",bbOO.T[0].S[bbOOff].f[bbOOjii+1],\
	bbOO.T[0].S[bbOOff].ci[bbOOjii+1]);
	} 
/* Now, Delta-J=-1 (R) for Delta-K=+1.  */ 
			bbOOflgb=0; 
			if((bbOOh+1)<(bbOOLO[0].L)) bbOOflgb=1;
			if((bbOOJ-1)>(bbOOh+1+bbOOSl)) bbOOflgb=1;
			if((bbOOJ-1)<(bbOOh+1-bbOOSl)) bbOOflgb=1;
			if((bbOOJ-1)<0) bbOOflgb=1;
			if(bbOOflgb==0){
				bbOOB=-1;
				bbOOell=bbOOel-1;
				bbOO.T[0].S[bbOOff].f[bbOOjii]=\
					bbOO.Rh[0].EJ[bbOOehh]-\
					bbOO.Rl[0].EJ[bbOOell]; 
				bbOOcascadetemp=bbOOtrprob*bbOOfcf*\
					bbOO.Rh[0].CJ[bbOOehh]*bbOOSATT*\
				       HL(bbOOB,bbOOJ,bbOOL,bbOOdL);	
				bbOO.T[0].S[bbOOff].ci[bbOOjii]+=\
					bbOOcascadetemp;
				bbOO.Rl[0].CJ[bbOOell]+=bbOOfcfr*\
					bbOOcascadetemp;
				bbOOjii--;
				}
if(bbOOdebugflag<DEBUG){
if(bbOOflgb!=0) fprintf(DBG,"M M\n");
else{
fprintf(DBG,"%18.12e %18.12e (1d)\n",bbOO.T[0].S[bbOOff].f[bbOOjii+1],\
	bbOO.T[0].S[bbOOff].ci[bbOOjii+1]);
	} 
}
if(bbOOflgb!=0) fprintf(BBOOFILE,"M M\n");
else{
fprintf(BBOOFILE,"%18.12e %18.12e\n",bbOO.T[0].S[bbOOff].f[bbOOjii+1],\
	bbOO.T[0].S[bbOOff].ci[bbOOjii+1]);
	} 
/* Now, again for Delta-K=+0.  Delta-J=+1 (P), then (Q) and (R) */ 
			bbOOflgb=0;
			if((bbOOh+0)<(bbOOLO[0].L)) bbOOflgb=1;
/* Relative position in the low state E/P array */
			bbOOel=(bbOOh+0)*bbOO.nlJ + bbOO.nlJ-bbOO.nhJ + \
				(int)(bbOOSh-bbOOSl) + 0 + bbOOi; 
if(bbOOdebugflag<DEBUG){
fprintf(DBG,"\t(bbOOh+0)*bbOO.nlJ=%d, bbOO.nlJ-bbOO.nhJ=%d\n",\
		(bbOOh+0)*bbOO.nlJ,bbOO.nlJ-bbOO.nhJ);
fprintf(DBG,"\t(int)(bbOOSh-bbOOSl)=%d, bbOOi=%d, bbOOel=%d\n",\
		(int)(bbOOSh-bbOOSl),bbOOi,bbOOel);
fflush(DBG);
} 
			if((bbOOel+1)>=(bbOO.Rl[0].kc*bbOO.nlJ)) bbOOflgb=1;
/* Then check if the target value for J+1 is less than or greater than the 
   allowed J's for K+0. */
			if((bbOOJ+1)>(bbOOh+0+bbOOSl)) bbOOflgb=1;
			if((bbOOJ+1)<(bbOOh+0-bbOOSl)) bbOOflgb=1;
/* Still OK?  Calculate transition */
			if(bbOOflgb==0){
				bbOOB=+1;
				bbOOell=bbOOel+1;
				bbOO.T[0].S[bbOOff].f[bbOOjii]=\
					bbOO.Rh[0].EJ[bbOOehh]-\
					bbOO.Rl[0].EJ[bbOOell]; 
				bbOOcascadetemp=bbOOtrprob*bbOOfcf*\
					bbOO.Rh[0].CJ[bbOOehh]*bbOOSATT*\
					HL(bbOOB,bbOOJ,bbOOL,bbOOdL);
				bbOO.T[0].S[bbOOff].ci[bbOOjii]+=\
					bbOOcascadetemp;
				bbOO.Rl[0].CJ[bbOOell]+=bbOOfcfr*\
					bbOOcascadetemp;
				bbOOjii--;
				} 
if(bbOOdebugflag<DEBUG){
fprintf(DBG,"%d ",bbOOh);
if(bbOOflgb!=0) fprintf(DBG," %d %.1f (1c) M M\n",(bbOOh+0),bbOOJ);
else{
fprintf(DBG,"%d %.1f %18.12e %18.12e (1c)\n",(bbOOh+0),\
	bbOOJ,bbOO.T[0].S[bbOOff].f[bbOOjii+1],\
	bbOO.T[0].S[bbOOff].ci[bbOOjii+1]);
	} 
}
fprintf(BBOOFILE,"%d ",bbOOh);
if(bbOOflgb!=0) fprintf(BBOOFILE," %d %.1f M M ",(bbOOh+0),bbOOJ);
else{
fprintf(BBOOFILE,"%d %.1f %18.12e %18.12e ",(bbOOh+0),\
	bbOOJ,bbOO.T[0].S[bbOOff].f[bbOOjii+1],\
	bbOO.T[0].S[bbOOff].ci[bbOOjii+1]);
	} 
/* Now, Delta-J=0 (Q), Delta-K=+0 */ 
			bbOOflgb=0;
			if((bbOOh+0)<(bbOOLO[0].L)) bbOOflgb=1;
/* Then check if the target value for J+1 is less than or greater than the 
   allowed J's for K+0. */
			if((bbOOJ+0)>(bbOOh+0+bbOOSl)) bbOOflgb=1;
			if((bbOOJ+0)<(bbOOh+0-bbOOSl)) bbOOflgb=1;
/* Still OK?  Calculate transition */
			if(bbOOflgb==0){
				bbOOB=0;
				bbOOell=bbOOel+0;
				bbOO.T[0].S[bbOOff].f[bbOOjii]=\
					bbOO.Rh[0].EJ[bbOOehh]-\
					bbOO.Rl[0].EJ[bbOOell]; 
				bbOOcascadetemp=bbOOtrprob*bbOOfcf*\
					bbOO.Rh[0].CJ[bbOOehh]*\
					HL(bbOOB,bbOOJ,bbOOL,bbOOdL);
				bbOO.T[0].S[bbOOff].ci[bbOOjii]+=\
					bbOOcascadetemp;
				bbOO.Rl[0].CJ[bbOOell]+=bbOOfcfr*\
					bbOOcascadetemp;
				bbOOjii--;
				} 
if(bbOOdebugflag<DEBUG){
if(bbOOflgb!=0) fprintf(DBG," M M\n");
else{
fprintf(DBG,"%18.12e %18.12e (1c)\n",bbOO.T[0].S[bbOOff].f[bbOOjii+1],\
	bbOO.T[0].S[bbOOff].ci[bbOOjii+1]);
	} 
}
if(bbOOflgb!=0) fprintf(BBOOFILE," M M ");
else{
fprintf(BBOOFILE," %18.12e %18.12e ",bbOO.T[0].S[bbOOff].f[bbOOjii+1],\
	bbOO.T[0].S[bbOOff].ci[bbOOjii+1]);
	} 
/* Now, Delta-J=-1 (R) for Delta-K=+0.  */ 
			bbOOflgb=0; 
			if((bbOOh+0)<(bbOOLO[0].L)) bbOOflgb=1;
			if((bbOOJ-1)>(bbOOh+0+bbOOSl)) bbOOflgb=1;
			if((bbOOJ-1)<(bbOOh+0-bbOOSl)) bbOOflgb=1;
			if((bbOOJ-1)<0) bbOOflgb=1;
			if(bbOOflgb==0){
				bbOOB=-1;
				bbOOell=bbOOel-1;
				bbOO.T[0].S[bbOOff].f[bbOOjii]=\
					bbOO.Rh[0].EJ[bbOOehh]-\
					bbOO.Rl[0].EJ[bbOOell]; 
				bbOOcascadetemp=bbOOtrprob*bbOOfcf*\
					bbOO.Rh[0].CJ[bbOOehh]*bbOOSATT*\
				       HL(bbOOB,bbOOJ,bbOOL,bbOOdL);	
				bbOO.T[0].S[bbOOff].ci[bbOOjii]+=\
					bbOOcascadetemp;
				bbOO.Rl[0].CJ[bbOOell]+=bbOOfcfr*\
					bbOOcascadetemp;
				bbOOjii--;
				}
if(bbOOdebugflag<DEBUG){
if(bbOOflgb!=0) fprintf(DBG,"M M\n");
else{
fprintf(DBG,"%18.12e %18.12e (1d)\n",bbOO.T[0].S[bbOOff].f[bbOOjii+1],\
	bbOO.T[0].S[bbOOff].ci[bbOOjii+1]);
	} 
}
if(bbOOflgb!=0) fprintf(BBOOFILE,"M M\n");
else{
fprintf(BBOOFILE,"%18.12e %18.12e\n",bbOO.T[0].S[bbOOff].f[bbOOjii+1],\
	bbOO.T[0].S[bbOOff].ci[bbOOjii+1]);
	} 
/* Check Delta-J=+1 (P), then Delta-J=0 (Q) and Delta-J=-1 (R) for Delta-K=-1.
   First, check to see if the target K exists: */ 
			bbOOflgb=0;
			if((bbOOh-1)<(bbOOLO[0].L)) bbOOflgb=1;
/* Relative position in the low state E/P array */
			bbOOel=(bbOOh-1)*bbOO.nlJ + bbOO.nlJ-bbOO.nhJ + \
				(int)(bbOOSh-bbOOSl) + 1 + bbOOi; 
if(bbOOdebugflag<DEBUG){
fprintf(DBG,"\t(bbOOh-1)*bbOO.nlJ=%d, bbOO.nlJ-bbOO.nhJ=%d\n",\
		(bbOOh-1)*bbOO.nlJ,bbOO.nlJ-bbOO.nhJ);
fprintf(DBG,"\t(int)(bbOOSh-bbOOSl)=%d, bbOOi=%d, bbOOel=%d\n",\
		(int)(bbOOSh-bbOOSl),bbOOi,bbOOel);
fflush(DBG);
} 
			if((bbOOel+1)>=(bbOO.Rl[0].kc*bbOO.nlJ)) bbOOflgb=1;
/* Then check if the target value for J+1 is less than or greater than the 
   allowed J's for K-1. */
			if((bbOOJ+1)>(bbOOh-1+bbOOSl)) bbOOflgb=1;
			if((bbOOJ+1)<(bbOOh-1-bbOOSl)) bbOOflgb=1;
/* Still OK?  Calculate transition */
			if(bbOOflgb==0){
				bbOOB=+1;
				bbOOell=bbOOel+1;
				bbOO.T[0].S[bbOOff].f[bbOOjii]=\
					bbOO.Rh[0].EJ[bbOOehh]-\
					bbOO.Rl[0].EJ[bbOOell]; 
				bbOOcascadetemp=bbOOtrprob*bbOOfcf*\
					bbOO.Rh[0].CJ[bbOOehh]*bbOOSATT*\
					HL(bbOOB,bbOOJ,bbOOL,bbOOdL);
				bbOO.T[0].S[bbOOff].ci[bbOOjii]+=\
					bbOOcascadetemp;
				bbOO.Rl[0].CJ[bbOOell]+=bbOOfcfr*\
					bbOOcascadetemp;
				bbOOjii--;
				} 
if(bbOOdebugflag<DEBUG){
fprintf(DBG,"%d ",bbOOh);
if(bbOOflgb!=0) fprintf(DBG," %d %.1f M M ",(bbOOh-1),bbOOJ);
else{
fprintf(DBG,"%d %.1f %18.12e %18.12e (1e)\n",(bbOOh-1),\
	bbOOJ,bbOO.T[0].S[bbOOff].f[bbOOjii+1],\
	bbOO.T[0].S[bbOOff].ci[bbOOjii+1]);
	}
}
fprintf(BBOOFILE,"%d ",bbOOh);
if(bbOOflgb!=0) fprintf(BBOOFILE," %d %.1f M M ",(bbOOh-1),bbOOJ);
else{
fprintf(BBOOFILE,"%d %.1f %18.12e %18.12e ",(bbOOh-1),\
	bbOOJ,bbOO.T[0].S[bbOOff].f[bbOOjii+1],\
	bbOO.T[0].S[bbOOff].ci[bbOOjii+1]);
	} 
/* Now, Delta-J=+0 (Q),  Delta-K=-1. */ 
			bbOOflgb=0;
			if((bbOOh-1)<(bbOOLO[0].L)) bbOOflgb=1;
/* Then check if the target value for J+1 is less than or greater than the 
   allowed J's for K-1. */
			if((bbOOJ+0)>(bbOOh-1+bbOOSl)) bbOOflgb=1;
			if((bbOOJ+0)<(bbOOh-1-bbOOSl)) bbOOflgb=1;
/* Still OK?  Calculate transition */
			if(bbOOflgb==0){
				bbOOB=0;
				bbOOell=bbOOel+0;
				bbOO.T[0].S[bbOOff].f[bbOOjii]=\
					bbOO.Rh[0].EJ[bbOOehh]-\
					bbOO.Rl[0].EJ[bbOOell]; 
				bbOOcascadetemp=bbOOtrprob*bbOOfcf*\
					bbOO.Rh[0].CJ[bbOOehh]*bbOOSATT*\
					HL(bbOOB,bbOOJ,bbOOL,bbOOdL);
				bbOO.T[0].S[bbOOff].ci[bbOOjii]+=\
					bbOOcascadetemp;
				bbOO.Rl[0].CJ[bbOOell]+=bbOOfcfr*\
					bbOOcascadetemp;
				bbOOjii--;
				} 
if(bbOOdebugflag<DEBUG){
if(bbOOflgb!=0) fprintf(DBG," M M ");
else{
fprintf(DBG,"%18.12e %18.12e \n",bbOO.T[0].S[bbOOff].f[bbOOjii+1],\
	bbOO.T[0].S[bbOOff].ci[bbOOjii+1]);
	}
}
if(bbOOflgb!=0) fprintf(BBOOFILE," M M ");
else{
fprintf(BBOOFILE,"%18.12e %18.12e ",bbOO.T[0].S[bbOOff].f[bbOOjii+1],\
	bbOO.T[0].S[bbOOff].ci[bbOOjii+1]);
	} 
/* Now, Delta-J=-1 (R) for Delta-K=-1.  */ 
			bbOOflgb=0; 
			if((bbOOh-1)<(bbOOLO[0].L)) bbOOflgb=1;
			if((bbOOJ-1)>(bbOOh-1+bbOOSl)) bbOOflgb=1;
			if((bbOOJ-1)<(bbOOh-1-bbOOSl)) bbOOflgb=1;
			if((bbOOJ-1)<0) bbOOflgb=1;
			if(bbOOflgb==0){
				bbOOB=-1;
				bbOOell=bbOOel-1;
				bbOO.T[0].S[bbOOff].f[bbOOjii]=\
					bbOO.Rh[0].EJ[bbOOehh]-\
					bbOO.Rl[0].EJ[bbOOell]; 
				bbOOcascadetemp=bbOOtrprob*bbOOfcf*\
					bbOO.Rh[0].CJ[bbOOehh]*\
				       HL(bbOOB,bbOOJ,bbOOL,bbOOdL);	
				bbOO.T[0].S[bbOOff].ci[bbOOjii]+=\
					bbOOcascadetemp;
				bbOO.Rl[0].CJ[bbOOell]+=bbOOfcfr*\
					bbOOcascadetemp; 
				}
if(bbOOdebugflag<DEBUG){
if(bbOOflgb!=0) fprintf(DBG,"M M\n");
else{
fprintf(DBG,"%18.12e %18.12e (1f)\n",bbOO.T[0].S[bbOOff].f[bbOOjii],\
	bbOO.T[0].S[bbOOff].ci[bbOOjii]);
	} 
}
if(bbOOflgb!=0) fprintf(BBOOFILE,"M M\n");
else{
fprintf(BBOOFILE,"%18.12e %18.12e\n",bbOO.T[0].S[bbOOff].f[bbOOjii],\
	bbOO.T[0].S[bbOOff].ci[bbOOjii]);
	} 
			} /* close if bbOOflga is zero. */
		} /* close upper-state J-of-K loop */ 
	} /* close upper-state K loop */ 
	fclose(BBOOFILE);
		} /* close if high-state v is cascade populated */ 
	} /* close low-state v */
	} /* close high-state v */
return; 
}

