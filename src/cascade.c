#include "rvesim.h"

/****************  case_a_cascade *****************/

/* This function calculates energy levels only when needed for destination
   states in a transition.  This function is called from a-a.c and a-b.c.
   For multiplet non-Sigma states, Hund's Case (a). 

   The input variables are (1) the rotset where the energy levels should
   be stored, (2) the position in the case-a info array, (3) the value of
   Omega for this set of J's, (4) the vibrational quantum number, (5) the
   current highest value of J calculated, (6) the highest value needed */
void 
case_a_cascade(rotset *cacR,int cacA,double cacO,int cacV,int cacJc,int cacJn){
int cacb=0,cacL=0;
double cacJ=0,cacDv=0,cacBe=0,cacAe=0,cacwe=0,cacwexe=0;
double cacbeta=0,cacTe=0,cacDe=0,cacBv=0,cacTeVib=0;

if(DEBUG>cacdebugflag){
fprintf(DBG,"\n\nTop of case_a_cascade\n\n");\
fflush(DBG);
} 
cacBe=CA[cacA].Be;
cacAe=CA[cacA].ae; 
cacTe=CA[cacA].Te; 
cacwe=CA[cacA].we; 
cacwexe=CA[cacA].wexe; 
cacL=CA[cacA].L; 
cacbeta=CA[cacA].beta; 
if(DEBUG>cacdebugflag){
fprintf(DBG,"Be=%f ae=%f Te=%f we=%f ",cacBe,cacAe,cacTe,cacwe);
fprintf(DBG,"wexe=%f cacL=%d, cacbeta=%f\n",cacwexe,cacL,cacbeta);
fflush(DBG);
}
	
/* Calculate some stuff that will be useful later */
cacDe=4*pow(cacBe,3)/(cacwe*cacwe); /* Centrifugal distortion */
cacBv=cacBe-cacAe*(cacV+0.5);
cacDv=cacDe+cacbeta*(cacV+0.5); 
cacTeVib=cacTe+(cacV+0.5)*cacwe-(cacV+0.5)*(cacV+0.5)*cacwexe; 
if(DEBUG>cacdebugflag){
fprintf(DBG,"cacO=%f cacDe=%f cacBv=%f cacDv=%f cacTeVib=%f\n",\
		cacO,cacDe,cacBv,cacDv,cacTeVib);
fflush(DBG);
} 

/* There aren't any extra E's here.  If someone ever wants to add in
   a splitting factor for the +/- wavefunctions, then double the space 
   available (cacJn) and change the indexing and energy calculation
   in the loops below.*/
cacR[0].J=(double*)realloc(cacR[0].J,(cacJn*sizeof(double)));
cacR[0].EJ=(double*)realloc(cacR[0].EJ,(cacJn*sizeof(double))); 
cacR[0].CJ=(double*)realloc(cacR[0].CJ,(cacJn*sizeof(double))); 
for(cacb=cacJc;cacb<cacJn;cacb++){
	cacJ=cacb+cacO;
	cacR[0].J[cacb]=cacJ;
	cacR[0].EJ[cacb]=cacTeVib+cacBv*(cacJ*(cacJ+1)-cacO*cacO)-\
		cacDv*cacJ*cacJ*(cacJ+1)*(cacJ+1); 	
	cacR[0].CJ[cacb]=0;
if(DEBUG>cacdebugflag){
fprintf(DBG,\
"cacb=%d cacR[0].J[cacb]=%f cacR[0].EJ[cacb]=%f cacR[0].CJ[cacb]=%f\n",\
		cacb,cacJ,cacR[0].EJ[cacb],cacR[0].CJ[cacb]);
fflush(DBG);
}
	} 
cacR[0].jc=cacJn; 
return; 
}

/**************** case_b_cascade_b *****************/

/* This function calculates energy levels only when needed for destination
   states in a transition.  This function is called from b-b.c. 

   The input variables are (1) the rotset where the energy levels should
   be stored, (2) the position in the case-b info array, (3) the vibrational 
   quantum number, (4) the current highest value of K calculated, (5) the 
   highest value needed, (6) the multiplicity (J's per K) */
void 
case_b_cascade(rotset *cbcR,int cbcA,int cbcV,int cbcKc,int cbcKn,int cbcM){
int cbcb=0,cbcc=0,cbcL=0;
double cbcK=0,cbcDv=0,cbcBe=0,cbcAe=0,cbcwe=0,cbcwexe=0;
double cbcbeta=0,cbcTe=0,cbcDe=0,cbcBv=0,cbcTeVib=0;

if(DEBUG>cbcdebugflag){
fprintf(DBG,"\n\nTop of case_b_cascade\n\n");\
fprintf(DBG,"cbcA=%d, cbcV=%d, cbcKc=%d, cbcKn=%d\n",\
		cbcA,cbcV,cbcKc,cbcKn);
fprintf(DBG,"cbcR[0].j=%d\n",cbcR[0].j);
fflush(DBG);
} 
cbcBe=CB[cbcA].Be;
cbcAe=CB[cbcA].ae; 
cbcTe=CB[cbcA].Te; 
cbcwe=CB[cbcA].we; 
cbcwexe=CB[cbcA].wexe; 
cbcL=CB[cbcA].L; 
cbcbeta=CB[cbcA].beta; 
if(DEBUG>cbcdebugflag){
fprintf(DBG,"Be=%f ae=%f Te=%f we=%f ",cbcBe,cbcAe,cbcTe,cbcwe);
fprintf(DBG,"wexe=%f cbcL=%d, cbcbeta=%f\n",cbcwexe,cbcL,cbcbeta);
fflush(DBG);
}
	
/* Calculate some stuff that will be useful later */
cbcDe=4*pow(cbcBe,3)/(cbcwe*cbcwe); /* Centrifugal distortion */
cbcBv=cbcBe-cbcAe*(cbcV+0.5);
cbcDv=cbcDe+cbcbeta*(cbcV+0.5); 
cbcTeVib=cbcTe+(cbcV+0.5)*cbcwe-(cbcV+0.5)*(cbcV+0.5)*cbcwexe; 
if(DEBUG>cbcdebugflag){
fprintf(DBG,"cbcDe=%f cbcBv=%f cbcDv=%f cbcTeVib=%f\n",\
		cbcDe,cbcBv,cbcDv,cbcTeVib);
fflush(DBG);
} 

/* There aren't any extra E's here.  If someone ever wants to add in
   a splitting factor for the +/- wavefunctions, then double the space 
   available (cbcKn) and change the indexing and energy calculation
   in the loops below.*/
cbcR[0].J=(double*)realloc(cbcR[0].J,(cbcM*cbcKn*sizeof(double)));
cbcR[0].EJ=(double*)realloc(cbcR[0].EJ,(cbcM*cbcKn*sizeof(double))); 
cbcR[0].CJ=(double*)realloc(cbcR[0].CJ,(cbcM*cbcKn*sizeof(double))); 
for(cbcb=cbcKc;cbcb<cbcKn;cbcb++){
	cbcK=cbcb;
	for(cbcc=0;cbcc<cbcM;cbcc++){
		cbcR[0].J[cbcb*cbcM+cbcc]=cbcK;
		cbcR[0].EJ[cbcb*cbcM+cbcc]=cbcTeVib+cbcBv*(cbcK*(cbcK+1))-\
			cbcDv*cbcK*cbcK*(cbcK+1)*(cbcK+1); 	
		cbcR[0].CJ[cbcb*cbcM+cbcc]=0;
if(DEBUG>cbcdebugflag){
fprintf(DBG,\
"cbcb=%d cbcR[0].J[cbcb]=%f cbcR[0].EJ[cbcb]=%f cbcR[0].CJ[cbcb]=%f\n",\
		cbcb,cbcK,cbcR[0].EJ[cbcb],cbcR[0].CJ[cbcb]);
fflush(DBG);
}
		}
	} 
cbcR[0].kc=cbcKn; 
return;
}

