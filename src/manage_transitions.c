#include "rvesim.h"


/***************** manage_transitions *****************/

/* This function identifies each transition's type and calls the correct 
   function to do the transition. */

void manage_transitions(){

int mta=0,mtb=0,mtc=0,mtd=0,mte=0,mtdone=1,mtcasesum=0,*mtz;
double *mtTeH,*mtTeL,mtEmax=0;
int mtcounter=0,mtsorth=0,mtsortdum=0;

if(DEBUG>mtdebugflag){
fprintf(DBG,"NUMMOL is %d\n",NUMMOL);
fflush(DBG);
}
/* begin loop over molecules */
for(mta=0;mta<NUMMOL;mta++){
/* begin loop over transitions */
	mtz=(int*)calloc(MOL[mta].trans,sizeof(int));
	mtTeH=(double*)calloc(MOL[mta].trans,sizeof(double));
	mtTeL=(double*)calloc(MOL[mta].trans,sizeof(double)); 
/* the following is a mechanism for sorting the transitions from
   highest energy to lowest energy -- this is necessary to get the
   cascade right.  The algorithm is probably primitive and inefficient,
   so if you want to improve it, go right ahead... */
	for(mtb=0;mtb<MOL[mta].trans;mtb++){ /* save Te's to temporary 
						arrays */
		mtz[mtb]=mtb;
		if(MOL[mta].s[MOL[mta].t[mtb].Hi].Ca>-1){
			mtTeH[mtb]=CA[MOL[mta].s[MOL[mta].t[mtb].Hi].Ca].Te;}
		if(MOL[mta].s[MOL[mta].t[mtb].Hi].Cb>-1){
			mtTeH[mtb]=CB[MOL[mta].s[MOL[mta].t[mtb].Hi].Cb].Te;}
		if(MOL[mta].s[MOL[mta].t[mtb].Lo].Ca>-1){
			mtTeL[mtb]=CA[MOL[mta].s[MOL[mta].t[mtb].Lo].Ca].Te;}
		if(MOL[mta].s[MOL[mta].t[mtb].Lo].Cb>-1){
			mtTeL[mtb]=CB[MOL[mta].s[MOL[mta].t[mtb].Lo].Cb].Te;}
if(DEBUG>mtdebugflag){
fprintf(DBG,"MOL[%d].s[%d].Ca (Hi) is %d",mta,MOL[mta].t[mtb].Hi,\
		MOL[mta].s[MOL[mta].t[mtb].Hi].Ca); 
fprintf(DBG,"MOL[%d].s[%d].Cb (Hi) is %d\n",mta,MOL[mta].t[mtb].Hi,\
		MOL[mta].s[MOL[mta].t[mtb].Hi].Cb); 
fprintf(DBG,"MOL[%d].s[%d].Ca (Lo) is %d",mta,MOL[mta].t[mtb].Lo,\
		MOL[mta].s[MOL[mta].t[mtb].Lo].Ca); 
fprintf(DBG,"MOL[%d].s[%d].Cb (Lo) is %d\n",mta,MOL[mta].t[mtb].Lo,\
		MOL[mta].s[MOL[mta].t[mtb].Lo].Cb); 
fprintf(DBG,"mtTeH[%d] is %f; mtTeL[%d} is %f\n",\
		mtb,mtTeH[mtb],mtb,mtTeL[mtb]);
fflush(DBG);
}
		} 
	mtdone=1;
	if(MOL[mta].trans==1){mtdone=0;}
if(DEBUG>mtdebugflag){
fprintf(DBG,"Before sort.  mtdone is %d\n ",mtdone);
fflush(DBG);
}
mtcounter=0;
	while(mtdone!=0){ /* this first loop orders the term energy for
			     the high state */
		for(mtb=0;mtb<(MOL[mta].trans-1);mtb++){ 
if(DEBUG>mtdebugflag){
fprintf(DBG,"First sort for-loop: mtb=%d, mtz[mtb]=%d, mtTeH[%d]=%f\n",\
		mtb,mtz[mtb],mtb,mtTeH[mtb]);
fflush(DBG);
}
			mtsorth=mtsortdum=mtz[mtb];
			mtEmax=mtTeH[mtz[mtb]];
			for(mtc=mtb;mtc<(MOL[mta].trans);mtc++){ 
				if(mtTeH[mtz[mtc]]>mtEmax) {
					mtsorth=mtz[mtc]; 
					mtEmax=mtTeH[mtz[mtc]];
					}
if(DEBUG>mtdebugflag){
fprintf(DBG,"scan mtb=%d, mtc=%d, mtz[mtc]=%d, mtTeH[%d]=%f\n",\
		mtb,mtc,mtz[mtc],mtz[mtc],mtTeH[mtz[mtc]]);
fprintf(DBG,"\t\t mtsorth=%d, mtEmax=%f\n",mtsorth,mtEmax);
fflush(DBG);
}
				}
			mtz[mtb]=mtsorth;

if(DEBUG>mtdebugflag){
	for(mtc=0;mtc<mtb+2;mtc++){
fprintf(DBG,"mtc=%d, mtz[mtc]=%d, mtTeH[%d]=%f\n",\
		mtc,mtz[mtc],mtz[mtc],mtTeH[mtz[mtc]]);
	}
fflush(DBG);
}

			for(mtc=(mtb+1);mtc<(MOL[mta].trans);mtc++){ 
				if(mtz[mtc]==mtsorth) {
					mtz[mtc]=mtsortdum; 
					}
if(DEBUG>mtdebugflag){
fprintf(DBG,"assign loop: mtb=%d, mtc=%d, mtz[mtb]=%d, mtz[mtc]=%d\n",\
		mtb,mtc,mtz[mtb],mtz[mtc]);
fflush(DBG);
}
				} 
			} 
		mtdone=0;
		for(mtb=0;mtb<(MOL[mta].trans-1);mtb++){ 
			if(mtTeH[mtz[mtb]]<mtTeH[mtz[mtb+1]]){
				mtdone=1;
				}
if(DEBUG>mtdebugflag){
fprintf(DBG,"sort check: mtz[%d]=%d, mtTeH[%d]=%f, mtTeH[%d]=%f\n",\
	mtb,mtz[mtb],mtz[mtb],mtTeH[mtz[mtb]],mtz[mtb+1],mtTeH[mtz[mtb+1]]);
fflush(DBG);
}
			}
if(DEBUG>mtdebugflag){
fprintf(DBG,"sort check: mtz[%d]=%d, mtTeH[%d]=%f\n",\
	mtb,mtz[mtb],mtz[mtb],mtTeH[mtz[mtb]]);
fflush(DBG);
}
		}
	for(mtd=0;mtd<MOL[mta].trans-1;mtd++){ /* this loop orders the 
					  term energy for the low state */
		mtdone=1;
		mte=mtd;
		while(mtdone!=0){
			if(mtTeH[mtz[mte]]!=mtTeH[mtz[mte+1]]){mtdone=0;}
			else mte++;
			}
		for(mtb=mtd;mtb<mte;mtb++){ 
if(DEBUG>mtdebugflag){
fprintf(DBG,"low sort mtb=%d, mtz[mtb]=%d, mtTeL[%d]=%f\n",\
		mtb,mtz[mtb],mtb,mtTeL[mtb]);
fprintf(DBG,"\t\t mtz[mtb+1]=%d, mtTeL[%d]=%f\n",\
		mtz[mtb+1],mtb+1,mtTeL[mtb+1]);
fflush(DBG);
}
			mtsorth=mtsortdum=mtz[mtb];
			mtEmax=mtTeL[mtz[mtb]];
			for(mtc=mtb;mtc<=mte;mtc++){ 
				if(mtTeL[mtz[mtc]]>mtEmax) {
					mtsorth=mtz[mtc]; 
					mtEmax=mtTeL[mtz[mtc]];
					}
if(DEBUG>mtdebugflag){
fprintf(DBG,"scan mtb=%d, mtc=%d, mtz[mtc]=%d, mtTeL[%d]=%f\n",\
		mtb,mtc,mtz[mtc],mtz[mtc],mtTeL[mtz[mtc]]);
fprintf(DBG,"\t\t mtsorth=%d, mtEmax=%f\n",mtsorth,mtEmax);
fflush(DBG);
}
				}
			mtz[mtb]=mtsorth;

if(DEBUG>mtdebugflag){
	for(mtc=mtd;mtc<mtb+2;mtc++){
fprintf(DBG,"mtc=%d, mtz[mtc]=%d, mtTeL[%d]=%f\n",\
		mtc,mtz[mtc],mtz[mtc],mtTeL[mtz[mtc]]);
	}
fflush(DBG);
}

			for(mtc=(mtb+1);mtc<=mte;mtc++){ 
				if(mtz[mtc]==mtsorth) {
					mtz[mtc]=mtsortdum; 
					}
if(DEBUG>mtdebugflag){
fprintf(DBG,"assign loop: mtb=%d, mtc=%d, mtz[mtb]=%d, mtz[mtc]=%d\n",\
		mtb,mtc,mtz[mtb],mtz[mtc]);
fflush(DBG);
}
				} 
			}
		mtdone=0;
		for(mtb=mtd;mtb<(mte);mtb++){ 
			if((mtTeH[mtz[mtb]]==mtTeH[mtz[mtb+1]])&&\
				(mtTeL[mtz[mtb]]<mtTeL[mtz[mtb+1]])){
				mtdone=1;
				}
if(DEBUG>mtdebugflag){
fprintf(DBG,"sort check low mtb=%d, mtz[mtb]=%d, mtTeH[%d]=%f, ",\
		mtb,mtz[mtb],mtz[mtb],mtTeH[mtz[mtb]]);
fprintf(DBG,"mtTeL[%d]=%f\n",mtz[mtb],mtTeL[mtz[mtb]]);
fprintf(DBG,"\t\t mtz[mtb+1]=%d, mtTeL[%d]=%f\n",\
		mtz[mtb+1],mtz[mtb+1],mtTeL[mtz[mtb+1]]);
fflush(DBG);
}
			}
if(DEBUG>mtdebugflag){
fprintf(DBG,"sort check low mtb=%d, mtz[mtb]=%d, mtTeH[%d]=%f, ",\
		mtb,mtz[mtb],mtz[mtb],mtTeH[mtz[mtb]]);
fprintf(DBG,"mtTeL[%d]=%f\n",mtz[mtb],mtTeL[mtz[mtb]]);
fflush(DBG);
} 
		mtd=mte;
		}
	for(mtc=0;mtc<MOL[mta].trans;mtc++){ 
		mtb=mtz[mtc]; /* do transitions in order determined in
				 the last loop */
/* decide which type of transition it is and call that one */
/* note:  lines for transitions involving Case c and d are included
   here for "forward compatibility" -- they can't be used yet.  But,
   having them in the code shouldn't cause a problem. */ 
		if(MOL[mta].s[MOL[mta].t[mtb].Hi].Ca>-1) mtcasesum=1;
		if(MOL[mta].s[MOL[mta].t[mtb].Hi].Cb>-1) mtcasesum=2;
		if(MOL[mta].s[MOL[mta].t[mtb].Hi].Cc>-1) mtcasesum=4;
		if(MOL[mta].s[MOL[mta].t[mtb].Hi].Cd>-1) mtcasesum=8;
if(DEBUG>mtdebugflag){
fprintf(DBG,"mtcasenum for mtb=%d (mtz[%d]=%d) is %d\n",\
		mtb,mtc,mtz[mtc],mtcasesum);
fflush(DBG);
}

		if(MOL[mta].s[MOL[mta].t[mtb].Lo].Ca>-1) mtcasesum+=1;
		if(MOL[mta].s[MOL[mta].t[mtb].Lo].Cb>-1) mtcasesum+=2;
		if(MOL[mta].s[MOL[mta].t[mtb].Lo].Cc>-1) mtcasesum+=4;
		if(MOL[mta].s[MOL[mta].t[mtb].Lo].Cd>-1) mtcasesum+=8;
if(DEBUG>mtdebugflag){ 
fprintf(DBG,"mtcasenum for mtb=%d (mtz[%d]=%d) is %d\n",\
		mtb,mtc,mtz[mtc],mtcasesum);
fflush(DBG);
}

		switch(mtcasesum){ 
			case 2:
if(DEBUG>mtdebugflag){
fprintf(DBG,"a_a transition called for\n\t MOL[%d].s[%d].Name = %s to ",\
	mta,(MOL[mta].t[mtb].Hi),MOL[mta].s[MOL[mta].t[mtb].Hi].Name); 
fprintf(DBG,"MOL[%d].s[%d].Name = %s\n",\
	mta,(MOL[mta].t[mtb].Lo),MOL[mta].s[MOL[mta].t[mtb].Lo].Name); 
fflush(DBG);
}
				a_a(mta,mtb); 
				break;
			case 3:
if(DEBUG>mtdebugflag){
fprintf(DBG,"a_b transition called for\n\t MOL[%d].s[%d].Name = %s to ",\
	mta,(MOL[mta].t[mtb].Hi),MOL[mta].s[MOL[mta].t[mtb].Hi].Name); 
fprintf(DBG,"MOL[%d].s[%d].Name = %s\n",\
	mta,(MOL[mta].t[mtb].Lo),MOL[mta].s[MOL[mta].t[mtb].Lo].Name); 
fflush(DBG);
}
				/*a_b(mta,mtb);*/
				break;
			case 5:
if(DEBUG>mtdebugflag){
fprintf(DBG,"a_c transition called for\n\t MOL[%d].s[%d].Name = %s to ",\
	mta,(MOL[mta].t[mtb].Hi),MOL[mta].s[MOL[mta].t[mtb].Hi].Name); 
fprintf(DBG,"MOL[%d].s[%d].Name = %s\n",\
	mta,(MOL[mta].t[mtb].Lo),MOL[mta].s[MOL[mta].t[mtb].Lo].Name); 
fflush(DBG);
}
				/*a_c(mta,mtb);*/
				break;
			case 9:
if(DEBUG>mtdebugflag){
fprintf(DBG,"a_d transition called for\n\t MOL[%d].s[%d].Name = %s to ",\
	mta,(MOL[mta].t[mtb].Hi),MOL[mta].s[MOL[mta].t[mtb].Hi].Name); 
fprintf(DBG,"MOL[%d].s[%d].Name = %s\n",\
	mta,(MOL[mta].t[mtb].Lo),MOL[mta].s[MOL[mta].t[mtb].Lo].Name); 
fflush(DBG);
}
				/*a_d(mta,mtb);*/
				break;
			case 4:
if(DEBUG>mtdebugflag){
fprintf(DBG,"b_b transition called for\n\t MOL[%d].s[%d].Name = %s to ",\
	mta,(MOL[mta].t[mtb].Hi),MOL[mta].s[MOL[mta].t[mtb].Hi].Name); 
fprintf(DBG,"MOL[%d].s[%d].Name = %s\n",\
	mta,(MOL[mta].t[mtb].Lo),MOL[mta].s[MOL[mta].t[mtb].Lo].Name); 
fprintf(DBG,"b_b transition called for\n\t MOL[%d].t[%d].Nhi = %s to ",\
	mta,mtb,MOL[mta].t[mtb].Nhi); 
fprintf(DBG,"b_b transition called for\n\t MOL[%d].t[%d].Nlo = %s \n",\
	mta,mtb,MOL[mta].t[mtb].Nlo); 
fflush(DBG);
}
				b_b(mta,mtb);
				break;
			case 6:
if(DEBUG>mtdebugflag){
fprintf(DBG,"b_c transition called for\n\t MOL[%d].s[%d].Name = %s to ",\
	mta,(MOL[mta].t[mtb].Hi),MOL[mta].s[MOL[mta].t[mtb].Hi].Name); 
fprintf(DBG,"MOL[%d].s[%d].Name = %s\n",\
	mta,(MOL[mta].t[mtb].Lo),MOL[mta].s[MOL[mta].t[mtb].Lo].Name); 
fflush(DBG);
}
				/*b_c(mta,mtb);*/
				break;
			case 10:
if(DEBUG>mtdebugflag){
fprintf(DBG,"b_d transition called for\n\t MOL[%d].s[%d].Name = %s to ",\
	mta,(MOL[mta].t[mtb].Hi),MOL[mta].s[MOL[mta].t[mtb].Hi].Name); 
fprintf(DBG,"MOL[%d].s[%d].Name = %s\n",\
	mta,(MOL[mta].t[mtb].Lo),MOL[mta].s[MOL[mta].t[mtb].Lo].Name); 
fflush(DBG);
}
				/*b_d(mta,mtb);*/
				break;
			case 8:
if(DEBUG>mtdebugflag){
fprintf(DBG,"c_c transition called for\n\t MOL[%d].s[%d].Name = %s to ",\
	mta,(MOL[mta].t[mtb].Hi),MOL[mta].s[MOL[mta].t[mtb].Hi].Name); 
fprintf(DBG,"MOL[%d].s[%d].Name = %s\n",\
	mta,(MOL[mta].t[mtb].Lo),MOL[mta].s[MOL[mta].t[mtb].Lo].Name); 
fflush(DBG);
}
				/*c_c(mta,mtb);*/
				break;
			case 12:
if(DEBUG>mtdebugflag){
fprintf(DBG,"c_d transition called for\n\t MOL[%d].s[%d].Name = %s to ",\
	mta,(MOL[mta].t[mtb].Hi),MOL[mta].s[MOL[mta].t[mtb].Hi].Name); 
fprintf(DBG,"MOL[%d].s[%d].Name = %s\n",\
	mta,(MOL[mta].t[mtb].Lo),MOL[mta].s[MOL[mta].t[mtb].Lo].Name); 
fflush(DBG);
}
				/*c_d(mta,mtb);*/
				break;
			case 16:
if(DEBUG>mtdebugflag){
fprintf(DBG,"d_d transition called for\n\t MOL[%d].s[%d].Name = %s to ",\
	mta,(MOL[mta].t[mtb].Hi),MOL[mta].s[MOL[mta].t[mtb].Hi].Name); 
fprintf(DBG,"MOL[%d].s[%d].Name = %s\n",\
	mta,(MOL[mta].t[mtb].Lo),MOL[mta].s[MOL[mta].t[mtb].Lo].Name); 
fflush(DBG);
}
				/*d_d(mta,mtb);*/
				break; 
			default:
printf("Error in manage_transitions at mta=%d, mtb=%d.  Exiting.\n",mta,mtb);
				exit(1);
			} /* close switch case */
		} /* close loop over transitions */
	free(mtz);
	} /* close loop over molecules */
return;
}


