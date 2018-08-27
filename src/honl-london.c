#include "rvesim.h"

/***************** HL *****************/
  
/* This function calculates the RELATIVE Holn-London factor for all
   transitions (called from all transition functions). The equations
   here are from Herzberg and expressed conveniently for emission.  
 
   For consistency with the currently written transition functions, the
   integer for P is +1, for Q is 0 and for R is -1.  Since Delta-Lambda
   isn't used in those functions, they will be counted from the lower
   state.....  confusion.....*/
double HL(int HB, double HJ, int HL, int HdL){

double Hfac=0;

if(HdL==0){ /* if Delta-Lambda =0 */
	if(HB==+1){ /* P-branch */
		Hfac=(HJ+1+HL)*(HJ+1-HL)/(2*HJ*HJ+3*HJ+1);
		}
	if(HB==0){ /* Q-branch */
		if(HJ==0) Hfac=0;
		else Hfac=HL*HL/(HJ*(HJ+1));
		}
	if(HB==-1){ /* R-branch */
		if(HJ==0) Hfac=0;
		else Hfac=(HJ+HL)*(HJ-HL)/(HJ*(2*HJ+1));
		}
	} 
if(HdL==+1){ /* if the low state is one less than the high state */
	if(HB==+1){ /* P-branch */
		Hfac=0.5*(HJ+1-HL)*(HJ+2-HL)/(2*HJ*HJ+3*HJ+1);
		}
	if(HB==0){ /* Q-branch */
		if(HJ==0) Hfac=0;
		else Hfac=0.5*(HJ+HL)*(HJ+1-HL)/(HJ*(HJ+1));
		}
	if(HB==-1){ /* R-branch */
		if(HJ==0) Hfac=0;
		else Hfac=0.5*(HJ+HL)*(HJ-1+HL)/(HJ*(2*HJ+1));
		}
	} 
if(HdL==-1){ /* if the low state is one more than the high state */
	if(HB==+1){ /* P-branch */
		Hfac=0.5*(HJ+1+HL)*(HJ+2+HL)/(2*HJ*HJ+3*HJ+1);
		}
	if(HB==0){ /* Q-branch */
		if(HJ==0) Hfac=0;
		else Hfac=0.5*(HJ-HL)*(HJ+1+HL)/(HJ*(HJ+1));
		}
	if(HB==-1){ /* R-branch */
		if(HJ==0) Hfac=0;
		else Hfac=0.5*(HJ-HL)*(HJ-1-HL)/(HJ*(2*HJ+1));
		}
	}
if(HLswitch==0) return Hfac;
else return 1;
}

