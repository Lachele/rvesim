/* These are the headers for the RVESIM.c program */

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h> 
#include <ctype.h>

#define SUB (double)0.5
#define PI 3.14159265358979323846264338
#define JKm 0.0504507857956213 /* constant for an equation trick */
#define JKb 3.04614091081434 /* another constant for an equation trick */
#define kB 1.38066e-23/1.9864e-23  /* Boltzmann's constant in wavenumbers
per kelvin.  The first number is Boltzmann's constant in J/K (joules per 
kelvin).  The second number is the conversion from joules to wavenumbers. */

/****************** Begin Section on States ******************/

typedef struct {
	char f[201]; /* filename */
	FILE *F; /* file pointer */
	} fileset;

typedef struct{
	int vlo, vnum; /* starting vibrational quantum number and number
			  of quantum numbers to consider */
	int vclo, vcnum; /* starting vib quant number for cascade population 
		and number of quanta to consider -- the population will all be 
		contained in the rotset, but need to know about it */
	double *p;  /* relative populations for each of the vnum levels */ 
	int *c; /* change flag for the relative population: (0) not
		   calculated, (1) calculated (-1) re-calculated or 
		   (-2) is awaiting recalculation. */
	} vset;

typedef struct{
	int c; /* change flag for this set */
	int j,k,jc,kc; /* number of J and/or K values in this set jc and kc 
		are there because there might be more E entries than there 
		are P entries -- there will also be more CJ/CK entries.*/
	double *J, *K; /* values of J and/or K for each level */
	double Jdissoc,Kdissoc; /* values of J/K where the molecule 
		   dissociates, according to the ratio of Bv to Dv. */
	double *EJ, *EK; /* the energies of this level (wavenumbers) */
	double *PJ, *PK, Jmaxp; /* relative populations of this level, not 
		including any contribution from cascade and the value of J 
		that was found to be the true maximum value -- 10 June 02:
		I just commented out all use of Jmaxp */
	double *CJ, *CK; /* cascade population of this level */
	} rotset;

/*** The `State' structure ***/
typedef struct {
	char Name[21]; /* The letter state designation (a', B, c, etc.). */
	double pop; /* relative population for this state */
	int cp, cC, cD, co, cT, cS, cr, cv; /* change flags for population 
		(cp), Case (cC), presence of rotational input file (cD), 
		number omegas (co), temperature (cT), spectroscopic 
		information (cS), rotation (cr) and vibration (cv) info */
	char Case[11]; /* Hund's Case.  a or b */ 
	int Ca; /* numerical position in stateinfo array if Case a */
	int Cb; /* numerical position in stateinfo array if Case b */
	int Cc, Cd;  /* similar to Ca and Cb */
	int Dist; /* is the distribution of states is from a file (-1), from
		     the temperature (1), or not relevant (0). */ 
	vset *v;  /* pointer to set of vibrational information 
		   MEMORY ALLOCATED in read_trans.c */ 
	rotset *r; /* pointer to set of rotational distribution info 
		    MEMORY ALLOCATED in rotation-distribution functions
		    (e.g., multi_nonzero_JKb.c) */ 
	int nr,nO,nV; /* number of rotsets, number Omegas and vnum (this will
		      change if different Omegas get different vnums...) */
	rotset *rc; /* pointer to rotsets for additional vibrational levels
		that are populated only by cascade -- the number of these 
		sets is 30 (or Onum*30) because it's hard to know which 
		levels will get populated from which states, and I'm too 
		lazy to shuffle all the array entries. */ 
	int pmsymm, Isymm; /* if Lambda==0:  Are even numbered K's + or - 
		(pmsymm); and, if this set is homonuclear, are even K's 
		symmetric or antisymmetric (Isymm). */
	int no;  /* How many times is this state an origination state? */
	int *o; /* the list of positions in the transition array where this
		   state is the origination state. 
		 MEMORY ALLOCATED in read_trans.c */
	int nd; /* How many times is this state a destination state? */
	int *d; /* the list of positions in the transition array where this
		   state is the destination state. 
		 MEMORY ALLOCATED in read_trans.c */
	fileset *f; /* file name and pointer for J/k distribution file(s) 
		     if there are any for this state 
		     MEMORY ALLOCATED in read_trans.c */
	double *T; /* Temperature of state.  Leave zero if this is only a
		       destination state or if temp is defined at a 
		       higher level.  There is one of these for each omega. 
		    MEMORY ALLOCATED in read_trans.c */
	} State;

/*** The `caseinfo' structures ***/
/* For now, only cases (a) and (b) are defined here. */ 

typedef struct {  
	int cTe, cwe, cwexe, cBe, cae, cbeta, cO; /* change flags for the 
		term energy (cTe), we and wexe (cwe, cwexe), the rotational 
		constant (cBe), the vibration-rotation interaction constant 
		(cae), the centrifugal distortion constant (cbeta), anything
	       	concerning Omegas (cO) */
	double Te, we, wexe, Be, ae, beta;  /* all units wavenumbers -- beta
					currently not used in the program */ 
	int L, O, g; /* Lambda (L).  The number of Omegas
	   to consider; use 0 (zero) to have the program use all 
	   possible.  If applicable (g), gerade (1), ungerade (-1); if 
	   not homonuclear, use 0 (zero). */ 
	double *lO, I, S, sflag;  /* if the user specifies certain Omegas to 
	   use, this is the list of desired Omegas (lO).  It has to be float 
	   since I don't want to work out some other way to have half-integral 
	   values for Omega.  Also, if homonuclear, the nuclear spin (I). 
	   The electronic spin, S, and a flag for looking up whether the spin
	   is integral (=0) or half-integral (=0.5).
	    MEMORY ALLOCATED in read_trans.c */
	double *pO; /* If the user specifies populations for each Omega, they 
	       go here.  MEMORY ALLOCATED in read_trans.c */
	int symm;  /* If homonuclear:  are even-plus (+1) or even-minus (-1) 
		 (even J-values) symmetric (default value 0)? -- This variable
		    will only be useful if someone includes the +/- splitting
		    in the program one day -- SEE PROGRAM_CHANGES if you are
		   the person wanting to do that */ 
	double Jmin, Jmax; /* Minimum and maximum values for J. */
	double *Jpop;  /* pointer to the array of J-state populations if they
		are entered from a file MEMORY ALLOCATED in read_rot.c */ 
	int a;  /* how many J-values long is the array? */
	} Case_a_stateinfo;

typedef struct {  
	int cTe, cwe, cwexe, cBe, cae, cbeta; /* change flags for the term 
		energy (cTe), we and wexe (cwe, cwexe), the rotational 
		constant (cBe), the vibration-rotation interaction constant 
		(cae), the centrifugal distortion constant (cbeta) */
	double Te, we, wexe, Be, ae, beta;  /* all units wavenumbers -- beta
	    is currently included in some calculations in the program, but 
	    since it isn't read in, the value is always zero. */ 
	int L, p, g; /* Lambda (L).  If L=0, indicate (p) plus (+1) or minus 
	    (-1); if L!=0, use zero (0).  If homonuclear, indicate (g) gerade 
	    (+1), ungerade (-1); if not homonuclear, use zero (0). */ 
	double I, S, sflag; /* If homonuclear, the nuclear spin.  The 
	   electronic spin, S and a flag for looking up whether the spin is
	   integral (=0) or half-integral (=0.5). */
	/*int symm;  [1] If L=0 and homonuclear:  are even (+1) or odd (1) 
		K-values symmetric?  [2] If Lambda>0 and homonuclear are 
		even-plus states (+1) or even-minus states (-1) symmetric?  
		(default 0) */ 
	double Kmin, Kmax, Jnum, *K, *J;  /* Minimum and maximum values for K,
	     the number of J's, and the actual values of J and K */
	double *Jpop;  /* pointer to the array of J-state populations if they
	  are entered from a file MEMORY ALLOCATED in read_rot.c */
	int a;  /* how many K-values long is the array? */
	} Case_b_stateinfo;


/****************** Begin Section on Transitions ******************/

/*** General structures for simulation parts ***/ 
typedef struct{
	int n; /* number of transitions here */
	double pop; /* relative population for this set */ 
	double *f; /* frequency of transition */
	double *ni; /* native intensity at that transition */
	double *ci; /* cascade intensity at that transition if needed */ 
	} simset;  /* set of information needed to make the simulation */ 
typedef struct{
	int n; /* number of points */
	double *x; /* position */
	double *y; /* relative intensity */
	} sim;  /* a simulation */

/*** A structure for Franck-Condon Factors ***/ 
typedef struct { 
	double fcfn[30]; /* An array of frequency-adjusted franck-condon 
			    factors from one vibrational level to a series 
			    of lower levels, for observed transitions. */
	double fcfc[30]; /* An array of frequency-adjusted franck-condon 
			    factors from one vibrational level to a series of 
			    lower levels for transitions that contribute to
			    cascade transitions. */
	double cm[30]; /* array of cm's for the origin of this transition */
	int vqn[30];  /* redundant, but just in case the vib quantum number 
			 isn't what's expected from the file -- these are the
		       lower state quantum numbers. */
	} vfcf;

/*** The `Transition' structure ***/ 
typedef struct {
	char Nhi[21], Nlo[21]; /* The letter designations for the 
		higher state (Nhi -- e.g., a', B, c, etc.), and the 
		lower state (Nlo). */
	int Hi, Lo; /* Numerical position in states array for info 
			     about the high state and the low state. */ 
	int Ohi,Olo; /* Number of Omegas specified for the high and low
			                        states. */
	int vs,vnh,vnl,vhlo,vllo; /* the actual number of v's found in the 
		file, the number of high and low v's simulated and the lowest
		values for those sets of v's */
	int cv,*cP;  /* change flag for FCF info and transition probability 
		      MEMORY ALLOCATED in read_trans.c */
	vfcf *v; /* pointer to the arrays of franck-condon factors.  
		       The declaration of 30 states is based on the maximum 
		       number output from Ervin's FCF program.  These array
		     spots are indexed by the upper vibrational q#. 
		  MEMORY ALLOCATED in read_FCF.c*/
	int nS, *fS; /* number of simulations for this transition and flags
		      to tell if the simulations should be included in the
		      overall simulation (won't be if spin-Sigma forbidden) 
		      MEMORY ALLOCATED in transition functions (e.g., a-a.c)*/ 
	simset *S; /* pointer to simulations for this transition 
		    MEMORY ALLOCATED in transition functions (e.g., a-a.c)*/ 
	double *P; /* overall electronic transition probability 
		    MEMORY ALLOCATED in read_trans.c */
	fileset *f; /* filenames and pointers for FCF files 
		     MEMORY ALLOCATED in read_trans.c */
	} Trans;

/* This is a structure made solely for the purpose of transmitting
   information between b_b() and the various sub-functions it calls */ 
typedef struct {
	int hvnlo,hvnnum,hvclo,hvcnum,hivlo,hivnum; /* Re high state: native
		low vib num, number vib levels; same, but for cascade; overall
		lowest vib num and overall number vib levels */
	int lvnlo,lvnnum,lvclo,lvcnum,lovlo,lovnum; /* like above, low state */
	int nT; /* maximum number transitions per K level */
	int nhJ,nlJ; /* 2S+1 for high and low states */
	int hiK; /* maximum number hi K might be */
	int m; /* position in MOL array */
	Trans *T; /* pointer to this transition structure */
	Case_b_stateinfo *Ch,*Cl; /* to high and low Case_b_info structures */
	rotset *Rh,*Rl; /* to high and low rotsets */ 
	State *Sh,*Sl; /* to high and low State structures */
	vset *Vh; /* high-state vib set */
	} BBTinfo;

/*** Included Atomic Transitions ***/ 
typedef struct {
	double p;  /* relative population for this file */
	fileset f; /* name and pointer for this file */
	double *x;  /* position of the transition */
	double *A;  /* Einstein A */
	double *pop; /* relative population for this transition */
	int M; /* multiplicity for the origination state -- not used now */
	simset S;  /* the part of the simulation due to this set */
	int c, n; /* change flag for this set of atomic transitions,
			number of sets of information */
	} Ainfo;  

/*** Included Transitions, `Other' `experimental', etc. ***/ 
typedef struct {
	double p; /* relative population for this set */ 
	double *x; /* position of transition */
	double *y; /* relative intensity at that position */
	int c, n; /* change flag for this set of data and number of
		   x,y pairs in this set. */
	fileset f; /* file where this info was found */
	} XYlist;

/* Structure for keeping up with molecule information */ 
typedef struct{
	int cp, cT; /* change flags for relative population and temperature */
	char Mol[41]; /* name of molecule */
	double pop; /* relative population for molecule */
	int states ; /* number of states for this molecule */ 
	State *s; /* state array 
		   MEMORY ALLOCATED in read_trans.c */
	int trans; /* number of transitions for this molecule */
	Trans *t; /* transition array 
		   MEMORY ALLOCATED in read_trans.c */
	double T; /* temperature if specified at molecule level */
	} Molecule;

/* these are the names for files containing information about the run */
char PREF[101],TMPFILE[300];
FILE *INMAIN,*DBG,*OUT,*PAR,*SCR,*SYS;
fileset *INOT,*INV,*OTR,*INEX,*INAT,INST,INTR;

/* these functions read in the various input files */
void read_state_and_transition_files(),read_FCF_file(),read_atomic_files();
void read_XY_file(XYlist*),read_rot_dist_files();
/* these are managerial and graphical functions */
/*void get_energies_intensities();*/
void manage_transitions(),manage_rotations();
void do_simulation(), graph_results(), ask_for_changes();
/* these functions calculate spectra for different transition types */
void a_a(int,int),a_b(int,int),a_c(int,int),a_d(int,int),b_b(int,int);
void b_c(int,int),b_d(int,int),c_c(int,int),c_d(int,int),d_d(int,int);
/* these functions calculate rotational levels and populations */
void singlet_JK(int,int),multiplet_zeroL_JK(int,int);
void multiplet_nonzeroL_JK_a(int,int),multiplet_nonzeroL_JK_b(int,int);
void case_a_cascade(rotset*,int,double,int,int,int);
void case_b_cascade(rotset*,int,int,int,int,int);
void bb_SigmaSigma(BBTinfo),bb_SigmaOther(BBTinfo),bb_Other(BBTinfo);
double HL(int,double,int,int);
char *new_mol_ck(char*,char*);

int DEBUG,NUMTR,NUMAT,NUMOT,NUMEX,NUMMOL,NUMST;
/* A note on the DEBUG variable.  If this is equal to -1, no functions
  will print any messages. Adjusting values in .rvesimconfig can cause
 messages about the execution of functions to print.  Use this feature 
 cautiously as it can produce enormous files. */ 
int Case_a_num,Case_b_num,SCANTYPE,DETECTTYPE,UNITTYPE;
int EFFTYPE,INTERACT,ROVIB;
double SIMMAX,*ATPOP,*OTINT,*EXINT,TEMP,JKCUT;
double MINV,MAXV,PSTEP,MAXINT,RESN;
Case_a_stateinfo *CA;
Case_b_stateinfo *CB;
Molecule *MOL;
Ainfo *AT;
XYlist *OT,*EX,*EF;
//char PROGRAM_NAME[200]="RVESIM";
/* debugging-level flags */
/* Functions that affect many functions */
int mdebugflag;/* main */
int tdebugflag; /* transition file */
/* Functions for calculating rotational populations and energies */
int mrdebugflag; /* manage rotations */
int mnadebugflag; /* multiplet non-zero-L Case a */
int mnbdebugflag;  /* multiplet non-zero-L Case b */
int msdebugflag; /* multiplet sigma */
int sdebugflag; /* singlet sigma */
/* Functions involved in calculating transitions */
int mtdebugflag; /* manage transitions */
int aadebugflag; /* Case a <-> case a transition */
int cacdebugflag; /* calculate energy levels in destination Case a state */
int bbdebugflag; /* manages case b <-> case b transitions */
int cbcdebugflag; /* calculate energy levels in destination Case b state */
int bbSSdebugflag; /* Case b Sigma<->Sigma transitions */
int bbSOdebugflag; /* Case b Sigma<->Non-Sigma transitions */
int bbOOdebugflag; /* Case b Non-Sigma<->Non-Sigma transitions */
/* The function that turns lines into spectra */
int ddebugflag; /* do_simulation */
/* Functions that read contents of files */
int rpdebugflag; /* read FCF files */
int radebugflag; /* read atomic files */
int rddebugflag; /* read rotational distribution files */
int XYdebugflag; /* read files in x-y format */
int HLswitch; /* use Honl-London factors? 1=no 0=yes */
double  bbSSSATT,bbSOSATT,bbOOSATT;

char *PROGRAM_NAME;
