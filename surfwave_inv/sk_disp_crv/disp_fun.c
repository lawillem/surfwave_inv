//Bram Willemsen.
//bramwillemsen -at- gmail.com
//May 2016

#include <math.h>
#include <stdlib.h>
#include <stdio.h>

double eval_rayleigh_disp_fun(int N, double* alphas, double* betas, double* rhos, double* ds, double* norm_facts, double omega, double C){
	//evaluate dispersion function for given C and omega using FIG 11
	//remember that C is 0 indexed unlike fortran. So indexes are different
	//ds = THKNES

	double CSQ, XK;
	double TWOBSQ, BETNSQ, ALPNSQ;
	double GAMMA1, GAM1M1, UKNP, VKNP, WKNP, RKNP, SKNP;
	double EPSIL0;
	double EPS1, EPS2, EPS3, EPS4, EPS15;
	double ZETA1,  ZETA2,  ZETA3,  ZETA4,  ZETA5,  ZETA6,  ZETA7, ZETA8, ZETA9;
	double THKKM, ARGALM, RALPHM, PM, QM, SINPM, SINQM, ARGBTM, RBETAM;
	double EXPPPM, EXPMPM, EXPPQM, EXPMQM, UKN, VKN, RALPHN, RBETAN, RALRBT, EPSILN;
	double XKNP, ZKNP, KKNP, LKNP, MAXVAL;
	double FRAYL;
	double *THKNES, *ALPMSQ, *BETMSQ, *EPS0, *EPS00;
	int MM, M, L;

	//allocate arrays
	EPS0   =(double *) malloc((N-1)*sizeof(double));
	EPS00  =(double *) malloc((N-1)*sizeof(double));
	ALPMSQ =(double *) malloc((N-1)*sizeof(double));
	BETMSQ =(double *) malloc((N-1)*sizeof(double));

	//fill arrays
	for(MM=0;MM<N-1;MM++){ //evaluate interface MM
		EPS0[MM]   = rhos[MM+1] / rhos[MM];
		EPS00[MM]  = 2.0e0*(betas[MM]*betas[MM] - EPS0[MM]*betas[MM+1]*betas[MM+1]);
		ALPMSQ[MM] = alphas[MM]*alphas[MM];
		BETMSQ[MM] = betas[MM]*betas[MM];
	}

	//INITIALIZING, TOP PAGE 127
	TWOBSQ = 2.0e0*betas[0]*betas[0];
	CSQ    = C*C;
	XK     = omega/C;
	ALPNSQ = alphas[N-1]*alphas[N-1];
	BETNSQ = betas[N-1]*betas[N-1];
	THKNES = ds;
	EPSIL0 = pow(-1,(N-1))*rhos[0]*rhos[0]/(2*BETNSQ*ALPNSQ*rhos[N-1]*rhos[N-1]);


	//initialize quantities (66) for land model
	GAMMA1 = TWOBSQ/CSQ;
	GAM1M1 = GAMMA1-1.0e0;
	UKNP   = -GAMMA1*GAM1M1;
	VKNP   = 0.0e0;
	WKNP   = GAM1M1*GAM1M1;
	RKNP   = GAMMA1*GAMMA1;
	SKNP   = 0.0e0;  //land case
	M      = N;
	L      = M;
	
	//COMPUTE THE ELEMENTS OF THE LEFT-HAND MATRIX IN (67) USING (64) AND (65)
	for(MM=0;MM<N-1;MM++){ //evaluate interface MM
		EPS15  = -EPS0[MM];
		EPS1   = EPS00[MM]/CSQ;
		EPS2   = EPS1-1.0e0;
		EPS3   = EPS1-EPS15;
		EPS4   = EPS2-EPS15;
		THKKM  = THKNES[MM]*XK;
		ARGALM = 1.0e0 - CSQ/ALPMSQ[MM];

		if(ARGALM >= 0.0e0) goto one_ninety; //goto statements... the horror

		RALPHM = sqrt(-ARGALM);
		PM     = THKKM*RALPHM;
		SINPM  = sin(PM);
		ZETA1  = cos(PM);
		ZETA3  = RALPHM*SINPM;

		one_eighty:

		ARGBTM=1.0e0-CSQ/BETMSQ[MM];
		if(ARGBTM>=0.0e0) goto two_hundred;
		RBETAM=sqrt(-ARGBTM);
		QM    =THKKM*RBETAM;
		SINQM =sin(QM);
		ZETA2 =cos(QM);
		ZETA5 =RBETAM*SINQM;
		goto two_ten;

		one_ninety:

		RALPHM = -sqrt(ARGALM);
		EXPPPM = 0.5e0*exp(THKKM*RALPHM);
		EXPMPM = 0.25e0/EXPPPM;
		SINPM  = EXPPPM-EXPMPM;
		ZETA1  = EXPPPM+EXPMPM;
		ZETA3  = -RALPHM*SINPM;
		goto one_eighty;

		two_hundred:

		RBETAM = -sqrt(ARGBTM);
		EXPPQM = 0.5e0*exp(THKKM*RBETAM);
		EXPMQM = 0.25e0/EXPPQM;
		SINQM  = EXPPQM-EXPMQM;
		ZETA2  = EXPPQM+EXPMQM;
		ZETA5  = -RBETAM*SINQM;

		two_ten:

		ZETA4  = SINPM/RALPHM;
		ZETA6  = SINQM/RBETAM;
		ZETA7  = ZETA1*ZETA2;
		ZETA8  = ZETA1*ZETA5;
		ZETA9  = ZETA1*ZETA6;
		UKN    = 2.0e0*UKNP;
		VKN    = VKNP;

		if(2*(MM/2)!=MM) goto two_twenty; //If odd. In original code the conditional is for even, but MM is zero indexed in C instead of 1 indexed in Fortran. So even in Fortran is odd in C

		XKNP = ZETA4*(ZETA2*VKNP+ZETA6*WKNP)-ZETA7*RKNP+ZETA9*SKNP;
		ZKNP = ZETA8*VKNP-ZETA7*WKNP+ZETA3*(ZETA5*RKNP+ZETA2*SKNP);
		UKNP =-(EPS1*EPS4+EPS2*EPS3)*UKNP+EPS2*EPS4*XKNP+EPS1*EPS3*ZKNP;
		VKNP = EPS15*(ZETA4*(ZETA5*VKNP-ZETA2*WKNP)-ZETA8*RKNP-ZETA7*SKNP);
		SKNP = EPS15*(-ZETA7*VKN-ZETA9*WKNP-ZETA3*(ZETA2*RKNP-ZETA6*SKNP));
		WKNP = EPS2*(EPS2*XKNP-EPS1*UKN)+EPS1*EPS1*ZKNP;
		RKNP = EPS4*(EPS4*XKNP-EPS3*UKN)+EPS3*EPS3*ZKNP;

		goto two_thirty;

		two_twenty:

		KKNP = ZETA9*VKNP+ZETA7*WKNP-ZETA4*(ZETA6*RKNP-ZETA2*SKNP);
		LKNP = ZETA3*(ZETA2*VKNP-ZETA5*WKNP)+ZETA7*RKNP+ZETA8*SKNP;
		UKNP = -(EPS1*EPS4+EPS2*EPS3)*UKNP+EPS2*EPS4*KKNP+EPS1*EPS3*LKNP;
		VKNP = EPS15*(ZETA3*(ZETA6*VKNP+ZETA2*WKNP)+ZETA9*RKNP-ZETA7*SKNP);
		SKNP = EPS15*(-ZETA7*VKN+ZETA8*WKNP+ZETA4*(ZETA2*RKNP+ZETA5*SKNP));
		WKNP = EPS4*(-EPS4*KKNP+EPS3*UKN)-EPS3*EPS3*LKNP;
		RKNP = EPS2*(-EPS2*KKNP+EPS1*UKN)-EPS1*EPS1*LKNP;

		two_thirty:

		//p132 says, determine max of U,V,S,W and R and normalize everything by that.
		//Store normalization so that when we do quadratic root search later we can equalize the normalization
		//Use absolute value. Get by getting root of square.
		MAXVAL = sqrt(UKNP*UKNP);
		if (sqrt(VKNP*VKNP) > MAXVAL) MAXVAL = sqrt(VKNP*VKNP);
		if (sqrt(SKNP*SKNP) > MAXVAL) MAXVAL = sqrt(SKNP*SKNP);
		if (sqrt(WKNP*WKNP) > MAXVAL) MAXVAL = sqrt(WKNP*WKNP);
		if (sqrt(RKNP*RKNP) > MAXVAL) MAXVAL = sqrt(RKNP*RKNP);

		//printf("TEMPORARILY NOT USING MAXVAL=%e\n",MAXVAL);
		//MAXVAL = 1e0;

		//normalize
		norm_facts[MM] = MAXVAL; //store
		UKNP = UKNP/MAXVAL;
		VKNP = VKNP/MAXVAL;
		SKNP = SKNP/MAXVAL;
		WKNP = WKNP/MAXVAL;
		RKNP = RKNP/MAXVAL;
	}
	//done with loop
	if (L > N) goto two_sixty;

	//Perform the matrix multiplication of eq (68)
	//two_forty:

	RALPHN = -sqrt(1.0e0-CSQ/ALPNSQ);
	RBETAN = -sqrt(1.0e0-CSQ/BETNSQ);
	RALRBT = RALPHN*RBETAN;
	EPSILN = -EPSIL0*CSQ*CSQ/RALRBT;
	if (2*(N/2) != N) goto two_fifty; //If odd. In original code the conditional is for even, but MM is zero indexed in C instead of 1 indexed in Fortran. So even in Fortran is odd in C
	FRAYL  = EPSILN*(-VKNP*RBETAN+WKNP-RKNP*RALRBT-SKNP*RALPHN);
	goto four_forty;

	two_fifty:
	FRAYL  = EPSILN*(VKNP*RALPHN+RKNP-WKNP*RALRBT+SKNP*RBETAN);
	goto four_forty;

	//FRAYL IS THE RETURN VALUE

	two_sixty:
	//Here would go code from figure 12 and 13 for those geometries with bottom layers in the core

	four_forty:

	//prepare to exit function
	//free heap allocated arrays
	free(EPS0);
	free(EPS00);
	free(ALPMSQ);
	free(BETMSQ);

	return FRAYL;
}
