LIBNAME SASDATA "/folders/myfolders/Applied_Statistics/Data";

/* Spearman's Rho Estimations */
%MACRO SpearmanRho(rho=);
	PROC IML;
		pi = constant("pi");
		tau=&rho;
		
		start innerGum(y) global(alpha, x); 
		   return(Exp(-((-Log(x))**alpha + (-Log(y))**alpha)**(1/alpha)));
		finish; 
		
		start outerGum(par) global(x,alpha); 
			x=par;
		   yinterval = 0 || 1;
		   /** evaluate inner integral for the parameter value, a=x **/ 
		   call quad(w, "innerGum", yinterval);
		   return (w);
		finish; 
		
		start finalGum(param) global(alpha, tau);
			alpha=param;
			xinterval= {0 1};
			call quad(v, "outerGum", xinterval); /** outer integral **/ 
			return(12*v-(3+tau));
		finish;
		
		intervalsGum = {1 100};        
		SGum = froot("finalGum", intervalsGum);
		
		start innerClay(y) global(alpha, x); 
			return((x**(-alpha)+y**(-alpha)-1)**(-1/alpha));
		finish; 
		
		start outerClay(par) global(x, alpha); 
			x=par;
			
			if(alpha>0) then yinterval= 0||1;
			else yinterval= (1-x**(-alpha))**(-1/alpha)||1;
		   	/** evaluate inner integral for the parameter value, a=x **/ 
		    call quad(w, "innerClay", yinterval);
		   	return (w);
		finish; 
		
		start finalClay(param) global(alpha, tau);
			alpha=param;
			xinterval= {0 1};
			call quad(v, "outerClay", xinterval); /** outer integral **/ 
			return(12*v-(3+tau));
		finish;
		                 
		intervalsClay = {-1 10};        
		SClay = froot("finalClay", intervalsClay);
		SGau=2*sin(pi*tau/6);
		SFGM=3*tau;
		
		start innerFrk(y) global(alpha, x); 
			return(-(1/alpha)*Log(1+(Exp(-alpha*x)-1)*(Exp(-alpha*y)-1)/(Exp(-alpha)-1)));
		finish; 
		
		start outerFrk(par) global(x, alpha); 
		x=par;
		   yinterval = 0 || 1;
		   /** evaluate inner integral for the parameter value, a=x **/ 
		   call quad(w, "innerFrk", yinterval);
		   return (w);
		finish; 
		
		start finalFrk(param) global(alpha, tau);
			alpha=param;
			xinterval= {0 1};
			call quad(v, "outerFrk", xinterval); /** outer integral **/ 
			return(12*v-(3+tau));
		finish;
         
		intervalsFrk = {-30 30};        
		SFrk = froot("finalFrk", intervalsFrk);
		
		CPAR=SGum||SClay||SFrk||SGau||SFGM;
		create EstSpearman from CPAR [colname={'Gumbel alpha','Clayton alpha','Frank alpha','Gaussian rho','FGM alpha'}]; 
			append from CPAR;       
		close EstSpearman;
	QUIT;
	
	PROC PRINT data=EstSpearman noobs;
		title "Spearman's rho estimations";
	RUN;
%MEND SpearmanRho;

/* Kendall's Tau Estimations */
%MACRO KendallTau(tau=);
	PROC IML;
		pi = constant("pi");
		tau=&tau;
		
		SGum=1/(1-tau);
		SClay=2*tau/(1-tau);
		SGau=sin(pi*tau/2);
		SFGM=9*(tau/2);
		
		start D(y);
			return(y/(Exp(y)-1));
		finish;
		
		*IF alpha>0 / tau>0;
		start FC(x) global(tau);
			dinterval=0||x;
			call quad(w, "D", dinterval);
			return(1-(4/x)*(1-(1/x)*w)-tau);
		finish;
		
		intervals = {.00001 20};        
		SFrk = froot("FC", intervals);
		
		CPAR=SGum||SClay||SFrk||SGau||SFGM;
		create EstKendall from CPAR [colname={'Gumbel alpha','Clayton alpha','Frank alpha','Gaussian rho','FGM alpha'}]; 
			append from CPAR;       
		close EstKendall;
	QUIT;
	
	PROC PRINT data=EstKendall noobs;
		title "Kendall's tau estimations";
	RUN;
%MEND KendallTau;

/* Gumbel's Copula Simulation */
%MACRO SIM_Gum_Uniform(alpha=, nsim=, seed=); 
	PROC IML;
		call streaminit(&seed);
		alpha=&alpha;
		
		do i=1 to &nsim by 1; 
			U1=rand('Uniform'); 
			U2=rand('Uniform');
			
			start Func(x) global(U1, U2, alpha); 
				return(
					Exp(-((-Log(x))**alpha + (-Log(U1))**alpha)**(1/alpha)) *
					((-Log(x))**alpha + (-Log(U1))**alpha)**(-1 + 1/alpha) *
					((-Log(U1))**(alpha -1))/U1-U2
				);
			finish;
			
			intervals = {.00001 1};
			U2C = froot("Func", intervals);
			X = X // U1; 
			Y = Y // U2C; 
			YI = YI // U2; 
		end;
		
		Total = X || Y || YI;
		create GumC from Total [colname={'X','Y','YI'}]; 
			append from Total;
		close GumC;
	QUIT;
%MEND SIM_Gum;
%MACRO SIM_Gum(alpha=, nsim=, seed=, dataset=, uvar=);
	PROC IML;
		call streaminit(&seed);
		alpha=&alpha;
		
		use &dataset;
			read all var{&uvar};
		close &dataset;
		U=&uvar;
		
		do i=1 to &nsim by 1;
			U1=U[i];
			U2=rand('Uniform');
			
			start Func(x) global(U1,U2,alpha);
				return(Exp(-((-Log(x))**alpha + (-Log(U1))**alpha)**(1/alpha))*((-Log(x))**alpha + (-Log(U1))**alpha)**(-1 + 1/alpha)*((-Log(U1))**(alpha-1))/U1-U2);
			finish;
			
			intervals = {0.00001 1};        
			U2C = froot("Func", intervals);
			
			X=X//U1;
			Y=Y//U2C;
			YI=YI//U2;
		end;
		
		Total=X||Y||YI;
		create GumC from Total [colname={'X','Y','YI'}]; 
			append from Total;       
		close GumC;
	QUIT;
%MEND SIM_Gum;

/* Gaussian Copula Simulation */
%MACRO SIM_GC(rho=, nsim=, seed=, dataset=, uvar=);
	PROC IML;
		call streaminit(&seed);
		rho=&rho;
		
		use &dataset;
			read all var{&uvar};
		close &dataset;
		U=&uvar;
		
		do i=1 to &nsim by 1;
			U1=U[i];
			U2=rand('Uniform');
			
			start Func(x) global(U1,U2,rho);
				return(CDF('Normal',quantile('NORMAL', x),rho*quantile('NORMAL',U1),(1-rho**2))-U2);
			finish;
			
			intervals = {0.00001 0.99999};        
			U2C = froot("Func", intervals);
			
			X=X//U1;
			Y=Y//U2C;
			YI=YI//U2;
		end;
		
		Total=X||Y||YI;
		create GC from Total [colname={'X','Y','YI'}]; 
			append from Total;       
		close GC;
	QUIT;
%MEND SIM_GC;
%MACRO SIM_GC_Uniform(rho=, nsim=, seed=);
	PROC IML;
		use marginals;
			read all var{U_IMP4};
		close marginals;
		
		U=U_IMP4;
		
		call streaminit(&seed);
		rho=&rho;
		
		do i=1 to &nsim by 1;
			U1=U[i];
			U1=rand('Uniform');
			U2=rand('Uniform');
			
			start Func(x) global(U1,U2,rho);
				return(
					CDF(
						'Normal',
						quantile('NORMAL', x),
						rho * quantile('NORMAL', U1),
						(1 - rho**2)
					) - U2
				);
			finish;
			
			intervals = {.00001 .99999};        
			U2C = froot("Func", intervals);
			
			X=X//U1;
			Y=Y//U2C;
			YI=YI//U2;
		end;
		
		Total=X||Y||YI;
		create GC from Total [colname={'X','Y','YI'}]; 
			append from Total;       
		close GC;
	QUIT;
%MEND SIM_GC;

/* Clayton's Copula Simulation */
%MACRO SIM_Clay(alpha=, nsim=, seed=, dataset=, uvar=);
	PROC IML;
		call streaminit(&seed);
		alpha=&alpha;
		
		use &dataset;
			read all var{&uvar};
		close &dataset;
		U=&uvar;
		
		do i=1 to &nsim by 1;
			U1=U[i];
			U2=rand('Uniform');
			
			start Func(x) global(U1,U2,alpha);
				return(U1**(-1 -alpha)*(x**(-alpha) + U1**(-alpha)-1)**(-1 - 1/alpha)-U2);
			finish;
			
			intervals = {0.001 1};        
			U2C = froot("Func", intervals);
			
			X=X//U1;
			Y=Y//U2C;
		end;
		
		Total=X||Y;
		create CC from Total [colname={'X','Y'}]; 
			append from Total;       
		close CC;
	QUIT;
%MEND SIM_Clay;
%MACRO SIM_Clay_Uniform(alpha=, nsim=, seed=);
	PROC IML;
		call streaminit(&seed);
		alpha=&alpha;
		
		do i=1 to &nsim by 1;
			U1=rand('Uniform');
			U2=rand('Uniform');
			
			start Func(x) global(U1,U2,alpha);
				return(
					U1**(-1 -alpha) *
					(x**(-alpha) + U1**(-alpha)-1)**(-1 - 1/alpha)
					-U2
				);
			finish;
			
			intervals = {0.001 1};        
			U2C = froot("Func", intervals);
			
			X=X//U1;
			Y=Y//U2C;
		end;
		
		Total=X||Y;
		create CC from Total [colname={'X','Y'}]; 
			append from Total;       
		close CC;
	QUIT;
%MEND SIM_Clay;

/* Frank's Copula Simulation */
%MACRO SIM_Frk(alpha=, nsim=, seed=, dataset=, uvar=);
	PROC IML;
		call streaminit(&seed);
		alpha=&alpha;
		
		use &dataset;
			read all var{&uvar};
		close &dataset;
		U=&uvar;
		
		do i=1 to &nsim by 1;
			U1=U[i];
			U2=rand('Uniform');
			
			start Func(x) global(U1,U2,alpha);
				return(
				(Exp(alpha)*(-1 + Exp(alpha*x))) /
				(-Exp(alpha) + Exp(alpha*(1+x)) - 
				Exp(alpha*(U1+x)) + Exp(alpha*(1 + U1))
				)-U2);
			finish;
			
			intervals = {0.00001 1};        
			U2C = froot("Func", intervals);
			
			X=X//U1;
			Y=Y//U2C;
		end;
		
		Total=X||Y;
		create FrkC from Total [colname={'X','Y'}]; 
			append from Total;       
		close FrkC;
	QUIT;
%MEND SIM_Frk;
%MACRO SIM_Frk_Uniform(alpha=, nsim=, seed=);
	PROC IML;
		call streaminit(&seed);
		alpha=&alpha;
		
		do i=1 to &nsim by 1;
			U1=rand('Uniform');
			U2=rand('Uniform');
			
			start Func(x) global(U1,U2,alpha);
				return(
					(Exp(alpha)*(-1 + Exp(alpha*x))) /
					(-Exp(alpha) + Exp(alpha*(1+x)) - 
						Exp(alpha*(U1+x)) + Exp(alpha*(1 + U1))
					)-U2
				);
			finish;
			
			intervals = {.00001 1};        
			U2C = froot("Func", intervals);
			
			X=X//U1;
			Y=Y//U2C;
		end;
		
		Total=X||Y;
		create FrkC from Total [colname={'X','Y'}]; 
			append from Total;       
		close FrkC;
	QUIT;
%MEND SIM_Frk;

/* FGM's Copula Simulation*/
%MACRO SIM_FGM(alpha=, nsim=, seed=, dataset=, uvar=);
	PROC IML;
		call streaminit(&seed);
		alpha=&alpha;
		
		use &dataset;
			read all var{&uvar};
		close &dataset;
		U=&uvar;
		
		do i=1 to &nsim by 1;
			U1=U[i];
			U2=rand('Uniform');
			
			start Func(x) global(U1,U2,alpha);
				return(x*(1 + alpha*(1 - x)*(1 - U1)) - alpha*(1 - x)*x*U1-U2);
			finish;
			
			intervals = {0.00001 1};        
			U2C = froot("Func", intervals);
			
			X=X//U1;
			Y=Y//U2C;
			YI=YI//U2;
		end;
		
		Total=X||Y||YI;
		create FGMC from Total [colname={'X','Y','YI'}]; 
			append from Total;       
		close FGMC;
	QUIT;
%MEND SIM_FGM;
%MACRO SIM_FGM_Uniform(alpha=, nsim=, seed=);
	PROC IML;
		call streaminit(&seed);
		alpha=&alpha;
		
		do i=1 to &nsim by 1;
			U1=rand('Uniform');
			U2=rand('Uniform');
			
			start Func(x) global(U1,U2,alpha);
				return(x*(1 + alpha*(1 - x)*(1 - U1)) - alpha*(1 - x)*x*U1-U2);
			finish;
			
			intervals = {.00001 1};        
			U2C = froot("Func", intervals);
			
			X=X//U1;
			Y=Y//U2C;
			YI=YI//U2;
		end;
		
		Total=X||Y||YI;
		create FGMC from Total [colname={'X','Y','YI'}]; 
			append from Total;       
		close FGMC;
	QUIT;
%MEND SIM_FGM;


/* QUESTION 3.1 */
%SIM_Gum_uniform(nsim=1000, alpha=5, seed=12345);

/* Question 3.1.A */
PROC CORR data=GumC spearman kendall;
	var X Y;
RUN;
/* For both tests the null hypothesis (X and Y ar independent) */
/* is rejected with a p-value of <.0001. */

 
/* Question 3.1.B */
PROC IML;
	tau = 1 - 1 / 5;
	print tau;
QUIT;


/* Question 3.1.C & 3.1.D */
%KendallTau(tau=.81311);
%SpearmanRho(rho=.95131);


/* Question 3.1.E */
%SIM_Gum_Uniform(nsim=1000, alpha=5, seed=12345);
PROC SGPLOT data=Gumc aspect=1;
	title "Original simulated data";
	scatter x=X y=Y / markerattrs=(symbol=circlefilled color='darkgreen' size=8 );
RUN;
%SIM_Gum_Uniform(nsim=1000, alpha=5.38, seed=6789);
PROC SGPLOT data=Gumc aspect=1;
	title "Estimated Gumbal copula";
	scatter x=X y=Y / markerattrs=(symbol=circlefilled color='darkgreen' size=8 );
RUN;
%SIM_Frk_Uniform(nsim=1000, alpha=19, seed=6789);
PROC SGPLOT data=Frkc aspect=1;
	title "Estimated Frank copula";
	scatter x=X y=Y / markerattrs=(symbol=circlefilled color='cyan' size=8);
RUN;


/* QUESTION 3.2 */
DATA Q3_2;
	SET SASDATA.IVF;
	IMP = IMP + (ranuni(1)-0.5);
RUN;

/* Question 3.2.A */
PROC TRANSPOSE out=Q3_2_wide(drop = _NAME_ _LABEL_)
			   data=Q3_2 prefix=IMP;
	by ID;
	id PER;
	var IMP;
RUN;
/* Remove missing data */
DATA Q3_2_wide;
	set Q3_2_wide;
	if cmiss(of _all_) then delete;
RUN;
PROC RANK data=Q3_2_wide out=Q3_2_rank;
	var IMP4 IMP18;
	ranks rank_IMP4 rank_IMP18;
RUN;
PROC MEANS data=Q3_2_rank n;
	var rank_IMP4 rank_IMP18;
RUN;


/* Question 3.2.B */
PROC CORR data=Q3_2_wide kendall spearman;
	var IMP4 IMP18;
RUN;
/* Neither test rejects the null hypothesis  */
/* (X and Y are independent). */


/* Question 3.2.C */

/* Original data */
DATA Q3_2_marginals;
	set Q3_2_rank;
	U_IMP4=rank_IMP4/237; /* or 236 if errors*/
	U_IMP18=rank_IMP18/237;
RUN;
PROC SGPLOT data=Q3_2_marginals aspect=1;
	title "Marginals IMP4 and IMP18";
	scatter x=U_IMP4 y=U_IMP18 / markerattrs=(symbol=circlefilled color='red' size=8 );
RUN;

/* Estimations using Kendall's tau: */
%KendallTau(tau=.05106);
%SIM_GC(nsim=236, rho=.080119, seed=12345, dataset=Q3_2_marginals, uvar=U_IMP4);
%SIM_Gum(nsim=236, alpha=1.05381, seed=12345, dataset=Q3_2_marginals, uvar=U_IMP4);
%SIM_Clay(nsim=236, alpha=.10761, seed=12345, dataset=Q3_2_marginals, uvar=U_IMP4);
%SIM_Frk(nsim=236, alpha=.46051, seed=12345, dataset=Q3_2_marginals, uvar=U_IMP4);
%SIM_FGM(nsim=236, alpha=.22977, seed=12345, dataset=Q3_2_marginals, uvar=U_IMP4);
PROC SGPLOT data=Gumc aspect=1;
	title "Simulated Gumbel copula";
	scatter x=X y=Y / markerattrs=(symbol=circlefilled color='darkgreen' size=8 );
RUN;
PROC SGPLOT data=Frkc aspect=1;
	title "Simulated Frank's copula";
	scatter x=X y=Y / markerattrs=(symbol=circlefilled color='cyan' size=8 );
RUN;
PROC SGPLOT data=Gc aspect=1;
	title "Simulated Gaussian copula";
	scatter x=X y=Y / markerattrs=(symbol=circlefilled color='blue' size=8 );
RUN;
PROC SGPLOT data=Cc aspect=1;
	title "Simulated Clayton copula";
	scatter x=X y=Y / markerattrs=(symbol=circlefilled color='navy' size=8 );
RUN;
PROC SGPLOT data=Fgmc aspect=1;
	title "Simulated FGM copula";
	scatter x=X y=Y / markerattrs=(symbol=circlefilled color='black' size=8 );
RUN;

/* Estimations using spearman's rho: */
%SpearmanRho(rho=.03411);
%SIM_GC(nsim=236, rho=.035718, seed=12345, dataset=Q3_2_marginals, uvar=U_IMP4);
%SIM_Gum(nsim=236, alpha=1.02331, seed=12345, dataset=Q3_2_marginals, uvar=U_IMP4);
%SIM_Clay(nsim=236, alpha=.046550, seed=12345, dataset=Q3_2_marginals, uvar=U_IMP4);
%SIM_Frk(nsim=236, alpha=.20477, seed=12345, dataset=Q3_2_marginals, uvar=U_IMP4);
%SIM_FGM(nsim=236, alpha=.10233, seed=12345, dataset=Q3_2_marginals, uvar=U_IMP4);
PROC SGPLOT data=Gumc aspect=1;
	title "Simulated Gumbel copula";
	scatter x=X y=Y / markerattrs=(symbol=circlefilled color='darkgreen' size=8 );
RUN;
PROC SGPLOT data=Frkc aspect=1;
	title "Simulated Frank's copula";
	scatter x=X y=Y / markerattrs=(symbol=circlefilled color='cyan' size=8 );
RUN;
PROC SGPLOT data=Gc aspect=1;
	title "Simulated Gaussian copula";
	scatter x=X y=Y / markerattrs=(symbol=circlefilled color='blue' size=8 );
RUN;
PROC SGPLOT data=Cc aspect=1;
	title "Simulated Clayton copula";
	scatter x=X y=Y / markerattrs=(symbol=circlefilled color='navy' size=8 );
RUN;
PROC SGPLOT data=Fgmc aspect=1;
	title "Simulated FGM copula";
	scatter x=X y=Y / markerattrs=(symbol=circlefilled color='black' size=8 );
RUN;


/* QUESTION 3.3 */
/* TODO */
DATA Q3_3;
	do i = 1 to 1000;
	   X = rand('Uniform');
	   Y = rand('Uniform');
	   output;
	end;
RUN;

/* Question 3.3.A */
/* TODO */
/* 0.7 * 0.7 = 0.49 */
PROC CORR data=Q3_3 kendall spearman 
	plots=scatter(ellipse=none); 
	var X Y;
RUN;


/* QUESTION 3.4 */
DATA RCT;
	set SASDATA.RCT;
/* 	RESP = RESP + (ranuni(1)-0.5); */
RUN;

/* Question 3.4.A */
PROC TRANSPOSE out=Q3_4_wide(drop = _NAME_ _LABEL_)
			   data=RCT prefix=RESP;
	by ID;
	id TIME;
	var RESP;
RUN;


/* Question 3.4.B */
PROC CORR data=Q3_4_wide pearson kendall spearman;
	var RESP1 RESP2;
RUN;


/* Question 3.4.D */
PROC CORR data=Q3_4_wide pearson fisher(biasadj=no);
	var RESP1 RESP2;
RUN;


/* Question 3.4.C & 3.4.E */
%KendallTau(tau=.35941);
%SpearmanRho(rho=.50477);
/* The estimated alpha values for the FGM copula above the */
/* supported range [-1,1], so the FGM copula is not appropriate. */


/* Question 3.4.F */
PROC IML;
/* 	τ = α / (2+α) */
/* 	α = (2t) / (1-t) */
	estimate_α = 1+1;
	print estimate_α;
run;
/* Estimate for α = 1.12212 */


/* Question 3.4.G */
%SpearmanRho(rho=.50477);
/* Estimate for α = 1.09357 */
/* Close to the previous estimate, so we cannot rule out */
/* that the Clayton's copula appropriately models the */
/* dependency. */


/* Question 3.4.H */
DATA Q3_4_wide_ranuni;
	set Q3_4_wide;
	RESP1 = RESP1 + 0.1*(ranuni(1)-0.5);
	RESP2 = RESP2 + 0.1*(ranuni(1)-0.5);
RUN;
PROC RANK data=Q3_4_wide_ranuni out=Q3_4_ranked;
	var RESP1 RESP2;
	ranks rank_RESP1 rank_RESP2;
RUN;
PROC MEANS data=Q3_4_ranked n;
	var rank_RESP1 rank_RESP1;
RUN;
DATA Q3_4_marginals;
	set Q3_4_ranked;
	U_RESP1=rank_RESP1/716; /* or 236 if errors */
	U_RESP2=rank_RESP2/716;
RUN;
PROC SGPLOT data=Q3_4_marginals aspect=1;
	title "Original RESP1 and RESP2";
	scatter x=U_RESP1 y=U_RESP2 / markerattrs=(symbol=circlefilled color='red' size=8 );
RUN;
/* %SIM_Clay(nsim=700, alpha=1.12212, seed=12345); */
%SIM_Clay(nsim=700, alpha=1.12212, seed=12345, dataset=Q3_4_marginals, uvar=U_RESP1);
%SIM_Frk(nsim=700, alpha=3.62642, seed=12345, dataset=Q3_4_marginals, uvar=U_RESP1);
PROC SGPLOT data=Frkc aspect=1;
	title "Simulated Frank's copula";
	scatter x=X y=Y / markerattrs=(symbol=circlefilled color='cyan' size=8 );
RUN;
PROC SGPLOT data=Cc aspect=1;
	title "Simulated Clayton copula";
	scatter x=X y=Y / markerattrs=(symbol=circlefilled color='navy' size=8 );
RUN;
/* Using a little uniform noise it is clear that Frank's copula */
/* describes the dependency sturture better than Clayton's. */


/* QUESTION 3.6 */
DATA Q3_6;
	set SASDATA.RCT;
	WHERE center = 1 and (time = 5 or time = 6);
RUN;


/* Question 3.6.A */
PROC TRANSPOSE out=Q3_6_wide(drop = _NAME_ _LABEL_) 
		data=Q3_6 prefix=RESP;
	by ID;
	id TIME;
	var RESP; 
RUN;
DATA Q3_6_wide;
	set Q3_6_wide;
	DIFF = RESP6 - RESP5;
	RATIO = RESP6/RESP5;
/* 	RATIO_0 = RATIO-1; */
	LDIFF = log(RATIO);
RUN;
PROC UNIVARIATE data=Q3_6_wide normaltest;
	var DIFF LDIFF RATIO;
	histogram /normal;
    ods select histogram;
RUN;
/* The shapes look symetrical and could even be */
/* normally distributed. */


/* Question 3.6.B */
PROC UNIVARIATE data=Q3_6_wide normal mu0=0 0 1;
	var DIFF LDIFF RATIO;
	ods select TestsForLocation;
RUN;
/* Cannot reject null hypotheses for all cases  */
/* in all tests. */


/* Question 3.6.C */
/* For the paired t-test, we assume normality.  */
/* Normality does not seem apparent for the ratio. */

/* For the sign test, no assumptions are made. */

/* For the Wilcoxon signed rank test, we need symmetric  */
/* distributions. The his- tograms suggest slightly  */
/* left-skewed distributions, but nothing extreme. */


/* Question 3.6.D */
/* T-test; Wilcoxon; Sign-test */
/* based on the power */


/* Question 3.6.E */
/* Only sign-test is directly appropriate to test  */
/* equal medians. */
/* The Wilcoxon-signed rank test cannot be used since */
/* the assumption of symmetric differences is violated. */
/* For the T-test it is possible to rely on the 
/* CLT to test equal means. */


/* QUESTION 3.7 */
DATA Q3_7;
	set SASDATA.IVF;
	where PER = 4 or PER = 18;
RUN;


/* Question 3.7.A */
PROC TRANSPOSE out=Q3_7_wide(drop = _NAME_ _LABEL_) 
		data=Q3_7 prefix=IMP;
	by ID;
	id PER;
	var IMP; 
RUN;
DATA Q3_7_wide;
	set Q3_7_wide;
	DIFF = IMP4 - IMP18;
	RATIO = IMP4/IMP18;
	RATIO_0 = RATIO-1;
	LDIFF = log(RATIO);
RUN;
PROC UNIVARIATE data=Q3_7_wide normaltest;
	var DIFF LDIFF RATIO;
	histogram /normal;
    ods select histogram;
RUN;
/* LDIFF and RATIO seem to approximate a */
/* normal distribution. */
/* DIFF does not approximate a normal */
/* distribution and is not symmetric. */


/* Question 3.7.B & 3.7.C */
PROC UNIVARIATE data=Q3_7_wide normal mu0=0 0 1;
	var DIFF LDIFF RATIO;
	ods select TestsForLocation;
RUN;
/* For DIFF only the Sign-test is reliable which  */
/* rejects that the medians are equal */
/* For LDIFF and RATIO all 3 tests are appropriate */
/* and indicate that the medians and means are not equal */


/* QUESTION3.8 */
DATA COAG; 
	input Patient C K@@;
	datalines;
		1 120 132 8 145 133 15 117 123
		2 114 116 9 120 123 16 125 108
		3 129 135 10 129 116 17 136 131 
		4 128 115 11 126 127 18 151 119
		5 155 134 12 136 140 19 130 129 
		6 105 56 13 135 140 20 136 124
		7 114 114 14 125 114 21 113 112
	;
RUN;


/* Question 3.8.A */
PROC CORR data=COAG pearson spearman kendall;
	var C K;
RUN;
/* For all tests, reject the null hypotheses  */
/* H0: Rho=0 and Tau=0 */


/* Question 3.8.B */
%KendallTau(tau=.51346);
%SpearmanRho(rho=.66134);
PROC RANK data=COAG out=COAG_rank;
	var C K;
	ranks rank_C rank_K;
RUN;
PROC MEANS data=COAG_rank n;
	var rank_C rank_K;
RUN;
DATA COAG_marginals;
	set COAG_rank;
	U_C=rank_C/21;
	U_K=rank_K/21;
RUN;
PROC SGPLOT data=COAG_marginals aspect=1;
	title "Marginals C and K";
	scatter x=U_C y=U_K/ markerattrs=(symbol=circlefilled color='red' size=8 );
RUN;
%SIM_Gum(nsim=21, alpha=2.0553295, seed=12345, dataset=COAG_marginals, uvar=U_C);
/* %SIM_Gum(nsim=21, alpha=1.929525, seed=12345, dataset=COAG_marginals, uvar=U_C); */
%SIM_FGM(nsim=21, alpha=2.31057, seed=12345, dataset=COAG_marginals, uvar=U_C);
/* %SIM_FGM(nsim=21, alpha=1.98402, seed=12345, dataset=COAG_marginals, uvar=U_C); */
PROC SGPLOT data=Gumc aspect=1;
	title "Simulated Gumbel copula";
	scatter x=X y=Y / markerattrs=(symbol=circlefilled color='darkgreen' size=8 );
RUN;
PROC SGPLOT data=Fgmc aspect=1;
	title "Simulated FGM copula";
	scatter x=X y=Y / markerattrs=(symbol=circlefilled color='black' size=8 );
RUN;
/* Very few observations, but the model approximating the */
/* FGM copula using Kendall's estimate seems to be the  */
/* closest to the original data. */


/* Question 3.8.C */
DATA COAG;
	set COAG;
	DIFF = C - K;
	RATIO = C/K;
	RATIO_0 = RATIO-1;
	LDIFF = log(RATIO);
RUN;
PROC UNIVARIATE data=COAG normaltest;
	var DIFF LDIFF RATIO;
	histogram /normal;
    ods select histogram;
RUN;
PROC UNIVARIATE data=COAG normal mu0=0 0 1;
	var DIFF LDIFF RATIO;
	ods select TestsForLocation;
RUN;


/* Question 3.8.D */
DATA COAG;
	set COAG;
	TT_C = (C < 120);
	TT_K = (K < 120);
RUN;
PROC CORR data=COAG pearson;
	var TT_C TT_K;
RUN;
/* PROC FREQ data=COAG; */
/* 	tables TT_C * TT_K /CHISQ; */
/* RUN; */


/* Question 3.8.E */
PROC FREQ data=COAG;
	tables TT_C * TT_K /* /agree */; 
	exact mcnem;
RUN;


/* QUESTION 3.9 */
DATA Q3_9;
	set SASDATA.IVF;
	where PER = 10 or PER = 18;
RUN;
PROC TRANSPOSE out=Q3_9_wide(drop = _NAME_ _LABEL_) 
		data=Q3_9 prefix=IMP;
	by ID;
	id PER;
	var IMP; 
RUN;
DATA Q3_9_wide;
	set Q3_9_wide;
	if cmiss(of _all_) then delete;
	IMP10_85 = (IMP10 < 85);
	IMP18_85 = (IMP18 < 85);
	IMP_85 = (IMP10 < 85 or IMP18 < 85);
RUN;
DATA Q3_9_wide;
	set Q3_9_wide;
	DIFF = IMP10 - IMP18;
	RATIO = IMP10 / IMP18;
	RATIO_0 = RATIO - 1;
	LDIFF = log(RATIO);
	SIGN = (IMP10 >= IMP18);
RUN;


/* Question 3.9.A */
/* TODO */
PROC MEANS data=Q3_9_wide;
	var IMP10 IMP18 SIGN;
RUN;
PROC FREQ data=Q3_9_wide;
	tables SIGN;
RUN;
PROC UNIVARIATE data=Q3_9_wide normal mu0=0 0 1;
	var DIFF LDIFF RATIO;
	ods select TestsForLocation;
RUN;
PROC IML;
	LPL = cdf("Binom", 112, 0.5, 242);
	UPL = 1 - cdf("Binom", 130, 0.5, 242);
	print LPL UPL;
QUIT;
PROC FREQ data=Q3_9_wide;
	tables SIGN/binomial(P=0.5);
	exact binomial;
RUN;


/* Question 3.9.B & 3.9.C & 3.9.D*/
PROC FREQ data=Q3_9_wide;
	tables IMP10_85 * IMP18_85 /* /agree */;
	exact mcnem;
RUN;


/* QUESTION 3.10 */
%MACRO kappatest(nsim=,ntab=,n00=,n01=,n11=);
	PROC IML;
		G0 = (1-&n00-&n01-&n11);
		G = {&n00 &n01 &n11};
		call randseed(1234);
		C = {"SIMID", "N00", "N01", "N11", "N10"};
		PROB = G||G0;
		X = RandMultinomial(&nsim, &ntab, PROB);
		ID = 1:&nsim;
		
		X = ID‘||X;
		create CT from X[c=c];
			append from X;
		close CT;
	QUIT;
	DATA CT;
		set CT;
		EPO = 1 - N11/&ntab-N00/&ntab;
		EPE = 1 - (N00+N10)*(N00+N01)/(&ntab**2)
			    - (N01+N11)*(N10+N11)/(&ntab**2);
		EKAPPA = 1 - EPO/EPE;
		PO = 1 - &n11 - &n00;
		PE = 1 - (&n00+(1-&n00-&n01-&n11))*(&n01+&n00)
			   - (&n01+&n11)*((1-&n00-&n01-&n11)+&n11);
		KAPPA = 1 - PO/PE;
	RUN;
%MEND;


/* Question 3.10.A */
%kappatest(nsim=1000, ntab=100, n00=.05, n01=.01, n11=.3);
/* Simmulation not working? */
PROC UNIVARIATE data=CT;
	var KAPPA;
	histogram /normal;
RUN;


/* Question 3.10.B */
/* TODO */


/* Question 3.10.C */
/* TODO */
%kappatest(nsim=1000, ntab=100, n00=.03, n01=.01, n11=.3);
/* not K approx .23 */
data CT;
	set CT;
	MCN = (N01-N10)**2/(N01+N10);
	PMCN_NAPPROX = 1 - cdf(’CHISQ ’,MCN ,1);
	PMCN_EXACT = 2*min(cdf(’BINOM ’,N01 ,0.5,N01+N10),
	1-cdf(’BINOM ’,N01 ,0.5,N01+N10));
run;


/* QUESTION 3.11 */
/* Drug A = 0, Drug B = 1 */
DATA Q3_11;
	input treatment high count drug;
	datalines;
	1 1 12 0
	0 1 11 0
	1 0 7 0
	0 0 70 0
	1 1 10 1
	0 1 13 1
	1 0 4 1
	0 0 73 1
	;
RUN;

/* Question 3.11.A */
PROC FREQ data=Q3_11;
	where drug = 0;
	tables treatment * high /* /agree */; 
	weight count;
	exact mcnem;
RUN;


/* Question 3.11.C */
PROC FREQ data=Q3_11;
	where drug = 0;
	tables treatment * high /* /agree */; 
	weight count;
	exact chisq;
RUN;


/* Question 3.11.D */
PROC FREQ data=Q3_11;
	table treatment * high /chisq;
	weight count ;
	exact chisq;
RUN;


/* QUESTION 3.12 */







































