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
		
		Total=X||Y||YI
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
		
		Total=X||Y
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
/* is rejected with a p-value of <.0001 */

 
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
DATA IVF;
	SET SASDATA.IVF;
	IMP = IMP + (ranuni(1)-0.5);
RUN;

/* Question 3.2.A */
PROC TRANSPOSE out=WIDE_IVFN(drop = _NAME_ _LABEL_)
			   data=IVFN prefix=IMP;
	by ID;
	id PER;
	var IMP;
RUN;
/* Remove missing data */
DATA WIDE_IVFN;
	set WIDE_IVFN;
	if cmiss(of _all_) then delete;
RUN;
PROC RANK data=WIDE_IVFN out=ranked_a;
	var IMP4 IMP18;
	ranks rank_IMP4 rank_IMP18;
RUN;
PROC MEANS data=ranked_a n;
	var rank_IMP4 rank_IMP18;
RUN;


/* Question 3.2.B */
PROC CORR data=WIDE_IVFN kendall spearman;
	var IMP4 IMP18;
RUN;
/* Neither test rejects the null hypothesis (X and Y are independent) */

/* Question 3.2.C */

/* Original data */
DATA marginals;
	set ranked_a;
	U_IMP4=rank_IMP4/237; /* or 236 if errors*/
	U_IMP18=rank_IMP18/237;
RUN;
PROC SGPLOT data=marginals aspect=1;
	title "Marginals IMP4 and IMP18";
	scatter x=U_IMP4 y=U_IMP18 / markerattrs=(symbol=circlefilled color='red' size=8 );
RUN;

/* Estimations using Kendall's tau */
%KendallTau(tau=.05106);
%SIM_GC(nsim=236, rho=.080119, seed=12345);
%SIM_Gum(nsim=236, alpha=1.05381, seed=12345);
%SIM_Clay(nsim=236, alpha=.10761, seed=12345);
%SIM_Frk(nsim=236, alpha=.46051, seed=12345);
%SIM_FGM(nsim=236, alpha=.22977, seed=12345);
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

/* Estimations using spearman's rho */
%SpearmanRho(rho=.03411);
%SIM_GC(nsim=236, rho=.035718	, seed=12345);
%SIM_Gum(nsim=236, alpha=1.02331	, seed=12345);
%SIM_Clay(nsim=236, alpha=.046550	, seed=12345);
%SIM_Frk(nsim=236, alpha=.20477, seed=12345);
%SIM_FGM(nsim=236, alpha=.10233, seed=12345);
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
DATA WEEK3_Q3;
	do i = 1 to 1000;
	   X = rand('Uniform');
	   Y = rand('Uniform');
	   output;
	end;
RUN;

/* Question 3.3.A */
/* TODO */
/* 0.7 * 0.7 = 0.49 */
PROC CORR data=WEEK3_Q3 kendall spearman 
	plots=scatter(ellipse=none); 
	var X Y;
RUN;


/* QUESTION 3.4 */
DATA RCT;
	SET SASDATA.RCT;
/* 	RESP = RESP + (ranuni(1)-0.5); */
RUN;

/* Question 3.4.A */
PROC TRANSPOSE out=WIDE_RCT(drop = _NAME_ _LABEL_)
			   data=RCT prefix=RESP;
	by ID;
	id TIME;
	var RESP;
RUN;


/* Question 3.4.B */
PROC CORR data=WIDE_RCT pearson kendall spearman;
	var RESP1 RESP2;
RUN;


/* Question 3.4.D */
PROC CORR data=WIDE_RCT pearson fisher(biasadj=no);
	var RESP1 RESP2;
RUN;


/* Question 3.4.E */
%KendallTau(tau=.35941);
%SpearmanRho(rho=.50477);

PROC RANK data=WIDE_RCT out=ranked_a;
	var RESP1 RESP2;
	ranks rank_RESP1 rank_RESP2;
RUN;
PROC MEANS data=ranked_a n;
	var rank_RESP1 rank_RESP1;
RUN;
DATA marginals;
	set ranked_a;
	RESP1=rank_RESP1/716; /* or 236 if errors*/
	RESP2=rank_RESP1/716;
RUN;
PROC SGPLOT data=WIDE_RCT aspect=1;
	title "Original RESP1 and RESP2";
	scatter x=RESP1 y=RESP2 / markerattrs=(symbol=circlefilled color='red' size=8 );
RUN;
%SIM_FGM(nsim=1000, alpha=1.61735, seed=12345);
PROC SGPLOT data=Fgmc aspect=1;
	title "Simulated FGM copula";
	scatter x=X y=Y / markerattrs=(symbol=circlefilled color='black' size=8 );
RUN;

proc print data=marginals;run;
/* 716 */














































