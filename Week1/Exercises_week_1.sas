LIBNAME SASDATA "/folders/myfolders/Applied_Statistics/Data";

DATA WEEK1;
	set SASDATA.IVF;
	where PER=4;
	drop IMP PER AGE;
RUN;


/* Question 1.1 */
PROC FREQ data = week1;
	tables TRT;
run;


/* Question 1.2 */

/* Question 1.2.A */
PROC IML;
	use WEEK1;
	read all var{agem};
	close WEEk1;
	
	alpha=0.05;
	Ybar=mean(agem);
	s=var(agem);
	n=nrow(agem);
	qT=quantile('t', alpha/2, n-1);
	UCL = Ybar - qT * sqrt(s/n);
	LCL = Ybar + qT * sqrt(s/n);
	UPL = Ybar - qT * sqrt((n+1)*s/n); 
	LPL = Ybar + qT * sqrt((n+1)*s/n); 
	A=Ybar||s||LCL||UCL||LPL||UPL;
	
	create DATA from A[colname={'mean' 'variance' 'LCL' 'UCL' 'LPL' 'UPL'}];
		append from A;
	close DATA;
QUIT;

PROC MEANS data=WEEK1 mean var;
	var AGEM;
RUN;


/* Question 1.2.B */
/* Mothers 40y and older */
PROC SQL;
	SELECT COUNT(AGEM)
	FROM WEEK1
	WHERE AGEM >= 40;
QUIT;


/* Question 1.2.C */
PROC IML;
	use WEEK1;
	read all var{agem};
	close WEEK1;
	
	alpha=0.05;
	s=var(agem);
	n=nrow(agem);
	qCL=quantile('chisquare',
	alpha/2,n-1);
	qCU=quantile('chisquare',
	1-alpha/2,n-1);
	UCL=(n-1)*s/qCL;
	LCL=(n-1)*s/qCU;
	sd=sqrt(s);
	UCLsd=sqrt((n-1)*s/qCL);
	LCLsd=sqrt((n-1)*s/qCU);
	
	A=(s||LCL||UCL)//(sd||LCLsd||UCLsd);
	create SD from A[colname={'statistic' 'LCL' 'UCL'}];
		append from A;
	close SD;
QUIT;

ods select BasicIntervals;
PROC UNIVARIATE data=WEEK1 cibasic(alpha=0.05);
   var AGEM;
RUN;

PROC UNIVARIATE data=WEEK1 cibasic;
	var AGEM;
	histogram AGEM/ normal;
RUN;


/* Question 1.3 */

/* Question 1.3.A */
/* Quartile confidence interval (CI)*/
/* Interquartile range */
PROC IML;
	use WEEK1;
	read all var{BW};
	close WEEK1;
	
	alpha=0.05;
	p=0.75;		/* 0.25 first quartile */
	s=p*(1-p);
	n=nrow(BW);
	z=quantile('Normal', 1-alpha/2);
	pU=p + z * sqrt(s/n);
	pL=p - z * sqrt(s/n);
	nU=min(floor(n*pU)+1, n);
	nL=max(1, floor(n*pL));
	call sort(BW);
	call qntl(pct, BW, p);
	LCL=BW[nL];		/* lower control limit */
	UCL=BW[nU];		/* upper control limit */
	
	A=(pct||LCL||UCL||nL||nU);
	create PCTL from A[colname={'pctl' 'LCL' 'UCL' 'LR' 'UR'}];
		append from A;
	close PCTL;
QUIT;

PROC UNIVARIATE data=WEEK1 cipctldf; /* cipctlnormal; */
	var BW;
RUN;


/* Question 1.3.B */
PROC UNIVARIATE data=WEEK1 cipctldf; /* cipctlnormal; */
	var BW;
RUN;	/* get median and interquartile range */

PROC SQL;
	SELECT COUNT(BW) / (SELECT count(*) FROM WEEK1)
	FROM WEEK1 
	WHERE BW BETWEEN 3365 - 870 AND 3365 + 870; 
QUIT;	/* meadian +- interquartile range */


/* Question 1.3.C */
DATA WEEK1BOXCOX;
	SET WEEK1;
	BWMINUS2= (-1/2)*(BW**-2 -1);
	BWMINUS1= (-1)*(BW**-1 -1);
	BWMINUS12= (-2)*(BW**-(0.5)-1);
	BW0= log(BW);
	BW13= (3)*(BW**(1/3) -1);
	BW12= (2)*(BW**(1/2) -1);
	BW2= (0.5)*(BW**(2) -1);
RUN;
PROC UNIVARIATE data=WEEK1BOXCOX;
	histogram BWMINUS2 /normal;
	histogram BWMINUS1 /normal;
	histogram BWMINUS12 /normal;
	histogram BW0 /normal;
	histogram BW13 /normal;
	histogram BW12 /normal;
	histogram BW2 /normal;
	histogram BW /normal;
RUN;


/* Question 1.3.D */
PROC IML;	
	use WEEK1BOXCOX;
	read all var{BW2};
	close WEEK1BOXCOX;
	
	alpha=0.05;
	Ybar=mean(BW2);
	s=var(BW2);
	n=nrow(BW2);
	qT=quantile('t', alpha/2, n-1);
	UCL = Ybar - qT * sqrt(s/n);
	LCL = Ybar + qT * sqrt(s/n);
	UPL = Ybar - qT * sqrt((n+1)*s/n); 
	LPL = Ybar + qT * sqrt((n+1)*s/n); 
	
	Ybar = sqrt(2*Ybar + 1);
	UCL = sqrt(2*UCL + 1);
	LCL = sqrt(2*LCL + 1);
	UPL = sqrt(2*UPL + 1);
	LPL = sqrt(2*LPL + 1);
	A=Ybar||LCL||UCL||LPL||UPL;
	
	create DATA from A[colname={'mean' 'LCL' 'UCL' 'LPL' 'UPL'}];
		append from A;
	close DATA;
QUIT;


/* Question 1.3.E */
PROC SQL;
	SELECT BW, SEX
	FROM WEEK1 
	WHERE BW = (SELECT MAX(BW) FROM WEEK1);
QUIT;


/* Question 1.3.F */
/* percentages */
PROC SQL;
	SELECT AVG(SEX) as BOY_PERC, 1 - AVG(SEX) as GIRL_PERC
	FROM WEEK1;
QUIT;

/* create boy and girl datasets */
PROC SQL;
	CREATE TABLE WEEK1BOY AS
		SELECT *
		FROM WEEK1
		WHERE SEX = 1;
QUIT;
PROC SQL;
	CREATE TABLE WEEK1GIRL AS
		SELECT *
		FROM WEEK1
		WHERE SEX = 0;
QUIT;
/* Data WEEK1GIRL; */
/* 	set WEEK1; */
/* 	where sex = 0; */
/* RUN; */

/* parameters */
PROC UNIVARIATE data=WEEK1BOY;
	var BW;
RUN;
PROC UNIVARIATE data=WEEK1GIRL;
	var BW;
RUN;


/* Question 1.3.G */
DATA WEEK1GIRLBOXCOX;
	SET WEEK1GIRL;
	BW2= (0.5)*(BW**(2) -1);
RUN;
PROC IML;	
	use WEEK1GIRLBOXCOX;
	read all var{BW2};
	close WEEK1GIRLBOXCOX;
/* 	BW=(0.5)*(BW**(2) -1); */
	
	alpha=0.05;
	Ybar=mean(BW2);
	s=var(BW2);
	n=nrow(BW2);
	qT=quantile('t', alpha/2, n-1);
	UCL = Ybar - qT * sqrt(s/n);
	LCL = Ybar + qT * sqrt(s/n);
	UPL = Ybar - qT * sqrt((n+1)*s/n); 
	LPL = Ybar + qT * sqrt((n+1)*s/n); 
	
	Ybar = sqrt(2*Ybar + 1);
	UCL = sqrt(2*UCL + 1);
	LCL = sqrt(2*LCL + 1);
	UPL = sqrt(2*UPL + 1);
	LPL = sqrt(2*LPL + 1);
	A=Ybar||LCL||UCL||LPL||UPL;
	
	create DATA from A[colname={'mean' 'LCL' 'UCL' 'LPL' 'UPL'}];
		append from A;
	close DATA;
QUIT;

DATA WEEK1BOYBOXCOX;
	SET WEEK1BOY;
	BW2= (0.5)*(BW**(2) -1);
RUN;
PROC IML;	
	use WEEK1BOYBOXCOX;
	read all var{BW2};
	close WEEK1BOYBOXCOX;
/* 	BW=(0.5)*(BW**(2) -1); */
	
	alpha=0.05;
	Ybar=mean(BW2);
	s=var(BW2);
	n=nrow(BW2);
	qT=quantile('t', alpha/2, n-1);
	UCL = Ybar - qT * sqrt(s/n);
	LCL = Ybar + qT * sqrt(s/n);
	UPL = Ybar - qT * sqrt((n+1)*s/n); 
	LPL = Ybar + qT * sqrt((n+1)*s/n); 
	
	Ybar = sqrt(2*Ybar + 1);
	UCL = sqrt(2*UCL + 1);
	LCL = sqrt(2*LCL + 1);
	UPL = sqrt(2*UPL + 1);
	LPL = sqrt(2*LPL + 1);
	A=Ybar||LCL||UCL||LPL||UPL;
	
	create DATA from A[colname={'mean' 'LCL' 'UCL' 'LPL' 'UPL'}];
		append from A;
	close DATA;
QUIT;


/* Question 1.4 */
DATA Q1_4; 
   INPUT Value; 
   DATALINES; 
	25.0 
	27.4 
	17.1 
	22.1 
	20.8 
	21.3 
	22.5 
	29.2 
	27.9 
	25.7 
	24.7 
	18.8
	;
RUN; 


/* Question 1.4.A */
PROC UNIVARIATE data=Q1_4;
	var Value;
RUN;


/* Question 1.4.B */
PROC UNIVARIATE data=Q1_4 cibasic;
	var Value;
RUN;

PROC IML;
	use Q1_4;
	read all var{Value};
	close Q1_4;
	
	alpha=0.05;
	Ybar=mean(Value);
	s=var(Value);
	n=nrow(Value);
	
	qT=quantile('t', alpha/2, n-1);
	meanUCL = Ybar - qT * sqrt(s/n);
	meanLCL = Ybar + qT * sqrt(s/n); 
	
	qCL=quantile('chisquare', alpha/2,n-1);
	qCU=quantile('chisquare', 1-alpha/2,n-1);
	sUCL=(n-1)*s/qCL;
	sLCL=(n-1)*s/qCU;
	sd=sqrt(s);
	sdUCL=sqrt((n-1)*s/qCL);
	sdLCL=sqrt((n-1)*s/qCU);
	
	A=Ybar||sd||s||meanLCL||meanUCL||sdLCL||sdUCL||sLCL||sUCL;
	
	create DATA from A[colname={
		'mean' 'sd' 'variance' 
		'mean_LCL' 'mean_UCL' 
		'sd_LCL' 'sd_UCL'
		'variance_LCL' 'variance_UCL'}];
		append from A;
	close DATA;
QUIT;


/* Question 1.4.C */
PROC IML;
	use Q1_4;
	read all var{Value};
	close Q1_4;
	
	alpha=0.01;
	Ybar=mean(Value);
	s=var(Value);
	n=nrow(Value);
	qT=quantile('t', alpha/2, n-1);
	UPL = Ybar - qT * sqrt((n+1)*s/n); 
	LPL = Ybar + qT * sqrt((n+1)*s/n); 
	A=LPL||UPL;
	
	create DATA from A[colname={'LPL' 'UPL'}];
		append from A;
	close DATA;
QUIT;


















































































