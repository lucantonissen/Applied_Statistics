LIBNAME SASDATA "/folders/myfolders/Applied_Statistics/Data";

DATA WEEK1;
	set SASDATA.IVF;
	where PER = 4;
	drop IMP PER AGE;
RUN;

PROC PRINT data=WEEK1;
RUN;

/* QUESTION 1.1 */
PROC FREQ data = week1;
	tables TRT;
run;


/* QUESTION 1.2 */

/* Question 1.2.A */
PROC IML;
	use WEEK1;
		read all var{agem};
	close WEEk1;
	
	alpha = 0.05;
	Ybar = mean(agem);
	s = var(agem);
	n = nrow(agem);
	qT = quantile('t', alpha/2, n-1);
	
	LCL = Ybar + qT * sqrt(s/n);
	UCL = Ybar - qT * sqrt(s/n);
	LPL = Ybar + qT * sqrt((n+1)*s/n);
	UPL = Ybar - qT * sqrt((n+1)*s/n); 
	
	A = Ybar||s||LCL||UCL||LPL||UPL;
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
	
	alpha = 0.05;
	s = var(agem);
	n = nrow(agem);
	qCU = quantile('chisquare', 1-alpha/2,n-1);
	qCL = quantile('chisquare', alpha/2,n-1);
	
	LCL = (n-1)*s/qCU;
	UCL = (n-1)*s/qCL;
	
	sd = sqrt(s);
	LCLsd = sqrt((n-1)*s/qCU);
	UCLsd = sqrt((n-1)*s/qCL);
	
	A = (s||LCL||UCL)//(sd||LCLsd||UCLsd);
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


/* QUESTION 1.3 */

/* Question 1.3.A */
/* Quartile confidence interval (CI)*/
/* Interquartile range */
PROC IML;
	use WEEK1;
		read all var{BW};
	close WEEK1;
	
	alpha = 0.05;
	p = 0.75;		/* 0.25 first quartile */
	s = p*(1-p);
	n = nrow(BW);
	z = quantile('Normal', 1-alpha/2);

	pL = p - z * sqrt(s/n);
	pU = p + z * sqrt(s/n);
	nL = max(1, floor(n*pL));
	nU = min(floor(n*pU)+1, n);

	call sort(BW);
	call qntl(pct, BW, p);
	UCL = BW[nU];		/* upper control limit */
	LCL = BW[nL];		/* lower control limit */
	
	A = (pct||LCL||UCL||nL||nU);
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
	BWMINUS2	= (-1/2)*(BW**-2 -1);
	BWMINUS1	= (-1)*(BW**-1 -1);
	BWMINUS12	= (-2)*(BW**-(0.5)-1);
	BW0			= log(BW);
	BW13		= (3)*(BW**(1/3) -1);
	BW12		= (2)*(BW**(1/2) -1);
	BW2			= (0.5)*(BW**(2) -1);
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
	
	alpha = 0.05;
	Ybar = mean(BW2);
	s = var(BW2);
	n = nrow(BW2);
	qT = quantile('t', alpha/2, n-1);
	
	LCL = Ybar + qT * sqrt(s/n);
	UCL = Ybar - qT * sqrt(s/n);
	LPL = Ybar + qT * sqrt((n+1)*s/n); 
	UPL = Ybar - qT * sqrt((n+1)*s/n); 
	
	Ybar = sqrt(2*Ybar + 1);
	LCL = sqrt(2*LCL + 1);
	UCL = sqrt(2*UCL + 1);
	LPL = sqrt(2*LPL + 1);
	UPL = sqrt(2*UPL + 1);
	
	A = Ybar||LCL||UCL||LPL||UPL;
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
	BW2 = (0.5)*(BW**(2) -1);
RUN;
PROC IML;	
	use WEEK1GIRLBOXCOX;
		read all var{BW2};
	close WEEK1GIRLBOXCOX;
	
	alpha = 0.05;
	Ybar = mean(BW2);
	s = var(BW2);
	n = nrow(BW2);
	qT = quantile('t', alpha/2, n-1);

	LCL = Ybar + qT * sqrt(s/n);
	UCL = Ybar - qT * sqrt(s/n);
	LPL = Ybar + qT * sqrt((n+1)*s/n); 
	UPL = Ybar - qT * sqrt((n+1)*s/n); 
	
	Ybar = sqrt(2*Ybar + 1);
	LCL = sqrt(2*LCL + 1);
	UCL = sqrt(2*UCL + 1);
	LPL = sqrt(2*LPL + 1);
	UPL = sqrt(2*UPL + 1);
	
	A = Ybar||LCL||UCL||LPL||UPL;
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
	
	alpha = 0.05;
	Ybar = mean(BW2);
	s = var(BW2);
	n = nrow(BW2);
	qT = quantile('t', alpha/2, n-1);

	LCL = Ybar + qT * sqrt(s/n);
	UCL = Ybar - qT * sqrt(s/n);
	LPL = Ybar + qT * sqrt((n+1)*s/n); 
	UPL = Ybar - qT * sqrt((n+1)*s/n); 
	
	Ybar = sqrt(2*Ybar + 1);
	LCL = sqrt(2*LCL + 1);
	UCL = sqrt(2*UCL + 1);
	LPL = sqrt(2*LPL + 1);
	UPL = sqrt(2*UPL + 1);
	
	A = Ybar||LCL||UCL||LPL||UPL;
	create DATA from A[colname={'mean' 'LCL' 'UCL' 'LPL' 'UPL'}];
		append from A;
	close DATA;
QUIT;


/* QUESTION 1.4 */
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
	
	alpha = 0.05;
	Ybar = mean(Value);
	s = var(Value);
	n = nrow(Value);
	
	qT = quantile('t', alpha/2, n-1);
	meanLCL = Ybar + qT * sqrt(s/n); 
	meanUCL = Ybar - qT * sqrt(s/n);
	
	qCL = quantile('chisquare', alpha/2,n-1);
	qCU = quantile('chisquare', 1-alpha/2,n-1);
	sLCL = (n-1)*s/qCU;
	sUCL = (n-1)*s/qCL;
	sd = sqrt(s);
	sdLCL = sqrt((n-1)*s/qCU);
	sdUCL = sqrt((n-1)*s/qCL);
	
	A = Ybar||sd||s||meanLCL||meanUCL||sdLCL||sdUCL||sLCL||sUCL;
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
	
	alpha = 0.01;
	Ybar = mean(Value);
	s = var(Value);
	n = nrow(Value);
	qT = quantile('t', alpha/2, n-1);
	
	LPL = Ybar + qT * sqrt((n+1)*s/n); 
	UPL = Ybar - qT * sqrt((n+1)*s/n); 
	
	A = LPL||UPL;
	create DATA from A[colname={'LPL' 'UPL'}];
		append from A;
	close DATA;
QUIT;


/* QUESTION 1.5 */

/* Question 1.5.A */
PROC IML;
	logmeanLCL = -0.137;
	logmeanUCL = 1.128;
	
	meanLCL = exp(1)**logmeanLCL;
	meanUCL = exp(1)**logmeanUCL;
	
	A = meanLCL||meanUCL;
	create DATA from A[colname={'LPL' 'UPL'}];
		append from A;
	close DATA;
QUIT;


/* Question 1.5.B */
/* Can't reject */


/* QUESTION 1.6 */
DATA Q1_6;
	SET WEEK1;
	GAlog44 = log(44 - GA);
RUN;


/* Question 1.6.A */
PROC IML;	
	use Q1_6;
		read all var{GAlog44};
	close Q1_6;
	
	alpha=0.05;
	Ybar=mean(GAlog44);
	s=var(GAlog44);
	n=nrow(GAlog44);
	qT=quantile('t', alpha/2, n-1);
	
	LPL = Ybar - qT * sqrt((n+1)*s/n); 
	UPL = Ybar + qT * sqrt((n+1)*s/n);
	
	LPL = 44 - exp(1)**(LPL);
	UPL = 44 - exp(1)**(UPL);
	
	A = LPL||UPL;
	create DATA from A[colname={'LPL' 'UPL'}];
		append from A;
	close DATA;
QUIT;


/* Question 1.6.B */
PROC IML;	
	use Q1_6;
		read all var{GAlog44};
	close Q1_6;
	
	alpha=0.05;
	Ybar=mean(GAlog44);
	s=var(GAlog44);
	n=nrow(GAlog44);
	qT=quantile('t', alpha/2, n-1);
	
	LCL = Ybar + qT * sqrt(s/n);
	UCL = Ybar - qT * sqrt(s/n);
	
	A = LCL||UCL;
	create DATA from A[colname={'LCL' 'UCL'}];
		append from A;
	close DATA;
QUIT;


/* Question 1.6.C */
ods select Quantiles;
PROC UNIVARIATE data=Q1_6 CIPCTLDF CIPCTLNORMAL;
	var GA;
RUN;


/* Question 1.6.D */
DATA Q1_6D;
	set WEEK1;
	if GA <= 38 then PRETERM = 1;
	else PRETERM = 0;
RUN;


/* Question 1.6.E */
PROC IML;	
	use Q1_6D where(PRETERM=1);
		read all var{FIS};
	close Q1_6D;
	
	avg = mean(FIS);
	
	A = avg;
	create DATA from A[colname={'percentage'}];
		append from A;
	close DATA;
RUN;
	
PROC SQL;
	SELECT avg(FIS)
	FROM Q1_6D
	WHERE PRETERM = 1;
QUIT;
PROC SQL;
	SELECT avg(FIS)
	FROM Q1_6D
	WHERE PRETERM = 0;
QUIT;


/* Question 1.6.F */
PROC FREQ data=Q1_6D; /* use this because small n */
	where PRETERM = 1;
	tables FIS /binomial(wald wilson exact level=2) alpha=0.1;
RUN; /* wald */
PROC FREQ data=Q1_6D;
	where PRETERM = 0;
	tables FIS /binomial(wald wilson exact level=2) alpha=0.1;
RUN;


/* Question 1.6.G */
/* Same as 1.6.F but use Clopper-Pearson (exact) */


/* QUESTION 1.7 */
%MACRO samples(dataset=,ns=,n=);
	PROC SURVEYSELECT data=&dataset NOPRINT
		method=urs n=&n out=FINAL;
	RUN;
	DATA FINAL;
		set FINAL;
		sampleno=1;
	RUN;
		%do sn = 2 %to &ns;
	PROC SURVEYSELECT data=&dataset NOPRINT
		method=urs n=&n out=SAMPLEI;
	RUN;
	DATA SAMPLEI;
		set SAMPLEI;
		sampleno= &sn;
	RUN;
	DATA FINAL;
		set Final SAMPLEI;
	RUN;
	%end;
	PROC DATASETS library=work NOPRINT;
		delete SAMPLEI;
	RUN;
%mend;

%samples(dataset=WEEK1, ns=1000, n=10);

/* Question 1.7.A */
PROC SQL;
	CREATE TABLE FINAL_A AS
		SELECT sampleno,
			mean(AGEM) as AGEM_mean,
			var(AGEM) as AGEM_variance
		FROM FINAL
		GROUP BY sampleno;
QUIT;


/* Question 1.7.B */
ods select histogram;
PROC UNIVARIATE data=FINAL_A cibasic;
	histogram AGEM_mean/ normal;
	histogram AGEM_variance/ normal;
RUN;

/* ????? */
PROC SQL;
	CREATE TABLE FINAL_B AS
		SELECT AGEM_mean, AGEM_variance,
			(AGEM_mean - 32.668355369) / sqrt(12.137436726) as AGEM_mean_standardized,
			AGEM_variance / 12.137436726 as AGEM_scaled_sample_variance
		FROM FINAL_A;
QUIT;


/* Question 1.7.C */
ods select histogram;
PROC UNIVARIATE data=FINAL_B cibasic;
	histogram AGEM_variance/ normal;
	histogram AGEM_scaled_sample_variance/ normal;
RUN;


/* QUESTION 1.8 */

/* Question 1.8.A */
PROC TTEST data=WEEK1 h0= 3200 sides=2 alpha=0.05;
	var BW;
RUN;
/* test statistic: 2.20 */
/* P value: 0.0285 */
/* ????? */
/* ¿reject because t > alpha? */


/* Question 1.8.B */
PROC TTEST data=WEEK1 h0=3200 sides=U alpha=0.05;
	var BW;
RUN;
/* P value: 0.0143 */


/* Question 1.8.C */
%samples(dataset=WEEK1, ns=100 , n=253);
PROC MEANS data=FINAL mean NOPRINT;
	var BW;
	by sampleno;
	output out=MEANSBW mean=BW_MEAN;
RUN;
PROC UNIVARIATE data=MEANSBW;
	hist BW_MEAN / normal;
RUN;
%mend;
/* Yes, it seems like 'n' is large enough to to assume that  */
/* the sample average 'y¯' follows a normal distribution */


/* Question 1.8.D */
/* The H0 hypothesis will be rejected faster */


/* Question 1.8.E */
PROC MEANS data=WEEK1 mean;
	var BW;
RUN;

PROC TTEST data=WEEK1 h0=3200 sides=U alpha=0.05;
	var BW;
RUN;
PROC TTEST data=WEEK1 h0=2800 sides=L alpha=0.05;
	var BW;
RUN;
/* 2 * min p val = 2 * 0.0143 = 0.0285 */


/* Question 1.8.F */
DATA WEEK1_BW;
	set WEEK1;
	keep BW;
RUN;
%samples(dataset=WEEK1_BW, ns=1000, n=50);

ods graphics off;
ods exclude all;
PROC TTEST data=FINAL h0=3200 sides=U alpha=0.05;
	var BW;
	by SAMPLENO;
	ods output ttests=POWER(keep=PROBT);
RUN;
ods exclude none;
ods graphics on;

PROC PRINT data=POWER;
RUN;


/* Question 1.8.G */
PROC SQL;
	SELECT count(*) / (SELECT count(*) FROM POWER)
	FROM POWER
	WHERE Probt <= 0.05;
QUIT;


/* Question 1.8.H */
/* n = 100 gives 0.327 */
/* n = 50  gives 0.225 */