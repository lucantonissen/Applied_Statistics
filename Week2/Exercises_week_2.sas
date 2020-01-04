LIBNAME SASDATA "/folders/myfolders/Applied_Statistics/Data";

DATA WEEK2;
	set SASDATA.IVF;
	where PER=4;
	drop IMP PER AGE;
RUN;


/* QUESTION 2.1 */

/* Question 2.1.A */
/* F-test (using T-test) */
PROC TTEST data=WEEK2;
	class FIS;
	var AGEM;
RUN;
/* Bartlett's & Levene's test */
PROC GLM data=WEEK2;
	class FIS;
	model AGEM = FIS;
	means FIS / hovtest=BARTLETT hovtest=levene;
RUN;

/* Question 2.1.B */
/* T-test from Q2.1.A */


/* QUESTION 2.2 */
PROC SQL;
	CREATE TABLE Q2_2 AS
		SELECT TRT, BW
		FROM WEEK2
		WHERE ID <= 100;
RUN;

/* Question 2.2.A */
ods output WilcoxonScores=Q2_2A_01 (keep= Class N SumOfScores);
PROC NPAR1WAY data=Q2_2 wilcoxon correct=NO;
	where trt < 2;
	class TRT;
	var BW;
	title 'TRT < 2';
	ods select WilcoxonTest; /* WilcoxonMC for monte carlo estimates */
RUN;
ods output WilcoxonScores=Q2_2A_02 (keep= Class N SumOfScores);
PROC NPAR1WAY data=Q2_2 wilcoxon correct=NO;
	where trt <> 1;
	class TRT;
	var BW;
	title 'TRT <> 1';
	ods select WilcoxonTest;
RUN;
ods output WilcoxonScores=Q2_2A_12 (keep= Class N SumOfScores);
PROC NPAR1WAY data=Q2_2 wilcoxon correct=NO;
	where trt > 0;
	class TRT;
	var BW;
	title 'TRT > 0';
	ods select WilcoxonTest;
RUN;

/* Question 2.2.B */
/* get p-value from wilcoxon-mann-whitney Q2_2A */


/* Question 2.2.C */
PROC IML;
	use Q2_2A_01;
		read all var{Class N SumOfScores};
	close Q2_2A_01;
	
	G = Num(Class); /*{0, 1};*/
	U = SumOfScores - N#(N+1)/2;
	P = U / prod(N);
	
	A = G||N||U||P;
	create MWU from A [colname={'Group' 'N' 'U' 'P'}];
		append from A;
	close MWU;
quit;
PROC IML;
	use Q2_2A_02;
		read all var{Class N SumOfScores};
	close Q2_2A_02;
	
	G = Num(Class); /*{0, 2};*/
	U = SumOfScores - N#(N+1)/2;
	P = U / prod(N);
	
	A = G||N||U||P;
	create MWU from A [colname={'Group' 'N' 'U' 'P'}];
		append from A;
	close MWU;
quit;


/* Question 2.2.D */
PROC NPAR1WAY data=Q2_2;
	where trt < 2;
	class TRT;
	var BW;
	exact KS /mc;
	title 'TRT < 2';
	ods select KSMC;
RUN;
PROC NPAR1WAY data=Q2_2;
	where trt <> 1;
	class TRT;
	var BW;
	exact KS /mc;
	title 'TRT <> 1';
	ods select KSMC;
RUN;
PROC NPAR1WAY data=Q2_2;
	where trt > 0;
	class TRT;
	var BW;
	exact KS /mc;
	title 'TRT > 0';
	ods select KSMC;
RUN;


/* Question 2.2.E */
DATA Q2_2E; 
	set Q2_2;
	BW_BOXCOX2 = (0.5)*(BW**(2) -1);
RUN;

PROC TTEST data=Q2_2E;
	where trt < 2;
	class TRT;
	var BW_BOXCOX2;
	ods select TTests;
RUN;
PROC TTEST data=Q2_2E;
	where trt <> 1;
	class TRT;
	var BW_BOXCOX2;
	ods select TTests;
RUN;
PROC TTEST data=Q2_2E;
	where trt > 0;
	class TRT;
	var BW_BOXCOX2;
	ods select TTests;
RUN;
/* Use pooled statistics for equal variances */


/* Question 2.2.F */
/* Perform Q2_2 A, D, E on the whole dataset */










































