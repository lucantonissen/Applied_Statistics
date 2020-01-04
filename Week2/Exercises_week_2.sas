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


/* QUESTION 2.3 */

/* Question 2.3.A */
%MACRO mann_whitney_u(dataset, class, var);
	ods select none;
	PROC NPAR1WAY data=&dataset wilcoxon correct=no;
		var &var;
		class &class;
		ods output WilcoxonScores=OUT_SCORES(rename=(SumOfScores=S));
		ods output WilcoxonTest=OUT_TEST;
		ods output KruskalWallisTest=OUT_KRUS;
	RUN;
	ods select all;
	
	/* proc print data=OUT_TEST;run; */
	
	/* Chi square p-value */
	PROC SQL;
		CREATE TABLE P_KRUS AS
			SELECT Prob2 FROM OUT_KRUS;
	RUN;
	
	/* Wilcoxon p-value */
	PROC SQL;
		CREATE TABLE P_TABLE AS
			SELECT Prob2 FROM OUT_TEST;
	RUN;
	
	DATA OUT_SCORES;
		set OUT_SCORES; 
		CLASS_ID = _N_ - 1;
	RUN;
	
	PROC TRANSPOSE data=OUT_SCORES
			out=OUT_N(drop=_NAME_) prefix=N;
		id CLASS_ID;
		var N;
	RUN;
	
	PROC TRANSPOSE data=OUT_SCORES
			out=OUT_S(drop=_NAME_ _LABEL_) PREFIX=S;
		id CLASS_ID;
		var S;
	RUN;
	
	DATA RESULT;
		merge OUT_N OUT_S P_TABLE P_KRUS;
		P_VALUE = Prob2;
		P_KRUS = Prob;
		U0 = S0 - N0 * (N0+1)/2;
		U1 = S1 - N1 * (N1+1)/2;
		P0 = U0 / (N0*N1);
		P1 = U1 / (N0*N1);
	RUN;
	
	title "Mann Whitney U test";
	PROC PRINT data=OUT_SCORES label noobs;
		var CLASS_ID CLASS; 
		label CLASS_ID="class"
		CLASS="group identifier";
	RUN;

	PROC PRINT data=RESULT label;
		var P_VALUE P_KRUS U0 U1 P0 P1;
		label P_VALUE = "p-value Wilcoxon Test"
		P_KRUS= "p-value Kruskal-Wallis Test"
		U0="statistic (U0)" 
		U1="statistic (U1)" 
		P0="P(class0 > class1)" 
		P1="P(class0 <= class1)";
	RUN;
	title;
%MEND;


/* Question 2.3.B */
DATA Q2_3B_01;
	set Q2_2; 
	if TRT < 2;
RUN;
%mann_whitney_u(Q2_3B_01, TRT, BW);
DATA Q2_3B_02;
	set Q2_2; 
	if TRT <> 1;
RUN;
%mann_whitney_u(Q2_3B_02, TRT, BW);
DATA Q2_3B_12;
	set Q2_2; 
	if TRT > 0;
RUN;
%mann_whitney_u(Q2_3B_12, TRT, BW);


/* QUESTION 2.4 */

/* Question 2.4.A & 2.4.C*/
DATA Q2_4;
	set WEEK2;
	heavy = BW > 4000;
	late = GA > 41; 
	keep GA BW heavy late;
RUN;


/* Question 2.4.B */
PROC TTEST data=Q2_4;
	class heavy;
	var GA;
RUN;


/* Question 2.4.D */
PROC TTEST data=Q2_4;
	class late;
	var BW;
RUN;


/* Question 2.4.E */
/* Most likely reliable (CLT) */


/* QUESTION 2.5 */
DATA Q5;
	input treatment$ high$ count;
	datalines;
	0 0 77
	0 1 23
	1 0 81
	1 1 19
	;
RUN;
PROC FREQ data=Q5;
	table treatment * high /chisq;
	weight count ;
	exact chisq;
RUN;






































