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
	CREATE TABLE WEEK2_Q2 AS
		SELECT TRT, BW
		FROM WEEK2
		WHERE ID <= 100;
RUN;

/* Question 2.2.A */
PROC NPAR1WAY data=WEEK2_Q2 wilcoxon correct=NO;
	where trt < 2;
	class TRT;
	var BW;
	title 'TRT < 2';
	ods select WilcoxonTest;
RUN;
PROC NPAR1WAY data=WEEK2_Q2 wilcoxon correct=NO;
	where trt <> 1;
	class TRT;
	var BW;
	title 'TRT <> 1';
	ods select WilcoxonTest;
RUN;
PROC NPAR1WAY data=WEEK2_Q2 wilcoxon correct=NO;
	where trt > 0;
	class TRT;
	var BW;
	title 'TRT > 1';
	ods select WilcoxonTest;
RUN;

/* Question 2.2.B */

	
	























































