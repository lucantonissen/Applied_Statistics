LIBNAME SASDATA "/folders/myfolders/Applied_Statistics/Data";

DATA WEEK1;
	set SASDATA.IVF;
	where PER=4;
RUN;


/* View data */
PROC PRINT data=WEEK1 NOOBS;


/* Check normal distribution */
/* ods select HISTOGRAM, QQPLOT, PPPLOT; */
PROC UNIVARIATE data=WEEK1;
   VAR AGEM;
   HISTOGRAM AGEM/NORMAL; /* change distributions */
   QQPLOT    AGEM/NORMAL;
   PPPLOT 	 AGEM/NORMAL;
RUN;


/* Frequency */
PROC FREQ data=WEEK1;
	tables AGEM;
Run;


/* Prediction Interval Regression*/
proc reg data=WEEK1 ;
	model AGEM= / cli alpha=0.05;
run;


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