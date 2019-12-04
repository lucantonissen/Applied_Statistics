LIBNAME SASDATA "/folders/myfolders/Applied_Statistics/Data";

DATA WEEK1;
	set SASDATA.IVF;
	where PER=4;
RUN;


/* Question 1.1 */
PROC FREQ data = week1;
	tables TRT;
run;


/* Question 1.2.A */
PROC MEANS data=WEEK1 mean var;
	var AGEM;
RUN;

PROC IML;
use WEEK1;
read all var{agem};
close Week1;

alpha=0.05;
Ybar=mean(agem);
s=var(agem);
n=nrow(agem);
qT=quantile('t', alpha/2, n-1);
UCL = Ybar - qT * sqrt(s/n);
LCL = Ybar + qT * sqrt(s/n);
UPL = Ybar - qT * sqrt((n+1)*s/n); 
LPL = Ybar + qT * sqrt((n+1)*s/n); 
A=Ybar||LCL||UCL||LPL||UPL;

create DATA from A[colname={'mean' 'LCL' 'UCL' 'LPL' 'UPL'}];
append from A;
close DATA;
quit;