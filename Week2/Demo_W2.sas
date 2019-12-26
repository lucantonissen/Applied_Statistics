LIBNAME SASDATA "/folders/myfolders/Applied_Statistics/Data";

data IVF_DATASET;
	set SASDATA.IVF;
	where PER=4;
	drop IMP PER AGE;
run;

data Q16;
	SET IVF_DATASET;
	LGA=LOG(44-GA);
run;

ods select histogram;
proc univariate data=Q16;
histogram GA /normal;
histogram LGA/normal;
run;
	
*Q1a;

PROC IML;
use Q16;
read all var{LGA};
close Q16;

alpha=0.05;
Ybar=mean(LGA);
s=var(LGA);
n=nrow(LGA);
qT=quantile('t',alpha/2,n-1);
UPL=Ybar-qT*sqrt((n+1)*s/n);
LPL=Ybar+qT*sqrt((n+1)*s/n);
A=Ybar||LPL||UPL;

create PIa from A[colname={'mean' 'LPL' 'UPL'}]; 
append from A;       
close PIa;
quit;

data PIa;
set PIa;
LPLGA=44-exp(UPL);
UPLGA=44-exp(LPL);
run;

*b;
proc univariate data=Q16 cibasic;
	var LGA;
run;
*1.44064	1.55192;

*c;
proc univariate data=Q16 cipctldf;
var GA;
run;

