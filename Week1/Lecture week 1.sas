
LIBNAME SASDATA "/folders/myfolders/Applied_Statistics/Data";

/*Question 1 */ 


DATA WEEK1;
 SET SASDATA.IVF;
 WHERE PER=4;
 AGEMB= (AGEM<30);
 KEEP AGEM AGEMB;
RUN;


*CODE LECTURES
*W1-1 example;
proc ttest data=WEEK1 h0= 32 sides=2 alpha=0.05;
	var AGEM;
run;


*W1-2 CI for mean;
proc univariate data=WEEK1 cibasic;
	var AGEM;
	histogram AGEM/ normal;
run;

PROC IML;
use WEEK1;
read all var{agem};
close WEEK1;

alpha=0.05;

Ybar=mean(agem);
s=var(agem);
n=nrow(agem);
qT=quantile('t',alpha/2,n-1);
UCL=Ybar-qT*sqrt(s/n);
LCL=Ybar+qT*sqrt(s/n);

A=Ybar||LCL||UCL;

print(Ybar||LCL||UCL);

create DATA from A[colname={'mean' 'LCL' 'UCL'}]; 
append from A;       
close DATA;
quit;

*W1-2 CLT simulation - IVF data;
%macro samples(dataset=,ns=,n=);
                proc surveyselect data=&dataset NOPRINT
   			method=urs n=&n out=Final;
		run;
		
		Data Final;
		set final;
		sampleno=1;
		run;

               
               %do sn = 2 %to &ns;
                        proc surveyselect data=&dataset NOPRINT
   			method=urs n=&n out=SampleI;
			run;
			
			Data SampleI;
			Set SampleI;
			sampleno= &sn;
			run;

			Data Final;
			Set Final SampleI;
			run;
     		%end;
     		
     		proc datasets library=work NOPRINT;
     		delete SampleI;
     		run;
%MEND;

%samples(dataset=WEEK1,ns=1000,n=10);

proc means data=Final mean NOPRINT;
               var agem;
               by sampleno;
               output out=meansA(drop=_type_ _freq_) mean=agem_mean var=agem_var N=n;
run;

Data meansA NOPRINT;
set meansA;
norm_agem_var = (n-1)*agem_var/(3.46854**2);
*n=100;
agem_sta = (agem_mean - 32.37194)/sqrt(agem_var/n);
run;

proc univariate data=meansA;
               var agem_mean;
               hist agem_mean /normal(mu=est);
run;

*W1-2 CI proportion;
proc freq data=WEEK1;
	tables AGEMB /binomial(wald wilson exact level=2) alpha=0.05;
run; 


PROC IML;
use WEEK1;
read all var{agemb};
close WEEK1;

alpha=0.05;
Ybar=mean(agemb);
s=Ybar*(1-Ybar);
n=nrow(agemb);
z=quantile('Normal',1-alpha/2);
UCL=Ybar+z*sqrt(s/n);
LCL=Ybar-z*sqrt(s/n);
A=Ybar||LCL||UCL;

print(Ybar||LCL||UCL);

create PROP from A[colname={'p' 'LCL' 'UCL'}]; 
append from A;       
close PROP;
quit;

*W1-2 CI for variance;
proc univariate data=WEEK1 cibasic;
	var AGEM;
	histogram AGEM/ normal;
run;

PROC IML;
use WEEK1;
read all var{agem};
close WEEK1;

alpha=0.05;

s=var(agem);
n=nrow(agem);
qCL=quantile('chisquare',alpha/2,n-1);
qCU=quantile('chisquare',1-alpha/2,n-1);
UCL=(n-1)*s/qCL;
LCL=(n-1)*s/qCU;

sd=sqrt(s);
UCLsd=sqrt((n-1)*s/qCL);
LCLsd=sqrt((n-1)*s/qCU);

A=(s||LCL||UCL)//(sd||LCLsd||UCLsd);

print(sd||LCLsd||UCLsd);

create SD from A[colname={'statistic' 'LCL' 'UCL'}]; 
append from A;       
close SD;
quit;



*W1-2 CI percentile;
proc univariate data=WEEK1 cipctldf;
	var AGEM;
run;

PROC IML;
use WEEK1;
read all var{agem};
close WEEK1;

alpha=0.05;
p=0.1;
s=p*(1-p);
n=nrow(agem);
z=quantile('Normal',1-alpha/2);
pU=p+z*sqrt(s/n);
pL=p-z*sqrt(s/n);
nU=min(floor(n*pU)+1,n);
nL=max(1,floor(n*pL));

*http://support.sas.com/documentation/cdl/en/imlug/67502/HTML/default/viewer.htm#imlug_langref_sect431.htm;
call sort(agem);
call qntl(pct, agem, p);
LCL=agem[nL];
UCL=agem[nU];
A=(pct||LCL||UCL||nL||nU);


create PCTL from A[colname={'pctl' 'LCL' 'UCL' 'LR' 'UR'}]; 
append from A;       
close PCTL;
quit;




*W1-2 PI example;
proc reg data=WEEK1;
     model AGEM= / cli alpha=0.05;
run;

proc means data=WEEK1 mean std n;
	var AGEM;
	output out=agem_sumstat;
run;

proc transpose data=agem_sumstat out=agem_PI (DROP= _TYPE_ _FREQ_ _NAME_ _LABEL_);
	by _type_ _freq_;
	id _stat_;
run;

data agem_PI;
	set agem_PI;
	T    = QUANTILE("T", 1 - 0.05/2, N-1);
	LPL = MEAN - T * std*SQRT((N+1)/ N);
	UPL = MEAN + T * std*SQRT((N+1)/ N);
run;

PROC PRINT DATA=agem_PI;
	var LPL UPL;
run;

PROC IML;
use WEEK1;
read all var{agem};
close WEEK1;

alpha=0.05;
Ybar=mean(agem);
s=var(agem);
n=nrow(agem);
qT=quantile('t',alpha/2,n-1);
UPL=Ybar-qT*sqrt((n+1)*s/n);
LPL=Ybar+qT*sqrt((n+1)*s/n);
A=Ybar||LPL||UPL;

create DATA from A[colname={'mean' 'LPL' 'UPL'}]; 
append from A;       
close DATA;
quit;

*W1-2 BOXCOX transformation;
PROC UNIVARIATE DATA=WEEK1 NORMALTEST;
VAR AGEM;
HISTOGRAM AGEM /NORMAL;
RUN;

DATA WEEK1BOXCOX;
 	SET WEEK1;
 	AGEMMINUS2= (-1/2)*((AGEM)**-2-1);
 	AGEMMINUS1= (-1)*((AGEM)**-1-1);
 	AGEMMINUS12= (-2)*((AGEM)**-(0.5)-1);
 	AGEM0= log(AGEM);
 	AGEM13= (3)*((AGEM)**(1/3)-1);
 	AGEM12= (2)*((AGEM)**(1/2)-1);
 	AGEM2= (0.5)*((AGEM)**(2)-1);
RUN;

proc univariate data=WEEK1BOXCOX noprint;
  	histogram AGEMMINUS2 /normal;
 	histogram AGEMMINUS1 /normal;
 	histogram AGEMMINUS12 /normal;
 	histogram AGEM0 /normal;
 	histogram AGEM13 /normal;
 	histogram AGEM12 /normal;
	histogram AGEM2 /normal;
run;

proc univariate data=WEEK1BOXCOX;
  	histogram AGEM/normal;
 run;


*Week 1-2 PI variance;
data agem_PI_var;
	set agem_PI;
	N_new = 50;
	F_U = QUANTILE('F', 1-0.01/2, N_new -1, N-1);
	F_L = QUANTILE('F',   0.01/2, N_new -1, N-1);
	LPL_var = F_L*std**2;
	UPL_var = F_U*std**2;
run;

proc print DATA=agem_PI_var;
	var LPL_var UPL_var;
run;

PROC IML;
use WEEK1;
read all var{agem};
close WEEK1;

alpha=0.01;
Ybar=mean(agem);
s=var(agem);
n=nrow(agem);
m=50;

qFL=quantile('F',alpha/2,m-1,n-1);
qFU=quantile('F',1-alpha/2,m-1,n-1);
LPL=s*qFL;
UPL=s*qFU;
A=s||LPL||UPL;

create PIvar from A[colname={'var' 'LPL' 'UPL'}]; 
append from A;       
close PIvar;
quit;

