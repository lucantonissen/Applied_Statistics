LIBNAME SAS "D:\Documents\PhD\Teaching\Applied statistics\2018-2019\SAS code";


DATA RCT;
 SET SAS.Cellsaver;
RUN;

DATA RCT1;
SET RCT;
where TIME=3;
run;

PROC MIXED DATA=RCT1 METHOD=TYPE3 CL;
	CLASS TRT;
	MODEL RESP = TRT /SOLUTION CL;
RUN;

proc mixed data=RCT method=TYPE3 cl;
	CLASS TIME;
	MODEL RESP = /SOLUTION cl;
	RANDOM TIME /SOLUTION;
run;

proc freq DATA=RCT;
	table RESP;
run;

proc npar1way DATA=RCT1 wilcoxon;
	class TRT;
	var RESP;
run;

*Week 4-2;

data CU;
	input MEASUREMENT@@;
	datalines;
102 104 101 97 95 106 103 98 96 99
;
run;

proc reg data=CU plots=none;
   model MEASUREMENT = /dwProb;
run;

proc timeseries data=CU OUTCORR=ACF;
Var MEASUREMENT;
run;

proc sort data=rct;
by center;
run;

proc sgplot data=rct;
where resp<12 and resp>2;
scatter X=center Y=RESP;
run;

%RUNSCUC(data=rct,var=RESP,alpha=0.05);

proc sort data=rct;
by time ID;
run;

proc mixed data=RCT method=TYPE3 cl;
	where ID<100;
	CLASS TIME;
	MODEL RESP = /SOLUTION cl OUTPM=RM OUTP=RC;
	RANDOM TIME /SOLUTION;
run;

proc means data=rct;
where ID<100;
var resp;
by time;
run;

DATA RM;
SET RM;
RESIDM=RESID;
PredM=Pred;
DROP RESID Pred;
RUN;

DATA PREDMERGE;
MERGE RM RC;
BY TIME ID;
RUN;

DATA PREDMERGE;
SET PREDMERGE;
EBLUP=PRED-PREDM;
RUN;


proc sort data=predmerge;
by EBLUP;
run;

Title 'Marginal residuals';	
proc sgplot data=predmerge;
scatter X=EBLUP Y=RESIDM;
run;

%RUNSCUC(data=predmerge,var=RESIDM,alpha=0.05);


Title 'Conditional residuals';	
proc sgplot data=predmerge;
scatter X=EBLUP Y=RESID;
run;
%RUNSCUC(data=predmerge,var=RESID,alpha=0.05);


