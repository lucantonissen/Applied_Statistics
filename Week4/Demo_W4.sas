libname SAS "/folders/myfolders/";

data RCT;
	set SAS.RCT;
run;

*a - create wide dataset;
proc sort data=rct;
by ID;
run;

PROC TRANSPOSE OUT=WIDE_RCT(DROP = _NAME_) DATA=RCT PREFIX=RESP; 
BY ID;
ID TIME;
VAR RESP;
RUN; 

*b/d;
proc corr data=WIDE_RCT kendall spearman pearson fisher(biasadj=no);
VAR RESP1 RESP2;
run;

*0.47800;

*e;
proc iml;
alpha=(0.35941)*9/2;
print(alpha);
run;

*f/g;

%SpearmanRho(rho=0.50477);
%KendallTau(tau=0.35941);

*h;

*Create Uniform marginals that are non-discrete;
DATA WIDE_RCT;
SET WIDE_RCT;
RESP1=RESP1 + 0.1*(ranuni(1)-0.5);
RESP2=RESP2 + 0.1*(ranuni(1)-0.5);
run;


DATA WIDE_RCT;
SET WIDE_RCT;
if cmiss(of RESP1 RESP2) then delete;
run;

proc rank data=WIDE_RCT out=ranked3;
      var RESP1 RESP2;
      ranks rank_RESP1 rank_RESP2;
 run;

proc means data=ranked3 N;
var  rank_RESP1 rank_RESP1;
run;

data marginals;
set ranked3;
U_RESP1=rank_RESP1/(712+1);
U_RESP2=rank_RESP2/(712+1);
run;

 proc sgplot data=marginals aspect=1;
title "Marginals RESP1 and RESP2 adding little noise";
scatter x=U_RESP1 y=U_RESP2 / markerattrs=(symbol=circlefilled color='red' size=8 );
run;

%SIM_Clay(nsim=700, alpha=1.1, seed=12345);
proc sgplot data=Cc aspect=1;
title "Simulated Clayton copula";
scatter x=X y=Y / markerattrs=(symbol=circlefilled color='navy' size=8 );
run;

%SIM_Frk(nsim=700, alpha=3.5, seed=12345);

proc sgplot data=Frkc aspect=1;
title "Simulated Frank's copula";
scatter x=X y=Y / markerattrs=(symbol=circlefilled color='cyan' size=8 );
run;
