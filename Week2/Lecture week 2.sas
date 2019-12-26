
*W2-1;
DATA WEEK2;
 SET IVF;
 WHERE TRT=0 AND PER=4;
 IMPSHIFT=IMP-10;
 label IMPSHIFT = "Neurological score - 10";
 IF AGEM>35 then IMPSO=IMP-10;
 else IMPSO=IMP;
 label IMPSO = "Neurological score AGEM>35 - 10";
 IMPU = 2*IMP-79.1363636;
 label IMPU = "Neurological score variance*4";
 RUN;

proc sgplot data=WEEK2;
   histogram IMP / fillattrs=graphdata1 transparency=0.7 binwidth=2;
   histogram IMPSO / fillattrs=graphdata2 transparency=0.5 binwidth=2;
   yaxis grid; xaxis values=(60 to 90 by 1);
run;



data CU;
	input BATCH$ MEASUREMENT@@;
	datalines;
1 102 1 104 1 102 1 97 1 99
1 101 1 103 1 98 1 96 1 97
2 99 2 97 2 99 2 100 2 99
2 96 2 99 2 98 2 97 2 98
	;
run;

proc iml;
use CU;
read all var{BATCH MEASUREMENT};
close CU;

x=MEASUREMENT;
g=BATCH;

n1=nrow(x[loc(g='1')]);
n2=nrow(x[loc(g='2')]);
s1=var(x[loc(g='1')]);
s2=var(x[loc(g='2')]);
sp=((n1-1)*s1+(n2-1)*s2)/(n1+n2-2);
C2=(1+(1/(n1-1)+1/(n2-1)-1/(n1+n2-2))/3);
gB=(n1+n2)*(sum((x[loc(g='1')]-median(x[loc(g='1')]))##4)+sum((x[loc(g='2')]-median(x[loc(g='2')]))##4))/((sum((x[loc(g='1')]-mean(x[loc(g='1')]))##2)+sum((x[loc(g='2')]-mean(x[loc(g='2')]))##2))**2);
C2a=C2*(1+(gB-3)/2);
Ba=((n1+n2-2)*log(sp)-(n1-1)*log(s1)-(n2-1)*log(s2))/C2a;
B=((n1+n2-2)*log(sp)-(n1-1)*log(s1)-(n2-1)*log(s2))/C2;
pa=1-CDF('CHISQ',Ba,1);
p=1-CDF('CHISQ',B,1);

FinalBC= Ba||pa||B||p;

create finalBc from FinalBC [colname={'Ba','p','B','pB'}]; 
append from finalBC;       
close finalBc;
quit;


proc freq data=cu;
table batch;
run;

proc ttest data=CU;
	class BATCH;
	var MEASUREMENT;
run;

data Fval;
 qU2=quantile('F',0.975,9,9);
 qL2=quantile('F',0.025,9,9);
 qU=quantile('F',0.95,9,9);
 qL=quantile('F',0.05,9,9);
 pU=1-cdf('F',5.360,9,9);
 pL=cdf('F',5.360,9,9);
 test=quantile('CHISQ',0.95,1);
run;

ods select Bartlett;
proc glm data=CU;
	class BATCH;
	model MEASUREMENT = BATCH;
	means BATCH / hovtest=BARTLETT;
run;

ods select HOVFTest;
proc glm data=CU;
	class BATCH;
	model MEASUREMENT = BATCH;
	means BATCH / hovtest=Levene;
run;

ods select HOVFTest;
proc glm data=CU;
	class BATCH;
	model MEASUREMENT = BATCH;
	means BATCH / hovtest=BF;
run;

ods output WilcoxonScores=WRS (keep= Class N SumOfScores);
proc npar1way data=CU correct=NO;
	class BATCH;
	var MEASUREMENT;
	exact wilcoxon;
run;

PROC IML;
use WRS;
read all var{N SumOfScores};
close WRS;
G={1 , 2};
U=SumOfScores-N#(N+1)/2;
P=U/prod(N);

A=G||N||U||P;


create MWU from A [colname={'Group' 'N' 'U' 'P'}]; 
append from A;       
close MWU;
quit;


proc npar1way data=CU correct=NO;
	class BATCH;
	var MEASUREMENT;
	exact wilcoxon /mc;
run;


proc npar1way DATA=CU;
	class BATCH;
	var MEASUREMENT;
	exact KS;
run;

*W2-2;
data CU;
	set CU;
	if MEASUREMENT<97 OR MEASUREMENT>103 then OSPECLIM=1;
	else OSPECLIM=0;
run;

proc freq data=CU;
	table BATCH*OSPECLIM / chisq;
	exact chisq;
RUN;

proc freq data=CU;
	table BATCH*OSPECLIM / chisq;
	exact fisher;
RUN;


proc format;
   value FIRSTFmt 1='Tea'
                0='Milk';
   value GUESSFmt 1='Tea'
                0='Milk';
run;
data LTT;
   input FIRST GUESS COUNT;
   datalines;
0 0  4
0 1  0
1 0  0
1 1  4
;
proc sort data=LTT;
   by descending FIRST descending GUESS;
run;

proc freq data=LTT order=data;
   format FIRST FIRSTFmt. GUESS GUESSFmt.;
   tables FIRST*GUESS / chisq;
   weight COUNT;
 run;

data RCT;
 set SAS.Cellsaver;
 where TIME=1;
 keep CENTER TRT;
run;

proc freq data=RCT;
	tables CENTER*TRT /chisq;
run;

