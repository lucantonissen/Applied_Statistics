libname SAS "/folders/myfolders/";

data IVF_DATASET;
	set SAS.IVF;
	where PER=4;
	drop IMP PER AGE;
run;

data Q22;
	set IVF_DATASET;
	WHERE ID LE 100;
run;

*(a);
ods output WilcoxonScores=WRS12 (keep= Class N SumOfScores);
proc npar1way data=Q22 correct=NO;
	WHERE TRT=1 OR TRT=2;
	class TRT;
	var BW;
	exact wilcoxon /MC;
run;

*(b);

PROC IML;
use WRS02;
read all var{N SumOfScores};
close WRS02;
G={0 , 1};
U=SumOfScores-N#(N+1)/2;
P=U/prod(N);

A=G||N||U||P;


create MWU02 from A [colname={'Group' 'N' 'U' 'P'}]; 
append from A;       
close MWU02;
quit;

*c;
*U_control=277;
*P(Y_TRT=1<Y_TRT=0)=277/(13*34)=0.6266968326;

*P(Y_TRT=2<Y_TRT=0)=538.5/(34*23)=0.6886189258;

*d;

proc npar1way data=Q22 correct=NO;
	WHERE TRT=1 OR TRT=2;
	class TRT;
	var BW;
	exact ks /MC;
run;

*test statistic = 0.3620 // exact MC p-value =0.1276 // approx p-value 0.1700
TRT=0 / TRT=1;

*e;

data Q22;
set Q22;
BW2= (0.5)*((BW)**(2)-1);
run;

ods select histogram;
proc univariate data=Q22;
histogram BW2/normal;
run;

proc ttest data=Q22;
where TRT=0 OR TRT=1;
class TRT;
var BW2;
run;

proc ttest data=Q22;
where TRT=0 OR TRT=2;
class TRT;
var BW2;
run;


