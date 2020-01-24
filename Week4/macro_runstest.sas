/*©R.A.J.Post 2017*/
Title ' Large sample approximation Runs test';	

DATA Lecture;
	INPUT  BATCH CONTENT @@;
	DATALINES;
	
1 	102   2   104   3 	102  4 	97 	5 	95  
6 	106   7   103   8 	98 	 9 	96 	10 	97	
;
RUN;


PROC SORT DATA=Lecture;
BY BATCH;
RUN;


%MACRO RUNSCUC(data=,var=,alpha=);
PROC IML;
use &data;
read all var {&var};
close &data;

X=&var;
n=nROW(X);
MED=median(X);

XC=X;
DO i=1 to n by 1;
	IF (XC[i] >= MED) then XC[i]=1;
	ELSE XC[i]=0;
END;

n1C=sum(XC);
n2C=n-n1C;

RC=1;
DO i=2 to n by 1;
	if(XC[i] ^= XC[i-1]) then RC=RC+1;
END;

MUC=1+(2*n1C*n2C)/(n1C+n2C);
VARC=2*n1C*n2C*(2*n1C*n2C-n1C-n2C)/((n1C+n2C-1)*(n1C+n2C)**2);

SC=(RC-MUC)/SQRT(VARC);
TC=QUANTILE('NORMAL',&alpha/2);
TCU=QUANTILE('NORMAL',1-&alpha/2);
PC=(1-CDF('NORMAL',abs(SC)))*2;

XUC=REPEAT(0,n-1,1);
TIES=0;
DO i=1 to (n-1) by 1;
	IF (X[i+1] > X[i]) then XUC[i]=1;
	IF (X[i+1] < X[i]) then XUC[i]=0;
	IF (X[i+1] = X[i]) then XUC[i]=XUC[i-1];
	IF (X[i+1] = X[i]) then TIES=TIES+1;
END;

RUC=1;
DO i=2 to (n-1) by 1;
	if(XUC[i] ^= XUC[i-1]) then RUC=RUC+1;
END;

MUUC=(2*(n-TIES)-1)/3;
VARUC=(16*(n-TIES)-29)/90;

SUC=(RUC-MUUC)/SQRT(VARUC);
TUC=QUANTILE('NORMAL',&alpha/2);
TUCU=QUANTILE('NORMAL',1-&alpha/2);
PUC=(1-CDF('NORMAL',abs(SUC)))*2;

PRINT("Median based (conditional) runs test");
PRINT(RC||MUC||sqrt(VARC)||PC||SC||TC||TCU||n);
PRINT("(unconditional) runst test for serial randomness");
PRINT(TIES);
PRINT(RUC||MUUC||sqrt(VARUC)||PUC||SUC||TUC||TUCU||(n-TIES));
quit;
%MEND;

%RUNSCUC(data=LECTURE,var=CONTENT,alpha=0.05);
