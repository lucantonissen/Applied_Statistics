
DATA IVF;
 SET SAS.IVF;
  KEEP PER IMP ID;
RUN;

PROC TRANSPOSE OUT=WIDE_IVF(DROP = _NAME_ _LABEL_) DATA=IVF PREFIX=IMP; 
BY ID;
ID PER;
VAR IMP;
RUN; 

ods graphics on;
proc corr data=WIDE_IVF  plots=scatter(ellipse=none) spearman kendall;
var IMP10 IMP4;
run;
ods graphics off;


data GRADES;
	input ID PRE EXAM @@;
	datalines;
1 8 8 2 7 6.9 3 7.3 7 4 6.5 6.4 5 7.5 6.8
6 7.1 7 7 7.5 6.9 8 8.3 8.4 9 7 7.3 10 7.2 6.5
11 7.1 7 12 7 7.4 13 6.9 7.4 14 6.8 6.5 15 7.5 7.2
	;
run;
proc corr DATA=GRADES pearson;
	var PRE EXAM;
run;

proc corr data=GRADES pearson fisher(biasadj=no);
	var PRE EXAM;
run;

data GRADES;
	set GRADES;
	PRE_A = (PRE>7);
	EXAM_A = (EXAM>7);
run;

proc corr data=GRADES pearson;
	var PRE_A EXAM_A;
run;

proc freq data=GRADES;
	tables PRE_A*EXAM_A/chisq;
run;

proc corr data=GRADES spearman;
	var PRE EXAM;
run;

proc corr data=GRADES kendall;
	var PRE EXAM;
run;

**Copula simulations;
%Macro SIM_GC(rho=, nsim=, seed=);
proc iml;
call streaminit(&seed);
rho=&rho;

do i=1 to &nsim by 1;
U1=rand('Uniform');
U2=rand('Uniform');

start Func(x) global(U1,U2,rho);
return(CDF('Normal',quantile('NORMAL', x),rho*quantile('NORMAL',U1),(1-rho**2))-U2);
finish;

intervals = {0.00001 0.99999};        
U2C = froot("Func", intervals);

X=X//U1;
Y=Y//U2C;
YI=YI//U2;
end;

Total=X||Y||YI;

create GC from Total [colname={'X','Y','YI'}]; 
append from Total;       
close GC;
quit;
%mend SIM_GC;

%SIM_GC(nsim=1000, rho=0.9, seed=12345);

proc sgplot data=GC;
scatter x=X y=Y;
run;

%Macro SIM_Gum(alpha=, nsim=, seed=);
proc iml;
call streaminit(&seed);
alpha=&alpha;

do i=1 to &nsim by 1;
U1=rand('Uniform');
U2=rand('Uniform');

start Func(x) global(U1,U2,alpha);
return(Exp(-((-Log(x))**alpha + (-Log(U1))**alpha)**(1/alpha))*((-Log(x))**alpha + (-Log(U1))**alpha)**(-1 + 1/alpha)*((-Log(U1))**(alpha-1))/U1-U2);
finish;

intervals = {0.00001 1};        
U2C = froot("Func", intervals);

X=X//U1;
Y=Y//U2C;
YI=YI//U2;
end;

Total=X||Y||YI;

create GumC from Total [colname={'X','Y','YI'}]; 
append from Total;       
close GumC;
quit;
%mend SIM_Gum;

%SIM_Gum(nsim=1000, alpha=5, seed=12345);

proc sgplot data=GumC aspect=1;
scatter x=X y=Y;
run;



%Macro SIM_Clay(alpha=, nsim=, seed=);
proc iml;
call streaminit(&seed);
alpha=&alpha;

do i=1 to &nsim by 1;
U1=rand('Uniform');
U2=rand('Uniform');

start Func(x) global(U1,U2,alpha);
return(U1**(-1 -alpha)*(x**(-alpha) + U1**(-alpha)-1)**(-1 - 1/alpha)-U2);
finish;

intervals = {0.001 1};        
U2C = froot("Func", intervals);

X=X//U1;
Y=Y//U2C;
end;

Total=X||Y;

create CC from Total [colname={'X','Y'}]; 
append from Total;       
close CC;
quit;
%mend SIM_Clay;

%SIM_Clay(nsim=1000, alpha=5, seed=12345);

proc sgplot data=CC;
scatter x=X y=Y;
run;

%Macro SIM_Frk(alpha=, nsim=, seed=);
proc iml;
call streaminit(&seed);
alpha=&alpha;

do i=1 to &nsim by 1;
U1=rand('Uniform');
U2=rand('Uniform');

start Func(x) global(U1,U2,alpha);
return((Exp(alpha)*(-1 + Exp(alpha*x)))/(-Exp(alpha) + Exp(alpha*(1+x)) - Exp(alpha*(U1+x)) + Exp(alpha*(1 + U1)))-U2);
finish;

intervals = {0.00001 1};        
U2C = froot("Func", intervals);

X=X//U1;
Y=Y//U2C;
end;

Total=X||Y;

create FrkC from Total [colname={'X','Y'}]; 
append from Total;       
close FrkC;
quit;
%mend SIM_Frk;


%SIM_Frk(nsim=1000, alpha=-5, seed=12345);

proc sgplot data=FrkC;
scatter x=X y=Y;
run;


%Macro SIM_FGM(alpha=, nsim=, seed=);
proc iml;
call streaminit(&seed);
alpha=&alpha;
do i=1 to &nsim by 1;
U1=rand('Uniform');
U2=rand('Uniform');

start Func(x) global(U1,U2,alpha);
return(x*(1 + alpha*(1 - x)*(1 - U1)) - alpha*(1 - x)*x*U1-U2);
finish;

intervals = {0.00001 1};        
U2C = froot("Func", intervals);

X=X//U1;
Y=Y//U2C;
YI=YI//U2;
end;

Total=X||Y||YI;

create FGMC from Total [colname={'X','Y','YI'}]; 
append from Total;       
close FGMC;
quit;
%mend SIM_FGM;

%SIM_FGM(nsim=240, alpha=0.5, seed=12345);

title "Simulated FGM copula with Uniform Marginals";
proc sgplot data=FGMC aspect=1;
scatter x=X y=Y;
run;

*Estimate copula parameters using Spearman rho;




%Macro SpearmanRho(rho=);
proc iml;
pi = constant("pi");
tau=&rho;


start innerGum(y) global(alpha, x); 
   return(Exp(-((-Log(x))**alpha + (-Log(y))**alpha)**(1/alpha)));
finish; 

start outerGum(par) global(x,alpha); 
	x=par;
   yinterval = 0 || 1;
   /** evaluate inner integral for the parameter value, a=x **/ 
   call quad(w, "innerGum", yinterval);
   return (w);
finish; 

start finalGum(param) global(alpha, tau);
alpha=param;
xinterval= {0 1};
call quad(v, "outerGum", xinterval); /** outer integral **/ 
return(12*v-(3+tau));
finish;


/*
t = do(1, 100, 1);          
F = j(1, ncol(t));
do i = 1 to ncol(t);
   F[i] = finalGum(t[i]);      
end;
title "Integral - 3+ tau";
call Series(t, F) grid={x y} label={"x" "F(x)"}
                  other="refline 0 / axis=y";  
*/

intervalsGum = {1 100};        
SGum = froot("finalGum", intervalsGum);
print(SGum);

start innerClay(y) global(alpha, x); 
	return((x**(-alpha)+y**(-alpha)-1)**(-1/alpha));
finish; 

start outerClay(par) global(x, alpha); 
	x=par;

if(alpha>0) then yinterval= 0||1;
else yinterval= (1-x**(-alpha))**(-1/alpha)||1;
   /** evaluate inner integral for the parameter value, a=x **/ 
   call quad(w, "innerClay", yinterval);
   return (w);
finish; 

start finalClay(param) global(alpha, tau);
alpha=param;
xinterval= {0 1};
call quad(v, "outerClay", xinterval); /** outer integral **/ 
return(12*v-(3+tau));
finish;

/*
t = do(-1, 3, 0.11);          
F = j(1, ncol(t));
do i = 1 to ncol(t);
   F[i] = finalClay(t[i]);      
   print(finalClay(t[i]));
end;
title "Integral - 3+ tau";
call Series(t, F) grid={x y} label={"x" "F(x)"}
                  other="refline 0 / axis=y";  
*/
                 
intervalsClay = {-1 10};        
SClay = froot("finalClay", intervalsClay);
print(SClay);

SGau=2*sin(pi*tau/6);
print(SGau);

SFGM=3*tau;
print(SFGM);

start innerFrk(y) global(alpha, x); 
return(-(1/alpha)*Log(1+(Exp(-alpha*x)-1)*(Exp(-alpha*y)-1)/(Exp(-alpha)-1)));
finish; 

start outerFrk(par) global(x, alpha); 
	x=par;
   yinterval = 0 || 1;
   /** evaluate inner integral for the parameter value, a=x **/ 
   call quad(w, "innerFrk", yinterval);
   return (w);
finish; 

start finalFrk(param) global(alpha, tau);
alpha=param;
xinterval= {0 1};
call quad(v, "outerFrk", xinterval); /** outer integral **/ 
return(12*v-(3+tau));
finish;

/*

t = do(-10, 10, 1.1);          
F = j(1, ncol(t));
do i = 1 to ncol(t);
   F[i] = finalFrk(t[i]);      
  
end;
title "Integral - 3+ tau";
call Series(t, F) grid={x y} label={"x" "F(x)"}
                  other="refline 0 / axis=y";  
*/
                 
intervalsFrk = {-30 30};        
SFrk = froot("finalFrk", intervalsFrk);
print(SFrk);




CPAR=SGum||SClay||SFrk||SGau||SFGM;

create EstSpearman from CPAR [colname={'Gumbel alpha','Clayton alpha','Frank alpha','Gaussian rho','FGM alpha'}]; 
append from CPAR;       
close EstSpearman;
quit;

%mend SpearmanRho;

%SpearmanRho(rho=0.3);

%Macro KendallTau(tau=);
proc iml;
pi = constant("pi");
tau=&tau;

SGum=1/(1-tau);
print(SGum);
SClay=2*tau/(1-tau);
print(SClay);
SGau=sin(pi*tau/2);
print(SGau);
SFGM=9*(tau/2);
print(SFGM);

start D(y);
return(y/(Exp(y)-1));
finish;

*IF alpha>0 / tau>0;
start FC(x) global(tau);
dinterval=0||x;
call quad(w, "D", dinterval);
return(1-(4/x)*(1-(1/x)*w)-tau);
finish;

intervals = {0.00001 20};        
SFrk = froot("FC", intervals);
print(SFrk);

/*

*IF alpha<0 / tau<0;
start FC(x) global(tau);
dinterval=0||-x;
call quad(w, "D", dinterval);
return(1-(4/x)*(1+(1/x)*w+0.5*x)-tau);
finish;

intervals = {-10 -0.00001};        
SFrk = froot("FC", intervals);
print(SFrk);

*/


CPAR=SGum||SClay||SFrk||SGau||SFGM;

create EstKendall from CPAR [colname={'Gumbel alpha','Clayton alpha','Frank alpha','Gaussian rho','FGM alpha'}]; 
append from CPAR;       
close EstKendall;
quit;

%mend KendallTau;

%KendallTau(tau=0.3);

*Estimate parameter Frank's copula using Kendall's tau;
proc iml;
tau=0.3;

start D(y);
return(y/(Exp(y)-1));
finish;

*IF alpha>0 / tau>0;
start FC(x) global(tau);
dinterval=0||x;
call quad(w, "D", dinterval);
return(1-(4/x)*(1-(1/x)*w)-tau);
finish;

intervals = {0.00001 20};        
alphaFrk = froot("FC", intervals);
print(alphaFrk);

/*

*IF alpha<0 / tau<0;
start FC(x) global(tau);
dinterval=0||-x;
call quad(w, "D", dinterval);
return(1-(4/x)*(1+(1/x)*w+0.5*x)-tau);
finish;

intervals = {-10 -0.00001};        
alphaFrk = froot("FC", intervals);
print(alphaFrk);

*/
quit;



 *W3-2: Paired samples;

data GRADES;
	input ID PRE EXAM @@;
	datalines;
1 8 8 2 7 6.9 3 7.3 7 4 6.5 6.4 5 7.5 6.8
6 7.1 7 7 7.5 6.9 8 8.3 8.4 9 7 7.3 10 7.2 6.5
11 7.1 7 12 7 7.4 13 6.9 7.4 14 6.8 6.5 15 7.5 7.2
	;
run;


data GRADES;
	set GRADES;
	ZD=EXAM-PRE;
	ZL=log(EXAM)-log(PRE);
	ZR=EXAM/PRE-1;
run;

PROC UNIVARIATE DATA=GRADES NORMAL;
	VAR ZD ZL ZR;
	HISTOGRAM ZD ZL ZR;
RUN;

data GRADES;
	set GRADES;
	B1 = (PRE>7);
	B2 = (EXAM>7);
run;


PROC FREQ DATA=GRADES;
	TABLES B1*B2; EXACT MCNEM;
RUN;


