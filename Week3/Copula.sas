LIBNAME SASDATA "/folders/myfolders/Applied_Statistics/Data";

DATA IVFN;
SET SASDATA.IVF;
IMP = IMP + (ranuni(1)-0.5);
KEEP ID PER IMP;
run;

PROC TRANSPOSE OUT=WIDE_IVFN(DROP = _NAME_ _LABEL_) DATA=IVFN PREFIX=IMP; 
BY ID;
ID PER;
VAR IMP;
RUN; 

DATA WIDE_IVFN;
SET WIDE_IVFN;
if cmiss(of _all_) then delete;
run;

proc rank data=WIDE_IVFN out=ranked_a;
      var IMP4 IMP10;
      ranks rank_IMP4 rank_IMP10;
   run;

proc means data=ranked_a N;
var rank_IMP4 rank_IMP10;
run;

data marginals;
set ranked_a;
U_IMP4=rank_IMP4/237;
U_IMP10=rank_IMP10/237;
run;
 

*Macro's;
%Macro SIM_GC(rho=, nsim=, seed=);
proc iml;

use marginals;
read all var{U_IMP4};
close marginals;

U=U_IMP4;

call streaminit(&seed);
rho=&rho;

do i=1 to &nsim by 1;
U1=U[i];
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

%Macro SIM_Gum(alpha=, nsim=, seed=);
proc iml;
call streaminit(&seed);
alpha=&alpha;

use marginals;
read all var{U_IMP4};
close marginals;

U=U_IMP4;

do i=1 to &nsim by 1;
U1=U[i];
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


%Macro SIM_Clay(alpha=, nsim=, seed=);
proc iml;
call streaminit(&seed);
alpha=&alpha;

use marginals;
read all var{U_IMP4};
close marginals;

U=U_IMP4;

do i=1 to &nsim by 1;
U1=U[i];
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


%Macro SIM_Frk(alpha=, nsim=, seed=);
proc iml;
call streaminit(&seed);
alpha=&alpha;

use marginals;
read all var{U_IMP4};
close marginals;

U=U_IMP4;

do i=1 to &nsim by 1;
U1=U[i];
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


%Macro SIM_FGM(alpha=, nsim=, seed=);
proc iml;
call streaminit(&seed);
alpha=&alpha;

use marginals;
read all var{U_IMP4};
close marginals;

U=U_IMP4;

do i=1 to &nsim by 1;
U1=U[i];
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

*Estimate parameters;
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


CPAR=SGum||SClay||SFrk||SGau||SFGM;

create EstKendall from CPAR [colname={'Gumbel alpha','Clayton alpha','Frank alpha','Gaussian rho','FGM alpha'}]; 
append from CPAR;       
close EstKendall;
quit;

%mend KendallTau;

proc corr data=marginals kendall spearman;
var U_IMP4 U_IMP10;
run;

proc corr data=marginals kendall spearman;
var IMP4 IMP10;
run;

%KendallTau(tau=0.18759);
%SpearmanRho(rho=0.27997);

*Run;

%SIM_GC(nsim=236, rho=0.293, seed=12345);
%SIM_Gum(nsim=236, alpha=1.235, seed=12345);
%SIM_Clay(nsim=236, alpha=0.47, seed=12345);
%SIM_Frk(nsim=236, alpha=1.76, seed=12345);
%SIM_FGM(nsim=236, alpha=0.85, seed=12345);

proc sgplot data=marginals aspect=1;
title "Marginals IMP4 and IMP10";
scatter x=U_IMP4 y=U_IMP10 / markerattrs=(symbol=circlefilled color='red' size=8 );
run;

proc sgplot data=Gumc aspect=1;
title "Simulated Gumbel copula";
scatter x=X y=Y / markerattrs=(symbol=circlefilled color='darkgreen' size=8 );
run;
proc sgplot data=Frkc aspect=1;
title "Simulated Frank's copula";
scatter x=X y=Y / markerattrs=(symbol=circlefilled color='cyan' size=8 );
run;

proc sgplot data=Gc aspect=1;
title "Simulated Gaussian copula";
scatter x=X y=Y / markerattrs=(symbol=circlefilled color='blue' size=8 );
run;

proc sgplot data=Cc aspect=1;
title "Simulated Clayton copula";
scatter x=X y=Y / markerattrs=(symbol=circlefilled color='navy' size=8 );
run;

proc sgplot data=Fgmc aspect=1;
title "Simulated FGM copula ";
scatter x=X y=Y / markerattrs=(symbol=circlefilled color='black' size=8 );
run;

*Simulate data using empirical F;
proc iml;
use WIDE_IVFN;
read all var{IMP4 IMP10};
close WIDE_IVFN;

use GUMC;
read all var{X Y};
close GUMC;

call qntl(IVFXgum, IMP4, X);
call qntl(IVFYgum, IMP10, Y);

use Cc;
read all var{X Y};
close Cc;

call qntl(IVFXclay, IMP4, X);
call qntl(IVFYclay, IMP10, Y);

use Frkc;
read all var{X Y};
close Frkc;

call qntl(IVFXfrk, IMP4, X);
call qntl(IVFYfrk, IMP10, Y);

use Gc;
read all var{X Y};
close Gc;

call qntl(IVFXgau, IMP4, X);
call qntl(IVFYgau, IMP10, Y);

use Fgmc;
read all var{X Y};
close Fgmc;

call qntl(IVFXfgm, IMP4, X);
call qntl(IVFYfgm, IMP10, Y);


plot=IVFXgum||IVFYgum||IVFXclay||IVFYclay||IVFXfrk||IVFYfrk||IVFXgau||IVFYgau||IVFXfgm||IVFYfgm;

create SimIVF from plot [colname={'Xgum','Ygum','Xclay','Yclay','Xfrk','Yfrk','Xgau','Ygau','Xfgm','Yfgm'}]; 
append from plot;       
close SimIVF;
quit;

proc sgplot data=marginals aspect=1;
title "IMP4 and IMP10";
scatter x=IMP4 y=IMP10 / markerattrs=(symbol=circlefilled color='red' size=8 );
run;

proc sgplot data=Simivf aspect=1;
title "Simulated Gumbel copula";
scatter x=Xgum y=Ygum / markerattrs=(symbol=circlefilled color='darkgreen' size=8 );
run;
proc sgplot data=Simivf aspect=1;
title "Simulated Frank's copula";
scatter x=Xfrk y=Yfrk / markerattrs=(symbol=circlefilled color='cyan' size=8 );
run;

proc sgplot data=Simivf aspect=1;
title "Simulated Gaussian copula";
scatter x=Xgau y=Ygau / markerattrs=(symbol=circlefilled color='blue' size=8 );
run;

proc sgplot data=Simivf aspect=1;
title "Simulated Clayton copula";
scatter x=Xclay y=Yclay / markerattrs=(symbol=circlefilled color='navy' size=8 );
run;

proc sgplot data=Simivf aspect=1;
title "Simulated FGM copula ";
scatter x=Xfgm y=Yfgm / markerattrs=(symbol=circlefilled color='black' size=8 );
run;
