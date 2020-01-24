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

