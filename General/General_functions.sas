LIBNAME SASDATA "/folders/myfolders/Applied_Statistics/Data";

DATA WEEK1;
	set SASDATA.IVF;
	where PER=4;
RUN;


/* View data */
PROC PRINT data=WEEK1 NOOBS;


/* Check normal distribution */
/* ods select HISTOGRAM, QQPLOT, PPPLOT; */
PROC UNIVARIATE data=WEEK1;
   VAR AGEM;
   HISTOGRAM AGEM/NORMAL;
   QQPLOT    AGEM/NORMAL;
   PPPLOT 	 AGEM/NORMAL;
RUN;