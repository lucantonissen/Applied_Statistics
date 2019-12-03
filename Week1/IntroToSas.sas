*Read in data manually;
DATA color; 
   INPUT region eyes$ hair$ count@@; 
   DATALINES; 
	1  blue  fair   23 1 blue  red  7
	1  blue  medium 24 1 blue  dark 11
	1  green fair   19 1 green red  7 
 ;
RUN; 

*Export .csv file to specified folder;
PROC export data = color
outfile=  "/folders/myfolders/color.csv"
dbms= csv
replace;
run;

*Import .csv file from specified folder;
PROC IMPORT OUT= WORK.color2
  DATAFILE= "/folders/myfolders/color.csv"
    DBMS=CSV REPLACE;
    GETNAMES=YES;
RUN;

*Set library;
LIBNAME SASDATA "/folders/myfolders";
*Save .sas7bdat file to library;
DATA SASDATA.SAVEDATA;
	SET COLOR;
RUN;

*Import .sas7bdat file from library;
DATA  Color3;
	SET SASDATA.SAVEDATA;
RUN;

PROC SORT DATA=COLOR OUT=SORT_COLOR;
	BY COUNT EYES;
RUN;

PROC PRINT DATA=SORT_COLOR;
 	ID  REGION EYES COUNT;
RUN; 

*Summary statistics + create output dataset;
PROC MEANS DATA=color MEAN MEDIAN STD KURT SKEW p10;
	VAR count;
	BY  eyes;
	OUTPUT OUT = WORK.OUTPUT;
RUN;

*Alternative way for summary statistics + compute percentiles;
PROC UNIVARIATE DATA=COLOR;
   VAR COUNT;
   OUTPUT OUT=outdata PCTLPTS = 10 PCTLPRE = PERCENT;
RUN;


*Data visualization;
PROC UNIVARIATE DATA = COLOR;
   VAR COUNT;
   HISTOGRAM COUNT/NORMAL;                     	QQPLOT    COUNT/NORMAL;
   PPPLOT 	 COUNT/NORMAL;
RUN;


*Subsets example;
DATA SUBDATA;
	SET COLOR;
	WHERE Hair ^= 'red';
	IF COUNT > 30 OR EYES = 'green' 	   	THEN DELETE;
	ELSE NEWCOUNT=COUNT**2;
RUN;

DATA SUBSET1;
	SET COLOR;
	KEEP REGION COUNT;
RUN;

DATA SUBSET2;
	SET COLOR;
	DROP HAIR;
RUN;

*Combine example;
DATA COMBINE;
	SET SUBSET1 SUBSET2;
  RUN;

*MERGE example;

  DATA CLASS;
   INPUT Name $ Year $  
               Major $ ;
   DATALINES;
   Jennifer         first
    Tom              third     Theater
    Elissa           fourth    Math
    Rachel          first       Math
 ; 
RUN;

        /* Sort data*/
PROC SORT DATA= CLASS;
            BY Name;
RUN;

         /* Print data*/
PROC PRINT DATA= CLASS;
RUN;

DATA TIME;
   INPUT Name$ DATE$  TIME
            ROOM@@;
    DATALINES; 
    Jennifer  14sep2000  10  103 
      Rachel  14sep2000  10  103 TOM 
    14sep2000 11 207 Elissa 15sep2000  
    10  105
   ; 
 RUN;

        /* Sort data*/
PROC SORT DATA= TIME;
            BY Name;
RUN;
        
           /* Print data*/
PROC PRINT DATA= TIME;
RUN;

DATA CL_TIME;
  MERGE CLASS (KEEP = NAME  YEAR) 
        TIME  (DROP = DATE		 RENAME=(ROOM= CLASS_NUMBER));
 	BY NAME;
RUN;

