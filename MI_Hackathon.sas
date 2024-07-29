* Topic: SIBS Program
*
* Title: Analysis of SIBS Program Myocardial Infarction Dataset
*
* Program name: MI_Hackathon.sas
*
* Author: Lou Lindsley and Keely Grice
*
* Date created: 4/20/24;

/*STAGE 1a: Importing data and libname*/
libname sibs "/home/u62079390/SIBS Hackathon/Data";

proc import file="/home/u62079390/SIBS Hackathon/Data/hack_1.csv" 
	out=work.mi_median
	dbms=csv;*dataset with missing values imputed with the median;
run;
proc import file="/home/u62079390/SIBS Hackathon/Data/hack_NN.csv" 
	out=work.mi_NN
	dbms=csv; *dataset with missing values imputed via nearest neighbor;
run;

data mi_median; *outcome of interest: REC_IM (relapse of MI);
	set work.mi_median;
	R_AB_3_n2 = put(R_AB_3_n, 3.);
	zab_leg_042 = put(zab_leg_04, 3.);
/* 	REC_IM2 = put(REC_IM, 3.); *converting complication vars; */
/* 	OTEK_LANC2 = put(OTEK_LANC, 3.); */
/* 	FIBR_PREDS2 = put(FIBR_PREDS, 3.); */
	
	drop ID PREDS_TAH JELUD_TAH FIBR_JELUD 
		 A_V_BLOK RAZRIV DRESSLER ZSN P_IM_STEN LET_IS
		 R_AB_3_n zab_leg_04; *dropping unnecissary outcome vars;
	
	rename R_AB_3_n2=R_AB_3_n;
	rename zab_leg_042=zab_leg_04;
/* 	rename REC_IM2=REC_IM; */
/* 	rename OTEK_LANC2=OTEK_LANC; */
/* 	rename FIBR_PREDS2=FIBR_PREDS; */

	if (FIBR_PREDS=1) OR (OTEK_LANC=1) OR (REC_IM=1)
		then COMPLIC=1;
		else COMPLIC=0; *New outcome var that's 1 if >0 outcomes are 1;
	if (STENOK_AN > 2) AND (STENOK_AN < 3) then STENOK_AN = 1; *replacing messed up imputation w/ median;
run;
data mi_NN; 
	set work.mi_NN;
	if (FIBR_PREDS=1) OR (OTEK_LANC=1) OR (REC_IM=1)
	then COMPLIC=1;
	else COMPLIC=0; *New outcome var that's 1 if >0 outcomes are 1;
	
	R_AB_3_n2 = put(R_AB_3_n, 3.);
	zab_leg_042 = put(zab_leg_04, 3.);
/* 	REC_IM2 = put(REC_IM, 3.); *converting complication vars; */
/* 	OTEK_LANC2 = put(OTEK_LANC, 3.); */
/* 	FIBR_PREDS2 = put(FIBR_PREDS, 3.); */
	
	drop ID PREDS_TAH JELUD_TAH FIBR_JELUD 
		 A_V_BLOK RAZRIV DRESSLER ZSN P_IM_STEN LET_IS
		 R_AB_3_n zab_leg_04; *dropping unnecissary outcome vars;
	
	rename R_AB_3_n2=R_AB_3_n;
	rename zab_leg_042=zab_leg_04;
/* 	rename REC_IM2=REC_IM; */
/* 	rename OTEK_LANC2=OTEK_LANC; */
/* 	rename FIBR_PREDS2=FIBR_PREDS; */
run;

proc contents data=mi_median order=varnum; *list of vars in data set in logical order;
proc contents data=mi_NN order=varnum;

/* proc means data=mi_cleaned2 median; */
/* 	var STENOK_AN; */
/* 	where (STENOK_AN <= 2) OR (STENOK_AN >= 3); */
/* run; */

/*STAGE 1b: Data Cleaning (DONE BY KEELY IN R!)*/
/*STAGE 2: Descriptive Stats and Visualization*/
/*STAGE 3: Modeling*/
	/*Model 1: Stepwise forward selection with 10-fold cross validation*/
proc glmselect data=mi_median plots=all seed=123; *best subset selection with 10-fold cross validation (max model size=8);
	class SEX--zab_leg_06(ref='0') O_L_POST--fibr_ter_08(ref='0') TIME_B_S--TRENT_S_n(ref='0') R_AB_3_n--zab_leg_04(ref='0');
	partition fraction(validate=0.1); *specifies that 10% of the data is used for test/validation;
	model REC_IM = Age--TRENT_S_n / selection=stepwise(choose=cv stop=sl slentry=0.15 slstay=0.15) 
	cvmethod=random(10);
run;
proc glmselect data=mi_median;   *Relapse of MI models with best Subset;
	class NA_R_3_n np_01 zab_leg_01;
	model REC_IM = NA_R_3_n|np_01|AGE|zab_leg_01 / 
		selection=stepwise(choose=bic stop=sl slentry=0.15 slstay=0.15); *Requests best 1 model(s) for each model size;
run;

	/*Model 2: Best Subset selection*/
/* proc logistic data=mi_median;   *Relapse of MI models with best Subset; */
/* 	class SEX--zab_leg_06(ref='0') O_L_POST--fibr_ter_08(ref='0') TIME_B_S--TRENT_S_n(ref='0') R_AB_3_n--zab_leg_04(ref='0'); */
/* 	model REC_IM = Age--TRENT_S_n /  */
/* 		selection=score best=1 stop=8; *Requests best 1 model(s) for each model size; */
/* run; */
/* proc glmselect data=mi_median;   *Relapse of MI models with best Subset; */
/* 	model REC_IM = STENOK_AN|np_01|zab_leg_02|K_SH_POST|GT_POST|L_BLOOD|R_AB_3_n|NA_R_2_n /  */
/* 		selection=stepwise(choose=bic stop=sl slentry=0.05 slstay=0.10); *Requests best 1 model(s) for each model size; */
/* run; */


	/*Model 3: LASSO With 10-fold cross validation*/
/* proc standard data=mi_median out=mi_median_standard; */
/* 	var AGE S_AD_ORIT D_AD_ORIT K_BLOOD Na_BLOOD ALT_BLOOD AST_BLOOD L_BLOOD ROE; */
/* run; */
/* proc glmselect data=mi_median_standard plots=all; *best subset selection with 10-fold cross validation (max model size=8); */
/* 	class SEX--zab_leg_06 O_L_POST--fibr_ter_08 TIME_B_S--NOT_NA_3_n; */
/* 	partition fraction(validate=0.1); *specifies that 10% of the data is used for test/validation; */
/* 	model REC_IM = Age--TRENT_S_n / selection=LASSO(choose=cv stop=none) cvmethod=random(10); */
/* run; *for some reason coefficients are not going to zero?; */




/*--------------------Nearest Neighbor Analyses-------------------*/

	/*Model 1: Stepwise forward selection with 10-fold cross validation*/
proc glmselect data=mi_NN plots=all seed=123; *best subset selection with 10-fold cross validation (max model size=8);
	class SEX--zab_leg_06(ref='0') O_L_POST--fibr_ter_08(ref='0') TIME_B_S--TRENT_S_n(ref='0') R_AB_3_n--zab_leg_04(ref='0');
	partition fraction(validate=0.1); *specifies that 10% of the data is used for test/validation;
	model REC_IM = Age--TRENT_S_n / selection=stepwise(choose=cv stop=sl slentry=0.15 slstay=0.15) 
	cvmethod=random(10);
run;
proc glmselect data=mi_NN;   *Relapse of MI models with best Subset;
	class NA_R_3_n np_01 zab_leg_01 GEPAR_S_n ritm_ecg_p_07;
	model REC_IM = NA_R_3_n|np_01|AGE|GEPAR_S_n|ritm_ecg_p_07 / 
		selection=stepwise(choose=bic stop=sl slentry=0.15 slstay=0.15); *Requests best 1 model(s) for each model size;
run;


/* 	Model 2: Best Subset selection */
/* proc logistic data=mi_NN;   *Relapse of MI models with best Subset; */
/* 	model REC_IM = Age--TRENT_S_n /  */
/* 		selection=score best=1 stop=8; *Requests best 1 model(s) for each model size; */
/* run; */
/* proc glmselect data=mi_median;   *Relapse of MI models with best Subset; */
/* 	model REC_IM = STENOK_AN np_01 zab_leg_02 K_SH_POST GT_POST L_BLOOD R_AB_3_n NA_R_2_n /  */
/* 		selection=score best=1 stop=8; *Requests best 1 model(s) for each model size; */
/* run; */
/* proc reg data=mi_median; */
/* 	model REC_IM = STENOK_AN np_01 zab_leg_02 K_SH_POST GT_POST L_BLOOD R_AB_3_n NA_R_2_n / */
/* 		vif collinoint; *requests collinearity analysis after adjusting out the intercept; */
/* run; */

/*List of categorical vars
SEX
INF_ANAM
STENOK_AN
FK_STENOK
IBS_POST
GB
SIM_GIPERT
ZSN_A
nr_11
nr_01
nr_02
nr_03
nr_04
nr_07
nr_08
np_01
np_04
np_05
np_07
np_08
np_09
np_10
endocr_01
endocr_02
endocr_03
zab_leg_01
zab_leg_02
zab_leg_03
zab_leg_06
O_L_POST
K_SH_POST
MP_TP_POST
SVT_POST
GT_POST
FIB_G_POST
ant_im
lat_im
inf_im
post_im
IM_PG_P
ritm_ecg_p_01
ritm_ecg_p_02
ritm_ecg_p_04
ritm_ecg_p_06
ritm_ecg_p_07
ritm_ecg_p_08
n_r_ecg_p_01
n_r_ecg_p_02
n_r_ecg_p_03
n_r_ecg_p_04
n_r_ecg_p_05
n_r_ecg_p_06
n_r_ecg_p_08
n_r_ecg_p_09
n_r_ecg_p_10
n_p_ecg_p_01
n_p_ecg_p_03
n_p_ecg_p_04
n_p_ecg_p_05
n_p_ecg_p_06
n_p_ecg_p_07
n_p_ecg_p_08
n_p_ecg_p_09
n_p_ecg_p_10
n_p_ecg_p_11
n_p_ecg_p_12
fibr_ter_01
fibr_ter_02
fibr_ter_03
fibr_ter_05
fibr_ter_06
fibr_ter_07
fibr_ter_08
TIME_B_S
NITR_S
LID_S_n
B_BLOK_S_n
ANT_CA_S_n
GEPAR_S_n
ASP_S_n
TIKL_S_n
TRENT_S_n
R_AB_1_n
NA_R_1_n
NOT_NA_1_n
R_AB_2_n
NA_R_2_n
NOT_NA_2_n
NA_R_3_n*/
