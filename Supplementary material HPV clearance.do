/* This do file is written in STATA 15, and contains data management and 
analyses of data from the H2M study, as used in the article by Wijstma et al.(2022) 
entitled "Approaches to estimating clearance rates for HHPV groupings: 
a systematic review and data examples".

** Chapter 1: Data management: preparing datasets including only incident or prevalent infections
		* 1.1 labelling raw data
		* 1.2 general data management
		* 1.3 preparing dataset including only incident infections
		* 1.4 preparing dataset including only prevalent infections

** Chapter 2: Data management: preparing datasets for the analysis of each specific approach
		* 2.1 dataset for approach A1
		* 2.2 dataset for approach A2
		* 2.3 dataset for approach A3
		* 2.4 dataset for approach B1
		* 2.5 dataset for approach B2 and B3
		* 2.6 dataset for approach C1
		* 2.7 dataset for approach C2 and C3

** Chapter 3: Estimating grouped clearance rates
		* 3.1 clearance rates according to approach A1, A2 and A3
		* 3.2 clearance rates according to approach B1, B2 and B3
		* 3.3 clearance rates according to approach C1, C2 and C3
		
*/

********************************************************************************
*** Chapter 1
********************************************************************************

********************* Chapter 1.1: labelling data

label variable persid "Participant ID"
label variable sampledat "Date of HPV sampling"
label variable round "Visit round (i.e., first, second, third, etc.)"
label variable Nrounds "the total number of visits the a participant has completed"

label variable Anal_6 "Anal HPV type 6 present (y/n)"
label variable Anal_11 "Anal HPV type 11 present (y/n)"
label variable Anal_16 "Anal HPV type 16 present (y/n)"
label variable Anal_18 "Anal HPV type 18 present (y/n)"
label variable Anal_31 "Anal HPV type 31 present (y/n)"
label variable Anal_33 "Anal HPV type 33 present (y/n)"
label variable Anal_34 "Anal HPV type 34 present (y/n)"
label variable Anal_35 "Anal HPV type 35 present (y/n)"
label variable Anal_39 "Anal HPV type 39 present (y/n)"
label variable Anal_40 "Anal HPV type 40 present (y/n)"
label variable Anal_42 "Anal HPV type 42 present (y/n)"
label variable Anal_43 "Anal HPV type 43 present (y/n)"
label variable Anal_44 "Anal HPV type 44 present (y/n)"
label variable Anal_45 "Anal HPV type 45 present (y/n)"
label variable Anal_51 "Anal HPV type 51 present (y/n)"
label variable Anal_52 "Anal HPV type 52 present (y/n)"
label variable Anal_53 "Anal HPV type 53 present (y/n)"
label variable Anal_54 "Anal HPV type 54 present (y/n)"
label variable Anal_56 "Anal HPV type 56 present (y/n)"
label variable Anal_58 "Anal HPV type 58 present (y/n)"
label variable Anal_59 "Anal HPV type 59 present (y/n)"
label variable Anal_66 "Anal HPV type 66 present (y/n)"
label variable Anal_6873 "Anal HPV type 68/73 present (y/n)"
label variable Anal_70 "Anal HPV type 70 present (y/n)"
label variable Anal_74 "Anal HPV type 74 present (y/n)"

******************** Chapter 1.2: General data management

* make time-variables: start- and stop dates for each visit interval
sort persid sampledat
by persid: gen tstart=sampledat[_n-1] 
label variable tstart "Starting date of observation period"
by persid: gen tstop=sampledat 
label variable tstop "Stopping date of observation period"
format tstart tstop %d

* make time-variables: duration of visit intervals; middate between visits
gen difvisdt=(tstop - tstart) 
label variable difvisdt "time in days between visit x and visit x-1"

gen middate=.
replace middate=tstart+0.5*difvisdt
label variable middate "Date between previous visit and midpoint to subsequent visit"
format middate %d

gen middays =(middate - tstart)
label variable middays "Days between previous visit and midpoint to subsequent visit"

* make time-variablesL time in days/months/years since first and subsequent visits
sort persid sampledat
by persid: gen t_y=(sampledat-sampledat[1]) / 365.25
label variable t_y "Time in years between first visit and subsequent visits"
by persid: gen t_m=(sampledat-sampledat[1]) / 30.5
label variable t_m "Time in months between first visit and subsequent visits"
by persid: gen t_d=(sampledat-sampledat[1])
label variable t_d "Time in days between first visit and subsequent visits"

* drop participants with <4 effective visits (see Methods)
drop if Nrounds ==1 | Nrounds ==2 | Nrounds ==3

* save dataset
save "H2M_all.dta", replace

******************* Chapter 1.3: Preparing dataset including only incident infections

/*Approaches A1, B1 and C1 include only incident infections. Thus, we create a 
dataset that includes only incident infections. Note that if participant X has an HPV
type (e.g., HPV16) at baseline, X can still contribute person-time to the analysis
later in time, when HPV16 is cleared and re-acquired.

In a loop, we will perform the following steps for each of the 25 HPV types:

Step 1: Clone the HPV type variable
Step 2: Create variable 'prevalent'. Change the value to 1 if an HPV type is
present at baseline
Step 3: If an infection is deemed prevalent: remove the corresponding positive 
value in the HPV type variable 
Step 4: Extrapolate this 'missing' value indicating the prevalent infection as
long as the infection remains prevalent (i.e. the previous (n-1) measure was 
deemed prevalent and the infection is also positive now.). Stop extrapolation 
if the prevalent infection clears, allowing for a new incident infection within 
the same person to contribute to CR estimation. */

*Note: the maximum number of follow-ups in the H2M study was 6

foreach hpv of varlist Anal_6 - Anal_74 {
  *step 1
  clonevar `hpv'_2 = `hpv' 
  
  *step 2
  generate prev_`hpv'  = 0 
  replace prev_`hpv'=1 if round == 1 & `hpv' ==1
  label variable prev_`hpv' "HPV type is prevalent at baseline" 
  *step3
  replace `hpv'=. if prev_`hpv' ==1 
  *step 4
  forvalues i = 2/5 {
  by persid: replace prev_`hpv'=1 if round== `i' & `hpv'_2==1 ///
			 & prev_`hpv'[_n-1]==1
  replace `hpv'=. if prev_`hpv' ==1
  }
  
  drop prev_`hpv' `hpv'_2
}

*save dataset
save "H2M_incident.dta", replace

******************* Chapter 1.4: Preparing dataset including only prevalent infections

/*Approaches A2, B2 and C2 include only prevalent infections. Thus, we create a 
dataset that includes only prevalent infections. 

Step 1: Clone the HPV type variable
Step 2: Create variable 'prev_AnalXX'. Change the value to 1 if an HPV type is
present at baseline
Step 3: Extrapolate values from 'prev_AnalXX' when someone:
	A) remains testing positive for the baseline-prevalent HPV type (step 3A)
	B) first tests negative for the baseline-prevalent HPV type (step 3B)		
Do not extrapolate values from 'prev_AnalXX'when somebody acquires an infection 
after clearing the prevalent infection (step 3C)
			
Step 4:  If an infection is NOT prevalent: remove the corresponding positive 
value in the HPV type variable 
*/

foreach hpv of varlist Anal_6 - Anal_74 {
 *step 1
 clonevar `hpv'_2 = `hpv' 
 *step 2
 generate prev_`hpv' =0 
 replace prev_`hpv'=1 if round == 1 & `hpv'==1
 
	*step 3A
	forvalues i = 2/5 {
	replace prev_`hpv'=1 if round== `i' & `hpv'==1 & `hpv'[_n-1]==1 & ///
	prev_`hpv'[_n-1]==1
	}	
	*step 3B
	forvalues i = 2/5 {
	replace prev_`hpv'=1 if round == `i' & `hpv'==0 & prev_`hpv'[_n-1]==1
	}	
	*step 3C
	forvalues i = 2/5 {
	replace prev_`hpv'=0 if round == `i' & `hpv'==1 & prev_`hpv'[_n-1]==0
	}
	
*step 4	
replace `hpv'=. if prev_`hpv'==0

drop prev_`hpv' `hpv'_2
}

* save dataset
save "H2M_prevalent.dta", replace

******************* Chapter 2.1 : preparing datasets for A1 ********************

*****************************************
/*The A1 approach measures clearance and PMO on the HPV type-level. 
Only incident infections are included.
The source dataset is 
"Dta\H2M_incident.dta"*/
*****************************************

use "H2M_incident.dta", clear

*for each HPV type, we should censor infections incident at the final visit, 
*since these cannot possibly be cleared.

foreach hpv of varlist Anal_6 - Anal_74 {
  replace `hpv'=. if round == Nrounds & `hpv'==1 & `hpv'[_n-1]==0
  }
  
*We define clearance events per HPV type (as a positive test followed by a 
*negative test, within an individual)
 foreach hpv of varlist Anal_6 - Anal_74 {
 by persid: gen event_`hpv'=1 if `hpv'==0 & `hpv'[_n-1]==1
}

*We define person-time. Person-time starts at the midpoint before the first positive
*visits. Person-time ends at the midpoint betwen the last positive and first 
*negative visit. We previously defined the 'middays' variable, as half the
*number of days between visits.

*preparation
*the data was previously structured broad (i.e., one row per person and timepoint)
*we will reshape, so that there is one row per person, timepoint AND HPV type.
*This way, we can count person-time at risk for clearance separately for each HPV type
*(as is defined in approaches of subset A)
gen uniquenumber =.
sort persid round
replace uniquenumber = _n
rename (Anal_6-Anal_74) Anal=
reshape long Anal, i(uniquenumber) j(HPVtype) string
rename Anal present

sort persid HPVtype round
tostring HPVtype, gen(HPVtype2)
tostring persid, gen(persid2)
gen persidtype = persid2 + HPVtype
label variable persidtype "HPV type X within participant number X"

	*part 1: first positive
	gen firstPos=.
	label variable firstPos "first positive test result"
	sort persidtype round firstPos
	by persidtype: replace firstPos = round if present==1 & present[_n-1]==0

	*part 2: first negative
	gen firstNeg=.
	label variable firstNeg "first negative test result after a positive result"
	sort persidtype round firstNeg
	by persidtype: replace firstNeg = round if present ==0 & present[_n-1]==1
	
	*part 3: person-days
	gen PDO=.
	label variable PDO "person-days of observation, at risk for clearance"
	replace PDO= middays if round == firstPos
	by persidtype: replace PDO = difvisdt if present==1 & present[_n-1]==1 
	replace PDO = middays if round == firstNeg

	*Part 4: convert to PMO
	gen PMO=.
	replace PMO= PDO/30.5 


*Now, we can make grouped event variables, and grouped person-time variabeles
*since this piece of code will be repeated for approach A1, A2 and A3, 
*we will write the following code into a program

program makegroupings 

		***** Grouped event-variables, per grouping
		*HPV grouping: Any HPV (including all 25 types identified by the LiPa array)
		gen events_Any=.
		replace events_Any = 1 if event_Anal_6 ==1
		replace events_Any = 2 if event_Anal_11 ==1
		replace events_Any = 3 if event_Anal_16 ==1
		replace events_Any = 4 if event_Anal_18 ==1
		replace events_Any = 5 if event_Anal_31 ==1
		replace events_Any = 6 if event_Anal_33 ==1
		replace events_Any = 7 if event_Anal_35 ==1
		replace events_Any = 8 if event_Anal_39 ==1
		replace events_Any = 9 if event_Anal_45 ==1
		replace events_Any = 10 if event_Anal_51 ==1
		replace events_Any = 11 if event_Anal_52 ==1
		replace events_Any = 12 if event_Anal_56 ==1
		replace events_Any = 13 if event_Anal_58 ==1
		replace events_Any = 14 if event_Anal_59 ==1
		replace events_Any = 15 if event_Anal_34 ==1
		replace events_Any = 16 if event_Anal_40 ==1
		replace events_Any = 17 if event_Anal_42 ==1
		replace events_Any = 18 if event_Anal_43 ==1
		replace events_Any = 19 if event_Anal_44 ==1
		replace events_Any = 20 if event_Anal_53 ==1
		replace events_Any = 21 if event_Anal_54 ==1
		replace events_Any = 22 if event_Anal_66 ==1
		replace events_Any = 23 if event_Anal_6873 ==1
		replace events_Any = 24 if event_Anal_70 ==1
		replace events_Any = 25 if event_Anal_74 ==1

		****High risk (HR) HPV
		*the 12 HR types are: 16, 18, 31, 33, 35, 39, 45, 51, 52, 56, 58, 59
		gen events_HR=.
		replace events_HR = 1 if  event_Anal_16 ==1
		replace events_HR = 2 if  event_Anal_18 ==1
		replace events_HR = 3 if  event_Anal_31 ==1
		replace events_HR = 4 if  event_Anal_33 ==1
		replace events_HR = 5 if  event_Anal_35 ==1
		replace events_HR = 6 if  event_Anal_39 ==1
		replace events_HR = 7 if  event_Anal_45 ==1
		replace events_HR = 8 if  event_Anal_51 ==1
		replace events_HR = 9 if  event_Anal_52 ==1
		replace events_HR = 10 if  event_Anal_56 ==1
		replace events_HR = 11 if  event_Anal_58 ==1
		replace events_HR = 12 if  event_Anal_59 ==1

		**Low risk (LR) HPV
		*the 13 LR types are: 6,11,34,40,42,43,44,53,54,66,68/73,70,74
		gen events_LR=.
		replace events_LR = 1 if  event_Anal_6 ==1
		replace events_LR = 2 if  event_Anal_11 ==1
		replace events_LR = 3 if  event_Anal_34 ==1
		replace events_LR = 4 if  event_Anal_40 ==1
		replace events_LR = 5 if  event_Anal_42 ==1
		replace events_LR = 6 if  event_Anal_43 ==1
		replace events_LR = 7 if  event_Anal_44 ==1
		replace events_LR = 8 if  event_Anal_53 ==1
		replace events_LR = 9 if  event_Anal_54 ==1
		replace events_LR = 10 if  event_Anal_66 ==1
		replace events_LR = 11 if  event_Anal_6873 ==1
		replace events_LR = 12 if  event_Anal_70 ==1
		replace events_LR = 13 if  event_Anal_74 ==1

		**Bivalent (2v) HPV: types 16/18
		gen events_2v=.
		replace events_2v = 1 if  event_Anal_16 ==1
		replace events_2v = 2 if  event_Anal_18 ==1

		**Quadrivalent (4v) HPV: types 6/11/16/18
		gen events_4v=.
		replace events_4v = 1 if  event_Anal_6 ==1
		replace events_4v = 2 if  event_Anal_11 ==1
		replace events_4v = 3 if  event_Anal_16 ==1
		replace events_4v = 4 if  event_Anal_18 ==1

		**Nonavalent (9v) HPV
		*the 9v vaccine types include: 6, 11, 16, 18, 31, 33, 45, 52 and 58 
		gen events_9v=.
		replace events_9v = 1 if  event_Anal_6 ==1
		replace events_9v = 2 if  event_Anal_11 ==1
		replace events_9v = 3 if  event_Anal_16 ==1
		replace events_9v = 4 if  event_Anal_18 ==1
		replace events_9v = 5 if  event_Anal_31 ==1
		replace events_9v = 6 if  event_Anal_33 ==1
		replace events_9v = 7 if  event_Anal_45 ==1
		replace events_9v = 8 if  event_Anal_52 ==1
		replace events_9v = 9 if  event_Anal_58 ==1

		***** Grouped person-time variables, per grouping
		* Grouping: Any HPV (including all 25 hpv types identified by the LiPa-25)
		gen PMO_Any=.
		replace PMO_Any = PMO if HPVtype == "Anal_6" 
		replace PMO_Any = PMO if HPVtype == "Anal_11" 
		replace PMO_Any = PMO if HPVtype == "Anal_16" 
		replace PMO_Any = PMO if HPVtype == "Anal_18" 
		replace PMO_Any = PMO if HPVtype == "Anal_31" 
		replace PMO_Any = PMO if HPVtype == "Anal_33"
		replace PMO_Any = PMO if HPVtype == "Anal_35"
		replace PMO_Any = PMO if HPVtype == "Anal_39"
		replace PMO_Any = PMO if HPVtype == "Anal_45"
		replace PMO_Any = PMO if HPVtype == "Anal_51"
		replace PMO_Any = PMO if HPVtype == "Anal_52"
		replace PMO_Any = PMO if HPVtype == "Anal_56"
		replace PMO_Any = PMO if HPVtype == "Anal_58"
		replace PMO_Any = PMO if HPVtype == "Anal_59"
		replace PMO_Any = PMO if HPVtype == "Anal_34"
		replace PMO_Any = PMO if HPVtype == "Anal_40"
		replace PMO_Any = PMO if HPVtype == "Anal_42"
		replace PMO_Any = PMO if HPVtype == "Anal_43"
		replace PMO_Any = PMO if HPVtype == "Anal_44"
		replace PMO_Any = PMO if HPVtype == "Anal_53"
		replace PMO_Any = PMO if HPVtype == "Anal_54"
		replace PMO_Any = PMO if HPVtype == "Anal_66"
		replace PMO_Any = PMO if HPVtype == "Anal_6873"
		replace PMO_Any = PMO if HPVtype == "Anal_70"
		replace PMO_Any = PMO if HPVtype == "Anal_74"

		*HR HPV
		*the 12 HR types are: 16, 18, 31, 33, 35, 39, 45, 51, 52, 56, 58, 59
		gen PMO_HR=.
		replace PMO_HR = PMO if HPVtype == "Anal_16" 
		replace PMO_HR = PMO if HPVtype == "Anal_18" 
		replace PMO_HR = PMO if HPVtype == "Anal_33"
		replace PMO_HR = PMO if HPVtype == "Anal_35"
		replace PMO_HR = PMO if HPVtype == "Anal_39"
		replace PMO_HR = PMO if HPVtype == "Anal_45"
		replace PMO_HR = PMO if HPVtype == "Anal_51"
		replace PMO_HR = PMO if HPVtype == "Anal_52"
		replace PMO_HR = PMO if HPVtype == "Anal_56"
		replace PMO_HR = PMO if HPVtype == "Anal_58"
		replace PMO_HR = PMO if HPVtype == "Anal_59"

		*LR HPV
		*the 13 LR types are: 6,11,34,40,42,43,44,53,54,66,68/73,70,74
		gen PMO_LR=.
		replace PMO_LR = PMO if HPVtype == "Anal_6" 
		replace PMO_LR = PMO if HPVtype == "Anal_11" 
		replace PMO_LR = PMO if HPVtype == "Anal_34"
		replace PMO_LR = PMO if HPVtype == "Anal_40"
		replace PMO_LR = PMO if HPVtype == "Anal_42"
		replace PMO_LR = PMO if HPVtype == "Anal_43"
		replace PMO_LR = PMO if HPVtype == "Anal_44"
		replace PMO_LR = PMO if HPVtype == "Anal_53"
		replace PMO_LR = PMO if HPVtype == "Anal_54"
		replace PMO_LR = PMO if HPVtype == "Anal_66"
		replace PMO_LR = PMO if HPVtype == "Anal_6873"
		replace PMO_LR = PMO if HPVtype == "Anal_70"
		replace PMO_LR = PMO if HPVtype == "Anal_74"

		*2v HPV
		*2v vaccine types include 16 and 18
		gen PMO_2v=.
		replace PMO_2v = PMO if HPVtype == "Anal_16" 
		replace PMO_2v = PMO if HPVtype == "Anal_18" 


		*4v HPV
		*the 4v vaccine types include 6, 11, 16 & 18
		gen PMO_4v=.
		replace PMO_4v = PMO if HPVtype == "Anal_6" 
		replace PMO_4v = PMO if HPVtype == "Anal_11" 
		replace PMO_4v = PMO if HPVtype == "Anal_16"
		replace PMO_4v = PMO if HPVtype == "Anal_18"


		*9v HPV
		*the 9v vaccine types include: 6, 11, 16, 18, 31, 33, 45, 52 and 58 
		gen PMO_9v=.
		replace PMO_9v = PMO if HPVtype == "Anal_6" 
		replace PMO_9v = PMO if HPVtype == "Anal_11" 
		replace PMO_9v = PMO if HPVtype == "Anal_16"
		replace PMO_9v = PMO if HPVtype == "Anal_18"
		replace PMO_9v = PMO if HPVtype == "Anal_31"
		replace PMO_9v = PMO if HPVtype == "Anal_33"
		replace PMO_9v = PMO if HPVtype == "Anal_45"
		replace PMO_9v = PMO if HPVtype == "Anal_52"
		replace PMO_9v = PMO if HPVtype == "Anal_58"

		end
		
save "H2m_ApproachA1.dta", replace

******************* Chapter 2.2 : preparing datasets for A2 ********************

*****************************************
/*The A2 approach measures clearance and PMO on the HPV type-level. 
Only prevalent infections are included.
The source dataset is "Dta\H2M_prevalent.dta"*/
*****************************************

use "H2M_prevalent.dta", clear 

*We define clearance events per HPV type (as a positive test followed by a 
*negative test, within an individual)
 foreach hpv of varlist Anal_6 - Anal_74 {
 by persid: gen event_`hpv'=1 if `hpv'==0 & `hpv'[_n-1]==1
}

*We define person-time. Person-time starts at the first visit (since only prevalent
*infections are included). Person-time ends at the midpoint between the last 
*positive and first negative visit. 

*preparations
*the data was previously structured broad (i.e., one row per person and timepoint)
*we will reshape, so that there is one row per person, timepoint AND HPV type.
*This way, we can count person-time at risk for clearance separately for each HPV type
*(as is defined in approaches of subset A)
gen uniquenumber =.
sort persid round
replace uniquenumber = _n
rename (Anal_6-Anal_74) Anal=
reshape long Anal, i(uniquenumber) j(HPVtype) string
rename Anal present

sort persid HPVtype round
tostring HPVtype, gen(HPVtype2)
tostring persid, gen(persid2)
gen persidtype = persid2 + HPVtype
label variable persidtype "HPV type X within participant number X"

	*part 1: first positive (at visit round 1!)
	gen firstPos=.
	sort persidtype round firstPos
	by persidtype: replace firstPos=1 if round==1 & present==1
	
	*part 2 : first negative
	gen firstNeg=.
	sort persidtype round firstNeg
	by persidtype:  replace firstNeg= round if present==0 & present[_n-1]==1

				
	*part 3: person-days. person time consists of:
	* (i) the full difvisdt between two positive visits, and
	* (ii) the middays before the first negative visit 
	
	sort persidtype round difvisdt
	gen PDO=.
	by persidtype: replace PDO= difvisdt if present==1 & present[_n-1]==1 //(i)
	replace PDO= middays if round == firstNeg							  //(ii)

	*part 4: convert to PMO
	gen PMO=.
	replace PMO= PDO/30.5 

*Now, similar to chapter 2.1, we make grouped event variables and grouped PMO variables
*we can use the program we created, named 'makegroupings'
makegroupings 

save "H2m_ApproachA2.dta", replace

******************* Chapter 2.3 : preparing datasets for A3 ********************

*****************************************
/*The A3 approach measures clearance and PMO on the HPV type-level. 
Both incident and baseline-prevalent infections are included.
The source dataset is "Dta\H2M_all.dta"*/
*****************************************

use "H2M_all.dta", clear 

*censor observations that become positive at the last visit							
foreach hpv of varlist Anal_6 - Anal_74 {
	replace `hpv'=. if round == Nrounds & `hpv'==1 & `hpv'[_n-1]==0	
  }  
  
 *We define clearance events per HPV type (as a positive test followed by a 
*negative test, within an individual)
 foreach hpv of varlist Anal_6 - Anal_74 {
 by persid: gen event_`hpv'=1 if `hpv'==0 & `hpv'[_n-1]==1
}

*We define person-time. Person-time starts at the first visit (since only prevalent
*infections are included). Person-time ends at the midpoint between the last 
*positive and first negative visit. 

*preparations
*the data was previously structured broad (i.e., one row per person and timepoint)
*we will reshape, so that there is one row per person, timepoint AND HPV type.
*This way, we can count person-time at risk for clearance separately for each HPV type
*(as is defined in approaches of subset A)
gen uniquenumber =.
sort persid round
replace uniquenumber = _n
rename (Anal_6-Anal_74) Anal=
reshape long Anal, i(uniquenumber) j(HPVtype) string
rename Anal present

sort persid HPVtype round
tostring HPVtype, gen(HPVtype2)
tostring persid, gen(persid2)
gen persidtype = persid2 + HPVtype
label variable persidtype "HPV type X within participant number X"

	*part 1: first positive (at visit round 1!)
	gen firstPos=.
	sort persidtype round firstPos
	by persidtype: replace firstPos=1 if round==1 & present==1
	
	*part 2 : first negative
	gen firstNeg=.
	sort persidtype round firstNeg
	by persidtype:  replace firstNeg= round if present==0 & present[_n-1]==1
				
	*part 3: person-days. person time consists of:
	/*PMO consists of:
	- the middays before the first positive (in case of incident infections)
	- the difvisdt betweemm two positive visits (for incident or prevalent infection)
	- the middays before the first negative (for incident of prevalent infection)
	*/

	gen PDO=.
	sort persidtype round
	replace PDO= middays if round == firstPos
	by persidtype: replace PDO= difvisdt if present==1 & present[_n-1]==1 
	replace PDO= middays if round == firstNeg
		
	*part 4: convert to PMO
	gen PMO=.
	replace PMO= PDO/30.5 

	
*Now, similar to chapter 2.1, we make grouped event variables and grouped PMO variables
*we can use the program we created, named 'makegroupings'
makegroupings 

save "H2m_ApproachA3.dta", replace

******************* Chapter 2.4 : preparing datasets for B1 ********************

*****************************************
/*The B1 approach counts a maximum of one clearance event. When the
first specific HPV type has cleared, no more person-time is counted.
Only incident infections are included.
The source dataset is "Dta\H2M_incident.dta"*/
*****************************************

/*General strategy to code approach subset B:
The first few steps are identical to the steps in approach A:

using a foreach variable loop, the following variables will be created
for each HPV type:
- first positive test
- first negative test
- clearance event 
- Person-time

Therafter, we will calculate the clearance events and person-time per grouping.
This is where approach set B differs from approach set A:
We only want to save the results of the **first** clearance, not of any
possible second clearances. We also want to stop the counting of person-time 
after the first clearance

Thus, instead of summing up all the type-specific clearance events from every row
(like was done in subset A, which had no maximum of clearance events), we will 
let the grouped clearance variable equal "1", if a type-specific clearance event
occurred for one or more types (i.e. in one or more rows)

Thus, the grouped events variable would be indexed as (for example):

replace events_2v = 1 if event_Anal_16 == 1 | event_Anal_18 ==1 

We will censor the person-time that comes after the first clearance event,
and we will replace the grouped_PMO variable with the PMO_variable for the rows
between infection and clearance.

*/

/*we should censor infections that are incident at the final visit, since these
cannot be cleared */
foreach hpv of varlist Anal_6 - Anal_74 {
  replace `hpv'=. if round == Nrounds & `hpv'==1 & `hpv'[_n-1]==0
  }

 program makeappB 
 
*indicate the first positive test result of each HPV type
foreach hpv of varlist Anal_6 - Anal_74 {
	gen firstPos`hpv'=.
	by persid: replace firstPos`hpv'= round if `hpv'==1 & `hpv'[_n-1]==0								
  }
  
*indicate the first negative test result of each HPV type
foreach hpv of varlist Anal_6 - Anal_74 {
  gen firstNeg`hpv'=.
  by persid: replace firstNeg`hpv' = round if `hpv' ==0 & `hpv'[_n-1]==1
  }
  
 *Now, we move to the groupings:
 
********************* approach B1
*bivalent (2v) grouping: including types 16 & 18
*********************

********* Clearance event

**Step 1: first positive of the grouping

	*put the first positives of all separate 2v-types in a single variable
	gen firstPos_2v =.
	replace firstPos_2v = 1 if firstPosAnal_16 !=. | firstPosAnal_18 !=.

	 *keep only keep the FIRST first positive (note: max 5 visits in H2M dataset)
	by persid: replace firstPos_2v =. if firstPos_2v[_n-1]==1 | ///
	firstPos_2v[_n-2]==1 | firstPos_2v[_n-3]==1 | firstPos_2v[_n-4]==1 
			   
	*indicate at which round (i.e., visit number) the person was first positive 
	*for the grouping. Carry this value forward (so that later, this can be compared
	*to the round number of the first negative value) 
	replace firstPos_2v =round if firstPos_2v ==1 
	sort persid round firstPos_2v
	by persid: carryforward firstPos_2v if missing(firstPos_2v), replace


**Step 2: first negative of the grouping

	gen firstNeg_2v =.
	replace firstNeg_2v = 1 if firstNegAnal_16 !=. | firstNegAnal_18 !=.
	*keep only the FIRST first negative
	by persid: replace firstNeg_2v =. if firstNeg_2v[_n-1]==1 | ///
	firstNeg_2v[_n-2]==1 | firstNeg_2v[_n-3]==1
									  
	*indicate at which round (i.e., visit number) the person first turned negative
	replace firstNeg_2v =round if firstNeg_2v ==1	

** Step 3: grouped event variable
*Because we carried forward the values for 'first positive for the grouping' 
*(within an individual), we can index the timing of a clearance event as:
*the value for firstPos is not empty, and the value for firstNeg is not empty at time X.

gen event_2v=.
replace event_2v=1	if firstNeg_2v !=. & firstPos_2v !=.

*********** Person-time  
/*PMO consists of:
	(i) 	the 'middays' (half the number of days between visits) 
			before the first positive visit
	(ii) 	the total number of days ('difvisdt') between two positive visits	
	(iiI) 	the 'middays' before the first negative visit		*/
	
	gen PDO_2v=.
	*(i) middays before first positive
	replace PDO_2v = middays if round == firstPos_2v

	*(ii) difvisdt between two positive visits clearance:
	sort persid round firstNeg_2v
	by persid: carryforward firstNeg_2v if missing(firstNeg_2v), replace	
	by persid: replace PDO_2v = difvisdt if round < firstNeg_2v	& round > firstPos_2v	

	*(iii) middays before the first negative
	replace PDO_2v = middays if round == firstNeg_2v
	 
	*convert to months (PMO)
	gen PMO_2v=.
	replace PMO_2v = PDO_2v/30.5

*We will now repeat the same code structure, but for the groupings 4v, 9v, LR, HR, and Any.


*********** approach B1
*quadrivalent (4v) grouping, including types 6,11,16,18
***********	
	
**** Clearance

	**** Step 1: first positive of the grouping
	gen firstPos_4v =.
	replace firstPos_4v = 1 if firstPosAnal_6 !=. | firstPosAnal_11 !=. | ///
								firstPosAnal_16 !=. | firstPosAnal_18 !=. 
			
	 *keep only keep the FIRST first positive (note: max 5 visits in H2M dataset)
	by persid: replace firstPos_4v =. if firstPos_4v[_n-1]==1 | ///
			  firstPos_4v[_n-2]==1 | firstPos_4v[_n-3]==1 | firstPos_4v[_n-4]==1 
			   
	*let the first positive indicate the roundnumber, and be carried forward
	replace firstPos_4v =round if firstPos_4v ==1 
	sort persid round firstPos_4v
	by persid: carryforward firstPos_4v if missing(firstPos_4v), replace

	***** Step 2: first negative of the grouping
	gen firstNeg_4v =.
	replace firstNeg_4v = 1 if firstNegAnal_6 !=. | firstNegAnal_11 !=. | ///
							   firstNegAnal_16 !=. | firstNegAnal_18 !=. 
													   
	*keep only the FIRST first negative
	by persid: replace firstNeg_4v =. if firstNeg_4v[_n-1]==1 | ///
									  firstNeg_4v[_n-2]==1 | firstNeg_4v[_n-3]==1
									  
	*indicate roundnumber, don't carry forward 
	replace firstNeg_4v =round if firstNeg_4v ==1	

	***** Step 3: grouped event variable
	gen event_4v=.
	replace event_4v=1	if firstNeg_4v !=. & firstPos_4v !=.

**** PMO 

/*for incident infections requiring 1 neg test for clearance, PMO consists of:
	- the middays before the first positive visit
	- the middays before the first negative visit
	- the difvisdt between two positive visits	*/
	
	gen PDO_4v=.

	*middays before first positive
	replace PDO_4v = middays if round == firstPos_4v

	*difvisdt between two positive visits clearance:
	*i.e. the 1st pos has already occurred, but the 1st neg has not yet occurred
	sort persid round firstNeg_4v
	by persid: carryforward firstNeg_4v if missing(firstNeg_4v), replace	
	by persid: replace PDO_4v = difvisdt if round < firstNeg_4v ///
										   & round > firstPos_4v	
	*middays before the first negative
	replace PDO_4v = middays if round == firstNeg_4v

	*convert to PMO
	gen PMO_4v =.
	replace PMO_4v = PDO_4v / 30.5


*********** approach B1
*Nonavalent (9V) vaccine grouping, including types 6,11,16,18,31, 33, 45, 52 and 58 
***********	
	
**** Clearance

	**** Step 1: first positive of the grouping
	gen firstPos_9v =.
	replace firstPos_9v = 1 if firstPosAnal_6 !=. | firstPosAnal_11 !=. | ///
								firstPosAnal_16 !=. | firstPosAnal_18 !=. | ///
								firstPosAnal_31 !=. | firstPosAnal_33 !=. | ///
								firstPosAnal_45 !=. | firstPosAnal_52 !=. | ///
								firstPosAnal_58 !=. 
			
	 *keep only keep the FIRST first positive (note: max 5 visits in H2M dataset)
	by persid: replace firstPos_9v =. if firstPos_9v[_n-1]==1 | ///
			  firstPos_9v[_n-2]==1 | firstPos_9v[_n-3]==1 | firstPos_9v[_n-4]==1 
			   
	*let the first positive indicate the roundnumber, and be carried forward
	replace firstPos_9v =round if firstPos_9v ==1 
	sort persid round firstPos_9v
	by persid: carryforward firstPos_9v if missing(firstPos_9v), replace

	***** Step 2: first negative of the grouping
	gen firstNeg_9v =.
	replace firstNeg_9v = 1 if firstNegAnal_6 !=. | firstNegAnal_11 !=. | ///
							   firstNegAnal_16 !=. | firstNegAnal_18 !=. | ///
							   firstNegAnal_31 !=. | firstNegAnal_33 !=. | ///
							   firstNegAnal_45 !=. | firstNegAnal_52 !=. | ///
							   firstNegAnal_58 !=. 
	*keep only the FIRST first negative
	by persid: replace firstNeg_9v =. if firstNeg_9v[_n-1]==1 | ///
									  firstNeg_9v[_n-2]==1 | firstNeg_9v[_n-3]==1
									  
	*indicate roundnumber, don't carry forward 
	replace firstNeg_9v =round if firstNeg_9v ==1	

	***** Step 3: grouped event variable
	gen event_9v=.
	replace event_9v=1	if firstNeg_9v !=. & firstPos_9v !=.

***** PMO 

/*for incident infections requiring 1 neg test for clearance, PMO consists of:
	- the middays before the first positive visit
	- the middays before the first negative visit
	- the difvisdt between two positive visits	*/
	
	gen PDO_9v=.

	*middays before first positive
	replace PDO_9v = middays if round == firstPos_9v

	*difvisdt between two positive visits clearance:
	*i.e. the 1st pos has already occurred, but the 1st neg has not yet occurred
	sort persid round firstNeg_9v
	by persid: carryforward firstNeg_9v if missing(firstNeg_9v), replace	
	by persid: replace PDO_9v = difvisdt if round < firstNeg_9v ///
										   & round > firstPos_9v	
	*middays before the first negative
	replace PDO_9v = middays if round == firstNeg_9v

	*convert to PMO
	gen PMO_9v =.
	replace PMO_9v = PDO_9v / 30.5
	

*********** approach B1
* Any HPV grouping, including all types identified by LiPa-25
***********	
	
****Clearance

	**** Step 1: first positive of the grouping
	gen firstPos_Any =.
	replace firstPos_Any = 1 if firstPosAnal_6 !=. | firstPosAnal_11 !=. | ///
			firstPosAnal_16 !=. | firstPosAnal_18 !=. | firstPosAnal_31 !=. | ///
			firstPosAnal_33 !=. | firstPosAnal_34 !=. | firstPosAnal_35 !=. | ///
			firstPosAnal_39 !=. | firstPosAnal_40 !=. | firstPosAnal_42 !=. | ///
			firstPosAnal_43 !=. | firstPosAnal_44 !=. | firstPosAnal_45 !=. | ///
			firstPosAnal_51 !=. | firstPosAnal_52 !=. | firstPosAnal_53 !=. | ///
			firstPosAnal_54 !=. | firstPosAnal_56 !=. | firstPosAnal_58 !=. | ///
			firstPosAnal_59 !=. | firstPosAnal_66 !=. | firstPosAnal_6873 !=. | ///
			firstPosAnal_70 !=. | firstPosAnal_74 !=. 

	*keep only keep the FIRST first positive (note: max 5 visits in H2M dataset)
	by persid: replace firstPos_Any =. if firstPos_Any[_n-1]==1 | ///
			  firstPos_Any[_n-2]==1 | firstPos_Any[_n-3]==1 | firstPos_Any[_n-4]==1 

	*let the first positive indicate the roundnumber, and be carried forward
	replace firstPos_Any =round if firstPos_Any ==1 
	sort persid round firstPos_Any
	by persid: carryforward firstPos_Any if missing(firstPos_Any), replace

	***** Step 2: first negative of the grouping
	gen firstNeg_Any =.
	replace firstNeg_Any = 1 if firstNegAnal_6 !=. | firstNegAnal_11 !=. | ///
			firstNegAnal_16 !=. | firstNegAnal_18 !=. | firstNegAnal_31 !=. | ///
			firstNegAnal_33 !=. | firstNegAnal_34 !=. | firstNegAnal_35 !=. | ///
			firstNegAnal_39 !=. | firstNegAnal_40 !=. | firstNegAnal_42 !=. | ///
			firstNegAnal_43 !=. | firstNegAnal_44 !=. | firstNegAnal_45 !=. | ///
			firstNegAnal_51 !=. | firstNegAnal_52 !=. | firstNegAnal_53 !=. | ///
			firstNegAnal_54 !=. | firstNegAnal_56 !=. | firstNegAnal_58 !=. | ///
			firstNegAnal_59 !=. | firstNegAnal_66 !=. | firstNegAnal_6873 !=. | ///
			firstNegAnal_70 !=. | firstNegAnal_74 !=. 
	*keep only the FIRST first negative
	by persid: replace firstNeg_Any =. if firstNeg_Any[_n-1]==1 | ///
									  firstNeg_Any[_n-2]==1 | firstNeg_Any[_n-3]==1
									  
	*indicate roundnumber, don't carry forward 
	replace firstNeg_Any =round if firstNeg_Any ==1	

	***** Step 3: grouped event variable
	gen event_Any=.
	replace event_Any=1	if firstNeg_Any !=. & firstPos_Any !=. 

**** PMO
/*for incident infections requiring 1 neg test for clearance, PMO consists of:
	- the middays before the first positive visit
	- the middays before the first negative visit
	- the difvisdt between two positive visits	*/
	
	gen PDO_Any=.

	*middays before first positive
	replace PDO_Any = middays if round == firstPos_Any

	*difvisdt between two positive visits clearance:
	*i.e. the 1st pos has already occurred, but the 1st neg has not yet occurred
	sort persid round firstNeg_Any
	by persid: carryforward firstNeg_Any if missing(firstNeg_Any), replace	
	by persid: replace PDO_Any = difvisdt if round < firstNeg_Any ///
										   & round > firstPos_Any	
	*middays before the first negative
	replace PDO_Any = middays if round == firstNeg_Any

	*convert to PMO
	gen PMO_Any =.
	replace PMO_Any = PDO_Any / 30.5


*********** approach B1
*High risk (HR) HPV GROUPING
*the 12 HR types are: 16, 18, 31, 33, 35, 39, 45, 51, 52, 56, 58, 59
***********	
	
****Clearance

	**** Step 1: first positive of the grouping
	gen firstPos_HR =.
	replace firstPos_HR = 1 if firstPosAnal_16 !=. | firstPosAnal_18 !=. | ///
								firstPosAnal_31 !=. | firstPosAnal_33 !=. | ///
								firstPosAnal_35 !=. | firstPosAnal_39 !=. | ///
								firstPosAnal_45 !=. | firstPosAnal_51 !=. | ///
								firstPosAnal_52 !=. | firstPosAnal_56 !=. | ///
								firstPosAnal_58 !=. | firstPosAnal_59 !=. 
			
	 *keep only keep the FIRST first positive (note: max 5 visits in H2M dataset)
	by persid: replace firstPos_HR =. if firstPos_HR[_n-1]==1 | ///
			  firstPos_HR[_n-2]==1 | firstPos_HR[_n-3]==1 | firstPos_HR[_n-4]==1 
			   
	*let the first positive indicate the roundnumber, and be carried forward
	replace firstPos_HR =round if firstPos_HR ==1 
	sort persid round firstPos_HR
	by persid: carryforward firstPos_HR if missing(firstPos_HR), replace

	***** Step 2: first negative of the grouping
	gen firstNeg_HR =.
	replace firstNeg_HR = 1 if firstNegAnal_16 !=. | firstNegAnal_18 !=. | ///
							   firstNegAnal_31 !=. | firstNegAnal_33 !=. | ///
							   firstNegAnal_35 !=. | firstNegAnal_39 !=. | ///
							   firstNegAnal_45 !=. | firstNegAnal_51 !=. | ///
							   firstNegAnal_52 !=. | firstNegAnal_56 !=. | ///
							   firstNegAnal_58 !=. | firstNegAnal_59 !=. 
							   
	*keep only the FIRST first negative
	by persid: replace firstNeg_HR =. if firstNeg_HR[_n-1]==1 | ///
									  firstNeg_HR[_n-2]==1 | firstNeg_HR[_n-3]==1
									  
	*indicate roundnumber, don't carry forward 
	replace firstNeg_HR =round if firstNeg_HR ==1	

	***** Step 3: grouped event variable
	gen event_HR=.
	replace event_HR=1	if firstNeg_HR !=. & firstPos_HR !=.

***** PMO 

/*for incident infections requiring 1 neg test for clearance, PMO consists of:
	- the middays before the first positive visit
	- the middays before the first negative visit
	- the difvisdt between two positive visits	*/
	
	gen PDO_HR=.

	*middays before first positive
	replace PDO_HR = middays if round == firstPos_HR

	*difvisdt between two positive visits clearance:
	*i.e. the 1st pos has already occurred, but the 1st neg has not yet occurred
	sort persid round firstNeg_HR
	by persid: carryforward firstNeg_HR if missing(firstNeg_HR), replace	
	by persid: replace PDO_HR = difvisdt if round < firstNeg_HR ///
										   & round > firstPos_HR	
	*middays before the first negative
	replace PDO_HR = middays if round == firstNeg_HR

	*convert to PMO
	gen PMO_HR =.
	replace PMO_HR = PDO_HR / 30.5


*********** approach B1
*Low risk (LR) HPV GROUPING
*the 13 LR types are: 6,11,34,40,42,43,44,53,54,66,68/73,70,74
***********	
	
*****Clearance

	**** Step 1: first positive of the grouping
	gen firstPos_LR =.
	replace firstPos_LR = 1 if firstPosAnal_6 !=. | firstPosAnal_11 !=. | ///
								firstPosAnal_34 !=. | firstPosAnal_40 !=. | ///
								firstPosAnal_42 !=. | firstPosAnal_43 !=. | ///
								firstPosAnal_44 !=. | firstPosAnal_53 !=. | ///
								firstPosAnal_54 !=. | firstPosAnal_66 !=. | ///
								firstPosAnal_6873 !=. | firstPosAnal_70 !=. | ///
								firstPosAnal_74 !=.
			
	 *keep only keep the FIRST first positive (note: max 5 visits in H2M dataset)
	by persid: replace firstPos_LR =. if firstPos_LR[_n-1]==1 | ///
			  firstPos_LR[_n-2]==1 | firstPos_LR[_n-3]==1 | firstPos_LR[_n-4]==1 
			   
	*let the first positive indicate the roundnumber, and be carried forward
	replace firstPos_LR =round if firstPos_LR ==1 
	sort persid round firstPos_LR
	by persid: carryforward firstPos_LR if missing(firstPos_LR), replace

	***** Step 2: first negative of the grouping
	gen firstNeg_LR =.
	replace firstNeg_LR = 1 if firstNegAnal_6 !=. | firstNegAnal_11 !=. | ///
							   firstNegAnal_34 !=. | firstNegAnal_40 !=. | ///
							   firstNegAnal_42 !=. | firstNegAnal_43 !=. | ///
							   firstNegAnal_44 !=. | firstNegAnal_53 !=. | ///
							   firstNegAnal_54 !=. | firstNegAnal_66 !=. | ///
							   firstNegAnal_6873 !=. | firstNegAnal_70 !=. | ///
							   firstNegAnal_74 !=.
							   
	*keep only the FIRST first negative
	by persid: replace firstNeg_LR =. if firstNeg_LR[_n-1]==1 | ///
									  firstNeg_LR[_n-2]==1 | firstNeg_LR[_n-3]==1
									  
	*indicate roundnumber, don't carry forward 
	replace firstNeg_LR =round if firstNeg_LR ==1	

	***** Step 3: grouped event variable
	gen event_LR=.
	replace event_LR=1	if firstNeg_LR !=. & firstPos_LR !=.

***** PMO 

/*for incident infections requiring 1 neg test for clearance, PMO consists of:
	- the middays before the first positive visit
	- the middays before the first negative visit
	- the difvisdt between two positive visits	*/
	
	gen PDO_LR=.

	*middays before first positive
	replace PDO_LR = middays if round == firstPos_LR

	*difvisdt between two positive visits clearance:
	*i.e. the 1st pos has already occurred, but the 1st neg has not yet occurred
	sort persid round firstNeg_LR
	by persid: carryforward firstNeg_LR if missing(firstNeg_LR), replace	
	by persid: replace PDO_LR = difvisdt if round < firstNeg_LR ///
										   & round > firstPos_LR	
	*middays before the first negative
	replace PDO_LR = middays if round == firstNeg_LR

	*convert to PMO
	gen PMO_LR =.
	replace PMO_LR = PDO_LR / 30.5


*clean up datafile
foreach var of varlist firstPosAnal_6 - firstPosAnal_74{
	drop `var'
	}
foreach var of varlist firstNegAnal_6 - firstNegAnal_74{
	drop `var'
	}

drop firstPos_Any firstNeg_Any PDO_Any
drop firstPos_HR firstNeg_HR PDO_HR
drop firstPos_LR firstNeg_LR PDO_LR
drop firstPos_2v firstNeg_2v PDO_2v
drop firstPos_4v firstNeg_4v PDO_4v
drop firstPos_9v firstNeg_9v PDO_9v

end 

save "H2M_ApproachB1.dta", replace

******************* Chapter 2.5 : preparing datasets for B2 and B3 ********************

*****************************************
/*The B2 approach counts a maximum of one clearance event. When the
first specific HPV type has cleared, no more person-time is counted.
Only prevalent infections are included.
The source dataset is "Dta\H2M_prevalent.dta"*/
*****************************************

/*To prepare a dataset for approach B2, we can use the same code as we used for
approach subset B1. The only difference between these approaches is that B1 
includes incident infections, and B2 includes prevalent infections. So, we should
start with a different source dataset (which censored non-prevalent infections). 

Then, the only differences are:
> in counting person-time: the middays before the first positive visit should not be
 counted (given that there was no visit prior to the baseline visit, upon which 
 the middate can be based). In our dataset, the 'middays' variable is thus empty 
 at visit 1. Therefore, we can copy the previously used code without a problem
> there is no need to remove non-first 'firstPos' values, given that all first positives
  will be at the baseline visit. Still, we can copy the previously used code without a problem
  
There are no problems with re-running the code we wrote for approach B1 to make
a dataset for approach B2. Therefore, we run the program 'makeappB' */

use "H2M_prevalent.dta", clear
makeappB
save "H2M_ApproachB2", replace

*****************************************
/*The B3 approach counts a maximum of one clearance event. When the
first specific HPV type has cleared, no more person-time is counted.
Both incident and baseline-prevalent infections are included.
The source dataset is "Dta\H2M_all.dta"*/
*****************************************

*Noting the same arguments as provided above, we can appplpy the program 'makeappB'

use "H2M_all.dta", clear
makeappB
save "H2M_ApproachB3", replace

******************* Chapter 2.6 : preparing datasets for C1 ********************

*****************************************
/*The C1 approach counts a clearance event when a person was positive for
at least one HPV type within a grouping, and thereafter is positive for zero 
types from that grouping. Still, multiple clearances per person per group are
allowed. Person-time is counted only once per person.
Only incident infections are included.
The source dataset is "Dta\H2M_incident.dta"

The code we use to make the dataset for approach C1 can also be used to make
the datasets for approach C2 and C3. Therefore, we write this into a program 
called makeappC*/
*****************************************

use "H2M_incident.dta", clear

program makeappC

***** Grouping: Any HPV 
*including all 25 types identified by the LiPa-25 
sort persid sampledat
gen Any_HPV=0
label variable Any_HPV "infection with any type of the Any HPV grouping; binary outcome (0: no, zero HPV types; 1: yes, >zero HPV types)"

replace Any_HPV=1 if Anal_6==1 | Anal_11==1 | Anal_16==1 | ///
	Anal_18==1 | Anal_31==1 | Anal_33==1 | Anal_34==1 | ///
	Anal_35==1 | Anal_39==1 | Anal_40==1 | Anal_42==1 | ///
	Anal_43==1 | Anal_44==1 | Anal_45==1 | Anal_51==1 | ///
	Anal_52==1 | Anal_53==1 | Anal_54==1 | Anal_56==1 | ///
	Anal_58==1 | Anal_59==1 | Anal_66==1 | Anal_6873==1 | ///
	Anal_70==1 | Anal_74==1

*censor observations that become positive for the grounping at the last visit
replace Any_HPV=. if round == Nrounds & Any_HPV==1 & Any_HPV[_n-1]==0

**creating a variable for clearance
gen event_Any=.
label variable event_Any "clearance event for the Any HPV grouping"
by persid: replace event_Any=1 if Any_HPV==0 & Any_HPV[_n-1]==1

**creating a variable for person-time
*person time comprises middays before the first visit a person is positive for
*a grouping; total days between visits if a person remains positive for that
*grouping; and the middays before a person becomes negative for all types within
*that grouping.

sort persid sampledat

*part 1: first positive
	gen firstPosAny=.
	label variable firstPosAny "the visit-number of the first visit at which a subject is positive for the Any-HPV grouping. This can occur multiple times"
	by persid: replace firstPosAny= round if Any_HPV == 1 & Any_HPV[_n-1]==0
	by persid: carryforward firstPosAny if missing(firstPosAny), replace									

	*part 2: first negative	
	gen firstNegAny =.
	label variable firstNegAny "the visit-number of the first visit at which a subject is negative for an infection with the Any HPV grouping. This can occur multiple times"
	by persid: replace firstNegAny = round if Any_HPV ==0 & Any_HPV[_n-1]==1

	*part 3: generate variables that will count persontime
	gen PDO_Any=.
	label variable PDO_Any "the person-days of observation at risk for clearance of grouped Any HPV"
	replace PDO_Any = middays if round == firstPosAny
	by persid: replace PDO_Any = difvisdt if Any_HPV==1 & Any_HPV[_n-1]==1 
	replace PDO_Any = middays if round == firstNegAny

	*part 4: convert days to months
	gen PMO_Any =.
	replace PMO_Any = PDO_Any /30.5

*clean up
drop firstPosAny firstNegAny 
drop PDO_Any
		
		
*********Grouping: HR_HPV 
*the 12 HR types are: 16, 18, 31, 33, 35, 39, 45, 51, 52, 56, 58, 59

gen HR_HPV=0
label variable HR_HPV "infection with any type of the HR-HPV grouping; binary outcome (0: no, zero HR types; 1: yes, >zero HR types)"

replace HR_HPV=1 if Anal_16==1 | Anal_18==1 | Anal_31==1 | ///
	Anal_33==1 | Anal_35==1 | Anal_39==1 | Anal_45==1 | ///
	Anal_51==1 | Anal_52==1 | Anal_56==1 | Anal_58==1 | ///
	Anal_59==1 
		
*censor observations that become positive for the grounping at 
*the last visit
replace HR_HPV=. if round == Nrounds & HR_HPV==1 & HR_HPV[_n-1]==0

****** creating a variable for clearance events
gen event_HR=.
label variable event_HR "clearance event for the HR-HPV grouping"
by persid: replace event_HR=1 if HR_HPV==0 & HR_HPV[_n-1]==1

*****creating a variable for person-time
/*calculate the total person time; this is done once per person*/	
/*for incident infections requiring 1 neg test for clearance, PMO consists of:
	- the middays before the 1st positive visit
	- the difvisdt between two positive visits	
	- the middays before the first negative visit*/
	
*part 1: first positive
	gen firstPosHR=.
	label variable firstPosHR "the visit-number of the first visit at which a subject is positive for the HR-HPV grouping. This can occur multiple times"
	by persid: replace firstPosHR= round if HR_HPV == 1 & HR_HPV[_n-1]==0
	by persid: carryforward firstPosHR if missing(firstPosHR), replace									

	*part 2: first negative	
	gen firstNegHR =.
	label variable firstNegHR "the visit-number of the first visit at which a subject is negative for an infection with the HR HPV grouping. This can occur multiple times"
	by persid: replace firstNegHR = round if HR_HPV ==0 & HR_HPV[_n-1]==1

	*part 3: generate variables that will count persontime
	gen PDO_HR=.
	label variable PDO_HR "the person-days of observation at risk for clearance of grouped HR HPV"
	replace PDO_HR = middays if round == firstPosHR
	by persid: replace PDO_HR = difvisdt if HR_HPV==1 & HR_HPV[_n-1]==1 
	replace PDO_HR = middays if round == firstNegHR

	*part 4: convert days to months
	gen PMO_HR =.
	replace PMO_HR = PDO_HR /30.5

*clean up
drop firstPosHR firstNegHR 
drop PDO_HR


*******Grouping: LR HPV
*the 13 LR types are: 6,11,34,40,42,43,44,53,54,66,68/73,70,74

sort persid sampledat
gen LR_HPV=0
label variable LR_HPV "infection with any type of the LR-HPV grouping; binary outcome (0: no, zero LR types; yes: >zero LR types)"

replace LR_HPV=1 if Anal_6==1 | Anal_11==1 | Anal_34==1 | ///
	Anal_40==1 | Anal_42==1 | Anal_43==1 | Anal_44==1 | ///
	Anal_53==1 | Anal_54==1 | Anal_66==1 | Anal_6873==1 | ///
	Anal_70==1 | Anal_74==1

*censor observations that become positive for the grounping at 
*the last visit
replace LR_HPV=. if round == Nrounds & LR_HPV==1 & LR_HPV[_n-1]==0

****** creating a variable for clearance events
gen event_LR=.
label variable event_LR "clearance event for the LR-HPV grouping"
by persid: replace event_LR=1 if LR_HPV==0 & LR_HPV[_n-1]==1

*****creating a variable for person-time
/*calculate the total person time; this is done once per person*/	
/*for incident infections requiring 1 neg test for clearance, PMO consists of:
	- the middays before the 1st positive visit
	- the difvisdt between two positive visits	
	- the middays before the first negative visit*/
	
*part 1: first positive
	gen firstPosLR=.
	label variable firstPosLR "the visit-number of the first visit at which a subject is positive for the LR-HPV grouping. This can occur multiple times"
	by persid: replace firstPosLR= round if LR_HPV == 1 & LR_HPV[_n-1]==0
	by persid: carryforward firstPosLR if missing(firstPosLR), replace									

	*part 2: first negative	
	gen firstNegLR =.
	label variable firstNegLR "the visit-number of the first visit at which a subject is negative for an infection with the LR HPV grouping. This can occur multiple times"
	by persid: replace firstNegLR = round if LR_HPV ==0 & LR_HPV[_n-1]==1

	*part 3: generate variables that will count persontime
	gen PDO_LR=.
	label variable PDO_LR "the person-days of observation at risk for clearance of grouped LR HPV"
	replace PDO_LR = middays if round == firstPosLR
	by persid: replace PDO_LR = difvisdt if LR_HPV==1 & LR_HPV[_n-1]==1 
	replace PDO_LR = middays if round == firstNegLR

	*part 4: convert days to months
	gen PMO_LR =.
	replace PMO_LR = PDO_LR /30.5

*clean up
drop firstPosLR firstNegLR 
drop PDO_LR

			
******* Grouping: 2v HPV 
*2v vaccine types include 16 & 18
sort persid sampledat
gen bival_HPV=0
label variable bival_HPV "infection with any type of the bivalent vaccine HPV grouping; binary outcome (0: no, zero 2v types; 1:yes, >zero 2v types)"
replace bival_HPV=1 if Anal_16==1 | Anal_18==1 
		
*censor observations that become positive for the grounping at 
*the last visit
replace bival_HPV=. if round == Nrounds & bival_HPV==1 & bival_HPV[_n-1]==0

****** creating a variable for clearance events
gen event_2v=.
label variable event_2v "clearance event for the 2v-HPV grouping"
by persid: replace event_2v=1 if bival_HPV==0 & bival_HPV[_n-1]==1

*****creating a variable for person-time
/*calculate the total person time; this is done once per person*/	
/*for incident infections requiring 1 neg test for clearance, PMO consists of:
	- the middays before the 1st positive visit
	- the difvisdt between two positive visits	
	- the middays before the first negative visit*/
	
*part 1: first positive
	gen firstPos2v=.
	label variable firstPos2v "the visit-number of the first visit at which a subject is positive for the 2v-HPV grouping. This can occur multiple times"
	by persid: replace firstPos2v= round if bival_HPV == 1 & bival_HPV[_n-1]==0
	by persid: carryforward firstPos2v if missing(firstPos2v), replace									

	*part 2: first negative	
	gen firstNeg2v =.
	label variable firstNeg2v "the visit-number of the first visit at which a subject is negative for an infection with the 2v HPV grouping. This can occur multiple times"
	by persid: replace firstNeg2v = round if bival_HPV ==0 & bival_HPV[_n-1]==1

	*part 3: generate variables that will count persontime
	gen PDO_2v=.
	label variable PDO_2v "the person-days of observation at risk for clearance of grouped 2v HPV"
	replace PDO_2v = middays if round == firstPos2v
	by persid: replace PDO_2v = difvisdt if bival_HPV==1 & bival_HPV[_n-1]==1 
	replace PDO_2v = middays if round == firstNeg2v

	*part 4: convert days to months
	gen PMO_2v =.
	replace PMO_2v = PDO_2v /30.5

*clean up
drop firstPos2v firstNeg2v 
drop PDO_2v
	

****Grouping: 4v HPV 
*the 4v vaccine types include 6, 11, 16 & 18

sort persid sampledat
gen quadri_HPV=0
label variable quadri_HPV "infection with any type of the quadrivalent vaccine HPV grouping; binary outcome (0: no, zero 4v types; 1: yes, >zero 4v types)"

replace quadri_HPV=1 if Anal_6 ==1 | Anal_11 ==1 | Anal_16==1 | ///
	Anal_18==1 

*censor observations that become positive for the grounping at 
*the last visit
replace quadri_HPV=. if round == Nrounds & quadri_HPV==1 & quadri_HPV[_n-1]==0

****** creating a variable for clearance events
gen event_4v=.
label variable event_4v "clearance event for the 4v-HPV grouping"
by persid: replace event_4v=1 if quadri_HPV==0 & quadri_HPV[_n-1]==1

*****creating a variable for person-time
/*calculate the total person time; this is done once per person*/	
/*for incident infections requiring 1 neg test for clearance, PMO consists of:
	- the middays before the 1st positive visit
	- the difvisdt between two positive visits	
	- the middays before the first negative visit*/
	
*part 1: first positive
	gen firstPos4v=.
	label variable firstPos4v "the visit-number of the first visit at which a subject is positive for the 4v-HPV grouping. This can occur multiple times"
	by persid: replace firstPos4v= round if quadri_HPV == 1 & quadri_HPV[_n-1]==0
	by persid: carryforward firstPos4v if missing(firstPos4v), replace									

	*part 2: first negative	
	gen firstNeg4v =.
	label variable firstNeg4v "the visit-number of the first visit at which a subject is negative for an infection with the 4v HPV grouping. This can occur multiple times"
	by persid: replace firstNeg4v = round if quadri_HPV ==0 & quadri_HPV[_n-1]==1

	*part 3: generate variables that will count persontime
	gen PDO_4v=.
	label variable PDO_4v "the person-days of observation at risk for clearance of grouped 4v HPV"
	replace PDO_4v = middays if round == firstPos4v
	by persid: replace PDO_4v = difvisdt if quadri_HPV==1 & quadri_HPV[_n-1]==1 
	replace PDO_4v = middays if round == firstNeg4v

	*part 4: convert days to months
	gen PMO_4v =.
	replace PMO_4v = PDO_4v /30.5

*clean up
drop firstPos4v firstNeg4v 
drop PDO_4v

	
******Grouping 9v HPV 
*the 9v vaccine types include: 6, 11, 16, 18, 31, 33, 45, 52 and 58 

sort persid sampledat
gen nine_HPV=0
label variable nine_HPV "infection with any type of the nine-valent vaccine HPV grouping; binary outcome (0: no, zero 9v types; 1: yes, >zero 9v types)"

replace nine_HPV=1 if Anal_6 ==1 | Anal_11 ==1 | Anal_16==1 | ///
						Anal_18==1 | Anal_31==1 | Anal_33 ==1 | ///
						Anal_45==1 | Anal_52 ==1 | Anal_58 ==1

*censor observations that become positive for the grounping at 
*the last visit
replace nine_HPV=. if round == Nrounds & nine_HPV==1 & nine_HPV[_n-1]==0

****** creating a variable for clearance events
gen event_9v=.
label variable event_9v "clearance event for the 9v-HPV grouping"
by persid: replace event_9v=1 if nine_HPV==0 & nine_HPV[_n-1]==1

*****creating a variable for person-time
/*calculate the total person time; this is done once per person*/	
/*for incident infections requiring 1 neg test for clearance, PMO consists of:
	- the middays before the 1st positive visit
	- the difvisdt between two positive visits	
	- the middays before the first negative visit*/
	
*part 1: first positive
	gen firstPos9v=.
	label variable firstPos9v "the visit-number of the first visit at which a subject is positive for the 9v-HPV grouping. This can occur multiple times"
	by persid: replace firstPos9v= round if nine_HPV == 1 & nine_HPV[_n-1]==0
	by persid: carryforward firstPos9v if missing(firstPos9v), replace									

	*part 2: first negative	
	gen firstNeg9v =.
	label variable firstNeg9v "the visit-number of the first visit at which a subject is negative for an infection with the 9v HPV grouping. This can occur multiple times"
	by persid: replace firstNeg9v = round if nine_HPV ==0 & nine_HPV[_n-1]==1

	*part 3: generate variables that will count persontime
	gen PDO_9v=.
	label variable PDO_9v "the person-days of observation at risk for clearance of grouped 9v HPV"
	replace PDO_9v = middays if round == firstPos9v
	by persid: replace PDO_9v = difvisdt if nine_HPV==1 & nine_HPV[_n-1]==1 
	replace PDO_9v = middays if round == firstNeg9v

	*part 4: convert days to months
	gen PMO_9v =.
	replace PMO_9v = PDO_9v /30.5

*clean up
drop firstPos9v firstNeg9v 
drop PDO_9v

end
		
save "H2M_ApproachC1", replace	

******************* Chapter 2.7 : preparing datasets for C2 and C3  ********************

/*The C2 approach counts a clearance event when a person was positive for
at least one HPV type within a grouping, and thereafter is positive for zero 
types from that grouping. Still, multiple clearances per person per group are
allowed. Person-time is counted only once per person.
Only prevalent infections are included.
The source dataset is "Dta\H2M_prevalent.dta" */

use "H2M_prevalent.dta", clear
makeappC
save "Dta\H2M_ApproachC2", replace

/*The C3 approach counts a clearance event when a person was positive for
at least one HPV type within a grouping, and thereafter is positive for zero 
types from that grouping. Still, multiple clearances per person per group are
allowed. Person-time is counted only once per person.
Both incident and baseline-prevalent infections are included.
The source dataset is "Dta\H2M_prevalent.dta" */

use "H2M_all.dta", clear
makeappC
save "Dta\H2M_ApproachC3", replace

	
********************************************************************************
*** Chapter 3: CR analyses
********************************************************************************

**************** Chapter 3.1: CRs for approaches A1, A2 and A3 *****************
/*note that for Approach set A, multiple clearance events can occur per person,
even at the same timepoint. To account for clustering of clearance observations
within individuals, we model a constant-only generalised linear model and account
for clustering of observations on the 'persid'-level.*/

/*since we repeat this code for approach A2 and A3, and C1, C2 and C3, we make
a program out of this code, called 'analysisAC' */

**Approach A1
use "H2M_ApproachA1.dta", clear 

program analysisAC 

gen period =1 

	* Any HPV
	gen period=1
	stset PMO_Any, f(events_Any) 
	glm _d period, exposure(PMO_Any) link(log) f(poisson) eform vce(cluster persid)

	*HR HPV
	stset PMO_HR, f(events_HR) 
	glm _d period, exposure(PMO_HR) link(log) f(poisson) eform vce(cluster persid)

	*LR
	stset PMO_LR, f(events_LR) 
	glm _d period, exposure(PMO_LR) link(log) f(poisson) eform vce(cluster persid)

	*2-valent
	stset PMO_2v, f(events_2v) 
	glm _d period, exposure(PMO_2v) link(log) f(poisson) eform vce(cluster persid)

	*4-valent
	stset PMO_4v, f(events_4v) 
	glm _d period, exposure(PMO_4v) link(log) f(poisson) eform vce(cluster persid)

	*9-valent
	stset PMO_9v, f(events_9v) 
	glm _d period, exposure(PMO_9v) link(log) f(poisson) eform vce(cluster persid)

	end

analysisAC 

**Approach A2
use "H2M_ApproachA2.dta", clear 
analysisAC

**Approach A3
use "H2M_ApproachA3.dta", clear 
analysisAC

**************** Chapter 3.2: CRs for approaches B1, B2 and B3 *****************
/*note that for Approach set B, only one clearance event can occur per person.
even at the same timepoint. Observations can therefore not be clustered within
individuals. We therefore do not use 'vce(cluster persid)'.*/

/*since we repeat the code below for approach B2 and B3, we make
a program out of this code, called 'analysisB' */

use "H2M_ApproachB1.dta", clear

program analysisB

gen period=1

	* Any HPV
	stset PMO_Any, failure(event_Any==1)
	glm _d period, exposure(PMO_Any) link(log) f(poisson) eform

	*HR
	stset PMO_HR, failure(event_HR==1)
	glm _d period, exposure(PMO_HR) link(log) f(poisson) eform

	*LR
	stset PMO_LR, failure(event_LR==1)
	glm _d period, exposure(PMO_LR) link(log) f(poisson) eform

	*2-valent
	stset PMO_2v, failure(event_2v==1)
	glm _d period, exposure(PMO_2v) link(log) f(poisson) eform

	*4-valent
	stset PMO_4v, failure(event_4v==1)
	glm _d period, exposure(PMO_4v) link(log) f(poisson) eform

	*9-valent
	stset PMO_9v, failure(event_9v==1)
	glm _d period, exposure(PMO_9v) link(log) f(poisson) eform

end

analysisB

**Approach B2
use "H2M_ApproachB2.dta", clear 
analysisB

**Approach B3
use "H2M_ApproachB3.dta", clear 
analysisB

**************** Chapter 3.3: CRs for approaches C1, C2 and C3 *****************
/*note that for Approach set A, multiple clearance events can occur per person.
To account for clustering of clearance observations within individuals, 
we model a constant-only generalised linear model and account for clustering of 
observations on the 'persid'-level.*/

**Approach C1
use "H2M_ApproachC1.dta", clear 
analysisAC

**Approach C2
use "H2M_ApproachC2.dta", clear 
analysisAC

**Approach C3
use "H2M_ApproachC3.dta", clear 
analysisAC
