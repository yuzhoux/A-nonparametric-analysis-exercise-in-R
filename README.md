# A-nonparametric-analysis-exercise-in-R
This is an analysis of data from a study conducted to determine if a medication could be effective at reducing allergic symptoms.  These data represent symptom scores for a randomly selected day for each patient during a 10-week follow-up period.  Note that a placebo group, and three active dose groups were studied.  In addition, the current pollen count (the stimulus for symptoms) was also measured for the day that symptoms were collected.

Variables (in column order): <br>
obs      	patient ID <br>
itchy    	symptoms (1=none, 2=mild, 3=moderate, 4=severe) <br>
sneezy   	symptoms (1=none, 2=mild, 3=moderate, 4=severe) <br>
runny    	symptoms (1=none, 2=mild, 3=moderate, 4=severe) <br>
stuffy   	symptoms (1=none, 2=mild, 3=moderate, 4=severe) <br>
treat    	treatment (1 = placebo, 2 = 15 mg, 3 = 30 mg, 4 = 60 mg ) <br>
day      	study day <br>
pollen  	pollen count (particles per cubic meter) <br>
age      	age in years <br>
bmi		body mass index (BMI) <br>
Note: BMI between 25 and 30 is considered “overweight” and BMI over 30 is considered “obese” <br>
female     	gender (0=male, 1=female) <br>
logpollen	log10(pollen+1) <br>
