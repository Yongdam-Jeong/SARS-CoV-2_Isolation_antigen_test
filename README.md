# Designing isolation guidelines for COVID-19 patients with rapid antigen tests


# 1) Viral dynamics model

'symptomatic109.csv' and 'asymptomatic101.csv' are viral load data for symptomatic and asymptomatic patients, respectively. The data were extracted from published literature.

'model.txt' is the description for mathematical model of viral dynamics.

These files are used to estimate parameters of the viral dynamics model in MONOLIX2019R2.


# 2) Simulation for isolation guideline

'populationParameters_sym.txt' and 'populationParameters_asym.txt' are the estimated parameters of viral dynamics model for symptomatic and asymptomatic patients, respectively.

'Symptomatic.R' calculates numerical values of "Probability of prematurely ending isolation" and "Length of unnecessarily prolonged isolation" for symptomatic patients depicted in Figure 2, 3, and  4.

'Asymptomatic.R' calculates numerical values of "Probability of prematurely ending isolation" and "Length of unnecessarily prolonged isolation" for asymptomatic patients depicted in Figure 2, 3, and  4.

* Text files and R codes above should be in the same location.
* You can change conditions "Detection limit" and "Infectiousness threshold" in the R codes. 
