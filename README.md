# Designing isolation guidelines for COVID-19 patients with rapid antigen tests

Yong Dam Jeong, Keisuke Ejima, Kwang Su Kim, Woo Joohyeon, Shoya Iwanami, Yasuhisa Fujita, Il Hyo Jung, Kazuyuki Aihara, Kenji Shibuya, Shingo Iwami, Ana I. Bento and Marco Ajelli


## 1) Viral dynamics model

'symptomatic109.csv' and 'asymptomatic101.csv' are viral load data for symptomatic and asymptomatic patients, respectively. The data were extracted from published literatures.

'model.txt' is the description for mathematical model of viral dynamics.

These files are used to estimate parameters of the viral dynamics model in MONOLIX2019R2.


## 2) Simulation for isolation guideline

'populationParameters_sym.txt' and 'populationParameters_asym.txt' are the estimated parameters of viral dynamics model for symptomatic and asymptomatic patients, respectively.

'Symptomatic.R' and 'Asymptomatic.R' calculate numerical values of "Probability of prematurely ending isolation" and "Length of unnecessarily prolonged isolation" for symptomatic and asymptomatic patients, respectively.

* Text files and R codes above should be in the same location.
* You can change conditions "Detection limit" and "Infectiousness threshold" in the R codes. 
