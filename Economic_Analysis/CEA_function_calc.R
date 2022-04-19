
CEA <- function( input.list,case.files) {
  
  
  #-----------------------------------------------------------------------------------------------------------------------------------
  # Call in or input data and prepare it for analysis
  #-----------------------------------------------------------------------------------------------------------------------------------
  
  # unpack input list to create free objects (Economic Inputs)
  for (i in 1:length(input.list)) assign(names(input.list)[i],input.list[[i]])
  
  # Model Inputs
  # read in the baseline and chosen scenario incidence over time data
  # Series starts from 1st Jan 2021 to 30th June 2022: A matrix of 545 by 6 (age groups) :"0-19","20-49","50-59","60-69","70-79","80+"
  
  # Cases: -
  #- index scenario cases
  asymp_cases <- case.files[[1]]
  mild_cases <- case.files[[2]]
  severe_cases <- case.files[[3]]
  critical_cases <- case.files[[4]]
  deaths <- case.files[[5]]
  num.doses.procured <-  case.files[[6]]
  
  #- baseline scenario cases
  asymp_cases.baseline <-  case.files[[7]]
  mild_cases.baseline <-  case.files[[8]]
  severe_cases.baseline <-  case.files[[9]]
  critical_cases.baseline <-  case.files[[10]]
  deaths.baseline <- case.files[[11]]
  num.doses.procured.baseline <-  case.files[[12]]
  
  
  
  GDPdef <- read.csv(file = "R_models/Kenya_GDPdeflator_Data.csv",header=TRUE, sep=',')
  # -----------------------------------------------------------------------------------------------------------------------------------
  
  # DEFINE SOME KEY vaccination inputs, those that DONT vary and you wont need to so sensitivity analysis on here
  #doses=2
  
  #the following are the unit costs per dose for the sobol function
  #Vaccine delivery costs (2021 costs)
  #vaccine.delivery.unitcost_novacc <- 0
  #vaccine.delivery.unitcost_30 <-6.11 #  range of costs +/- 20% (5.99 - 6.23)
  #vaccine.delivery.unitcost_50 <-4.16 #  range of costs +/- 20% (4.08 - 4.23)
  #vaccine.delivery.unitcost_70 <-3.90 #  range of costs +/- 20% (3.82 - 3.98)
  
  #Treatment costs are patient costs per day (2020 costs)
  #unit.treat.cost.asymp <- 19 * 0.52% = 0.099   # Applying case detection rate of 0.52% #range of costs +/- 20% (0.079 - 0.119)
  #unit.treat.cost.mild <- 19  * 0.52% = 0.099   # Applying case detection rate of 0.52% #range of costs +/- 20% 
  #unit.treat.cost.severe <- 124.54 # range of costs +/- 20% (99.63 - 149.45)
  #unit.treat.cost.critical <- 599.51 # range of costs +/- 20% (479.61 - 719.41)
  
  
  
  #-----------------------------------------------------------------------------------------------------------------------------------
  # Calculations: Vaccination costs
  #----------------------------------------------------------------------------------------------------------------------------------
  #No vaccination costs for baseline-as it is the no vaccination.
  
  # Total vaccine costs over time (columns), iterations (rows)
  total.vaccine.cost <- num.doses.procured*vaccine.cost.per.dose
  total.vaccine.cost.baseline <- num.doses.procured.baseline*
    vaccine.cost.per.dose
  #Total supplies costs over time (columns), iterations (rows)
  total.supplies.cost <- num.doses.procured*supplies.cost.per.dose
  total.supplies.cost.baseline <- num.doses.procured.baseline*
                                  supplies.cost.per.dose
  # Delivery costs over time, iterations
  vaccine.delivery.cost <- num.doses.procured*vaccine.delivery.unitcost
  vaccine.delivery.cost.baseline <- num.doses.procured.baseline*
                                  vaccine.delivery.unitcost.baseline
  
  total.vacc.cost <- total.vaccine.cost + total.supplies.cost + vaccine.delivery.cost
  Net.total.vacc.cost <- sum(total.vacc.cost)
  total.vacc.cost.baseline <- total.vaccine.cost.baseline + 
                            total.supplies.cost.baseline + 
                            vaccine.delivery.cost.baseline
  Net.total.vacc.cost.baseline <- sum(total.vacc.cost.baseline)
  
  #--------------------------------------------------------------------------------------------------------------------------------------
  # Calculations: Treatment costs by COVID severity
  #--------------------------------------------------------------------------------------------------------------------------------------
  # Asymptomatic treatment costs
  asymp_cases.total <- sum(asymp_cases, na.rm=T)
  asymp.cost  <- asymp_cases.total * unit.treat.cost.asymp * dur.asymptomatic.covid *365
  
  asymp_cases.baseline.total <- sum(asymp_cases.baseline, na.rm=T)
  asymp.cost.baseline <- asymp_cases.baseline.total * unit.treat.cost.asymp * dur.asymptomatic.covid*365
  
  # Mild treatment costs
  mild_cases.total <- sum(mild_cases, na.rm=T)
  mild.cost  <- mild_cases.total * unit.treat.cost.mild * dur.mild.covid*365
  
  mild_cases.baseline.total <- sum(mild_cases.baseline, na.rm = T)
  mild.cost.baseline <- mild_cases.baseline.total * unit.treat.cost.mild * dur.mild.covid*365
  
  # Severe disease treatment costs
  severe_cases.total <- sum(severe_cases, na.rm=T)
  severe.cost <- severe_cases.total * unit.treat.cost.severe * length.hosp.severe*365
  
  severe_cases.baseline.total <- sum(severe_cases.baseline, na.rm=T)
  severe.cost.baseline <- severe_cases.baseline.total * unit.treat.cost.severe * length.hosp.severe*365
  
  # Critical disease treatment costs
  critical_cases.total <- sum(critical_cases, na.rm=T)
  critical.ICU.cost <- critical_cases.total * unit.treat.cost.critical *length.ICU.critical*365
  critical.hospital.cost <- critical_cases.total * unit.treat.cost.severe *length.hosp.severe*365
  critical.cost <- critical.ICU.cost + critical.hospital.cost
  
  critical_cases.baseline.total <- sum(critical_cases.baseline, na.rm=T)
  critical.ICU.baseline.cost <- critical_cases.baseline.total * unit.treat.cost.critical * length.ICU.critical*365
  critical.hospital.baseline.cost <- critical_cases.baseline.total *unit.treat.cost.severe *length.hosp.severe*365
  critical.cost.baseline <- critical.ICU.baseline.cost + critical.hospital.baseline.cost
  
  
  # Total treatment costs across across all syndrome and age groups
  Net.total.treatment.cost <- asymp.cost + mild.cost +severe.cost + critical.cost
  Net.total.treatment.cost_baseline <- asymp.cost.baseline + mild.cost.baseline + severe.cost.baseline + critical.cost.baseline
  Net.total.treatment.hospital.cost <- severe.cost + critical.cost
  Net.total.treatment.cost.nonhospital <- asymp.cost + mild.cost
  Net.total.treatment.hospital.cost_baseline <- severe.cost.baseline + critical.cost.baseline
  Net.total.treatment.cost.nonhospital_baseline <- asymp.cost.baseline + mild.cost.baseline
  
  #----------------------------------------------------------------------------------------------------------------------------------------------------------
  # Calculations: Disease burden
  #-----------------------------------------------------------------------------------------------------------------------------------------------------------
  #---------------------------------------------------------------------------------------------------------------
  
  # Deaths from COVID  disease
  # Deaths are not defined by symptom just a general infection-fatality-rate
  # These are not observable deaths rather predicted deaths not accounting for recording/observation probability
  total.deaths_across.age<- cbind(deaths, total=rowSums(deaths))  # sum over age group
  total.deaths_across.age <- subset(total.deaths_across.age,select = total)
  Net.total.deaths = sum(total.deaths_across.age ,na.rm=T)
  
  total.deaths_across.age.baseline<- cbind(deaths.baseline, total=rowSums(deaths.baseline))  # sum over age group
  total.deaths_across.age.baseline <- subset(total.deaths_across.age.baseline,select = total)
  Net.total.deaths_baseline = sum(total.deaths_across.age.baseline ,na.rm=T)
  
  
  #Cases across syndrome
  cases <- asymp_cases  + mild_cases  + severe_cases  + critical_cases
  cases_across.age <- cbind(cases, total=rowSums(cases))# sum over age group
  cases_across.age <- subset(cases_across.age, select = total)
  Net.total.cases = sum(cases_across.age,na.rm=T)
  
  cases.baseline <- asymp_cases.baseline + mild_cases.baseline + severe_cases.baseline + critical_cases.baseline
  cases.baseline_across.age <- cbind(cases.baseline, total=rowSums(cases.baseline)) # sum over age group
  cases.baseline_across.age <- subset(cases.baseline_across.age,select = total)
  Net.total.cases_baseline = sum(cases.baseline_across.age,na.rm=T)
  
  #Cases by hospitalization/non-hospitalization across syndrome
  hospital_cases <- severe_cases + critical_cases
  hospital_cases_across.age <-  cbind(hospital_cases, total=rowSums(hospital_cases))
  hospital_cases_across.age <- subset(hospital_cases_across.age,select = total)
  Net.total.hospital.cases = sum(hospital_cases_across.age,na.rm=T)
  
  nonhospital_cases <- asymp_cases  + mild_cases 
  nonhospital_cases_across.age <-  cbind(nonhospital_cases, total=rowSums(nonhospital_cases))
  nonhospital_cases_across.age <- subset(nonhospital_cases_across.age,select = total)
  Net.total.nonhospital.cases = sum(nonhospital_cases_across.age,na.rm=T)
  
  hospital_cases.baseline <- severe_cases.baseline + critical_cases.baseline
  hospital_cases_across.age.baseline <-  cbind(hospital_cases.baseline, total=rowSums(hospital_cases.baseline))
  hospital_cases_across.age.baseline <- subset(hospital_cases_across.age.baseline,select = total)
  Net.total.hospital.cases_baseline = sum(hospital_cases_across.age.baseline,na.rm=T)
  
  nonhospital_cases.baseline <- asymp_cases.baseline  + mild_cases.baseline
  nonhospital_cases_across.age.baseline <-  cbind(nonhospital_cases.baseline, total=rowSums(nonhospital_cases.baseline))
  nonhospital_cases_across.age.baseline <- subset(nonhospital_cases_across.age.baseline,select = total)
  Net.total.nonhospital.cases_baseline = sum(nonhospital_cases_across.age.baseline,na.rm=T)
  
  
  # Averted cases by syndrome and deaths
  averted.cases.asymp  <-  asymp_cases.baseline - asymp_cases
  averted.cases.mild  <-  mild_cases.baseline - mild_cases
  averted.cases.severe  <-  severe_cases.baseline - severe_cases
  averted.cases.critical <-  critical_cases.baseline - critical_cases
  averted.cases_across.age <-  cases.baseline_across.age - cases_across.age
  
  averted.deaths <- deaths.baseline - deaths
  averted.deaths_across.age <- cbind(averted.deaths, total=rowSums(averted.deaths))# sum over age group
  averted.deaths_across.age <- subset(averted.deaths_across.age,select = total)
  
  # -----------------------------------------------------------------------------------------------------------------------------------
  
  # Disability weights, life expectancy -------
  #Age group in  current model: "0-19","20-49","50-59","60-69","70-79","80+")
  # The mean AGES to reflect the mean ages of the above
  mean.age <- c(9.27, 31.75, 54.10, 63.85, 73.41, 86.00)        # Weighted Mean age in each the the six age groups (years) updated from KNBS
  life.expectancy <- c(64.10,40.94,24.71,17.64,11.41,4.84)      # The Kenya standard life table for age groups (abridged) in 2019 (WHO) was used, for death at the mean age (rounded to 9, 32, 54, 64, 73, 86) of each of the six age groups
  # Previous ref: Kenya standard life table for single age years in 2011 (WHO)
  
  #TIME
  time_days <- 1:547
  time <- time_days/365
  
  
  #---------------------------------------------------------------------------------------------------------------------------------------
  # Calculations: DALY - Disability-Adjusted Life Year
  #---------------------------------------------------------------------------------------------------------------------------------------
  #effect.discount=0.03
  # Discounting function for YLL and YLD
  # K is a modulating factor equaling one if age weighting is applied and zero otherwise; A and L represent, respectively, the age at onset and the duration
  # 'r' the discount rate, 'x' the concerned age, and 'a' the age to which the burden will be assigned
  effect.discount.fun <- function(x,K,C=0.1658,beta=0.04,r,a){
    K*C*x*exp(-beta*x)*exp(-r*(x-a)) + (1-K)*exp(-r*(x-a))
  }
  discount.burden<- function(A,L,K,r,a) {integrate(f=effect.discount.fun, lower=A,upper=A+L, K=K,r=r,a=a)$value}
  
  #references for calculating disability weights
  #https://www.ssph-journal.org/articles/10.3389/ijph.2021.619011/full
  #https://www.burden-eu.net/docs/covid19-bod-protocol.pdf 
  #https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8165085/
  #https://www.who.int/quantifying_ehimpacts/publications/en/9241546204chap3.pdf
  
  # Years of life lost due to premature mortality (YLL) ----------------------------------------------------------------------------------------------------------------------------------------------------------------
  YLL <-  mapply ('*', deaths,mapply(as.list(life.expectancy), FUN=function(x,y){discount.burden(A=y,L=x,K=0,r=effect.discount,a=y)}, as.list(mean.age),SIMPLIFY=F), SIMPLIFY =F)   # discount  YLL
  Net.total.YLL <- Reduce('+', YLL)
  
  YLL.baseline <- mapply ('*', deaths.baseline,mapply(as.list(life.expectancy), FUN=function(x,y){discount.burden(A=y,L=x,K=0,r=effect.discount,a=y)}, as.list(mean.age),SIMPLIFY=F), SIMPLIFY =F)   # discount  YLL baseline
  Net.total.YLL_baseline <- Reduce('+', YLL.baseline)
  
  
  # Years lived with disability (YLD) ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  YLD.asymp  <- mapply ('*', asymp_cases, mapply(as.list(dur.asymptomatic.covid), FUN=function(x,y){DW.asymptomatic*discount.burden(A=y,L=x,K=0,r=effect.discount,a=y)}, as.list(mean.age), SIMPLIFY=F),SIMPLIFY = F)
  YLD.mild  <- mapply ('*', mild_cases, mapply(as.list(dur.mild.covid), FUN=function(x,y){DW.mild.disease*discount.burden(A=y,L=x,K=0,r=effect.discount,a=y)}, as.list(mean.age), SIMPLIFY=F),SIMPLIFY = F)
  YLD.severe  <- mapply ('*', severe_cases, mapply(as.list(dur.severe.covid), FUN=function(x,y){DW.severe.disease*discount.burden(A=y,L=x,K=0,r=effect.discount,a=y)}, as.list(mean.age), SIMPLIFY=F),SIMPLIFY = F)
  YLD.critical  <- mapply ('*', critical_cases, mapply(as.list(dur.critical.covid), FUN=function(x,y){DW.critical.disease*discount.burden(A=y,L=x,K=0,r=effect.discount,a=y)}, as.list(mean.age), SIMPLIFY=F),SIMPLIFY = F)
  
  YLD.asymp.baseline  <- mapply ('*', asymp_cases.baseline, mapply(as.list(dur.asymptomatic.covid), FUN=function(x,y){DW.asymptomatic*discount.burden(A=y,L=x,K=0,r=effect.discount,a=y)}, as.list(mean.age), SIMPLIFY=F),SIMPLIFY = F)
  YLD.mild.baseline  <- mapply ('*', mild_cases.baseline, mapply(as.list(dur.mild.covid), FUN=function(x,y){DW.mild.disease*discount.burden(A=y,L=x,K=0,r=effect.discount,a=y)}, as.list(mean.age), SIMPLIFY=F),SIMPLIFY = F)
  YLD.severe.baseline  <- mapply ('*', severe_cases.baseline, mapply(as.list(dur.severe.covid), FUN=function(x,y){DW.severe.disease*discount.burden(A=y,L=x,K=0,r=effect.discount,a=y)}, as.list(mean.age), SIMPLIFY=F),SIMPLIFY = F)
  YLD.critical.baseline  <- mapply ('*', critical_cases.baseline, mapply(as.list(dur.critical.covid), FUN=function(x,y){DW.critical.disease*discount.burden(A=y,L=x,K=0,r=effect.discount,a=y)}, as.list(mean.age), SIMPLIFY=F),SIMPLIFY = F)
  
  #Total YLD
  Total.YLD.1 <- mapply ('+', YLD.asymp, YLD.mild, SIMPLIFY = F)
  Total.YLD.2 <- mapply ('+', Total.YLD.1, YLD.severe, SIMPLIFY = F)
  Total.YLD <- mapply ('+', Total.YLD.2, YLD.critical, SIMPLIFY = F)
  Net.total.YLD <- Reduce('+', Total.YLD)
  
  Total.YLD.baseline.1 <- mapply ('+', YLD.asymp.baseline, YLD.mild.baseline, SIMPLIFY = F)
  Total.YLD.baseline.2 <- mapply ('+', Total.YLD.baseline.1, YLD.severe.baseline, SIMPLIFY = F)
  Total.YLD.baseline <- mapply ('+', Total.YLD.baseline.2, YLD.critical.baseline, SIMPLIFY = F)
  Net.total.YLD_baseline <- Reduce('+', Total.YLD.baseline)
  
  # YLD by hospitalization/non-hospitalization
  YLD.hospital <- mapply ('+', YLD.severe,YLD.critical,SIMPLIFY = F)
  Net.total.YLD.hospital <- Reduce('+', YLD.hospital)
  YLD.nonhospital <- mapply ('+', YLD.asymp,YLD.mild,SIMPLIFY = F)
  Net.total.YLD.nonhospital <- Reduce('+', YLD.nonhospital)
  
  YLD.hospital.baseline <- mapply ('+', YLD.severe.baseline,YLD.critical.baseline,SIMPLIFY = F )
  Net.total.YLD.hospital_baseline <- Reduce('+', YLD.hospital.baseline)
  YLD.nonhospital.baseline <- mapply ('+', YLD.asymp.baseline,YLD.mild.baseline,SIMPLIFY = F)
  Net.total.YLD.nonhospital_baseline <- Reduce('+', YLD.nonhospital.baseline)
  
  
  # Total DALYs ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  Net.total.DALY <- Net.total.YLD + Net.total.YLL
  Net.total.DALY_baseline <- Net.total.YLD_baseline + Net.total.YLL_baseline
  
  
  Averted.DALY <-  Net.total.DALY_baseline - Net.total.DALY
  
  #------------------------------------------------------------------------------------------------------------------------------------
  # Calculations: Productivity losses costs
  #-----------------------------------------------------------------------------------------------------------------------------------
  testing.rate <- 0.52/100
  GDP.per.capita <- 1838.21 #per annum
  informal.sector <- 80/100 #Of those who isolate for asymp/mild disease, those likely to be unproductive are those in the informal sector
  dur.quarantine <- 14/365 #14 days, use duration of quarantine if the duration of disease is less than 14days(qurantine period)
  
  #Productivity losses due to morbidity
  prod.cost.asymp <- asymp_cases.total* testing.rate * GDP.per.capita * dur.quarantine * informal.sector
  prod.cost.mild <- mild_cases.total * testing.rate * GDP.per.capita * dur.quarantine * informal.sector
  prod.cost.severe <- severe_cases.total * GDP.per.capita * dur.quarantine
  prod.cost.critical <- critical_cases.total * GDP.per.capita * dur.critical.covid
  #sum across syndrome
  prod.cost.cases <- prod.cost.asymp + prod.cost.mild + prod.cost.severe + prod.cost.critical 
  
  
  prod.cost.asymp.baseline <- asymp_cases.baseline.total * testing.rate * GDP.per.capita * dur.quarantine * informal.sector
  prod.cost.mild.baseline <- mild_cases.baseline.total * testing.rate * GDP.per.capita * dur.quarantine * informal.sector
  prod.cost.severe.baseline <- severe_cases.baseline.total * GDP.per.capita * dur.quarantine
  prod.cost.critical.baseline <- critical_cases.baseline.total * GDP.per.capita * dur.critical.covid
  #sum across syndrome
  prod.cost.cases.baseline <- prod.cost.asymp.baseline + prod.cost.mild.baseline + prod.cost.severe.baseline + prod.cost.critical.baseline
  
  #Productivity losses due to mortality
  prod.costs.deaths <- Net.total.YLL * GDP.per.capita
  prod.costs.deaths.baseline <- Net.total.YLL_baseline * GDP.per.capita
  
  #Total Productivity Costs
  Net.total.prod.costs <- prod.cost.cases + prod.costs.deaths
  Net.total.prod.costs_baseline <- prod.cost.cases.baseline + prod.costs.deaths.baseline
  
  
  #------------------------------------------------------------------------------------------------------------------------------------
  # Calculations: outputs
  #-------------------------------------------------------------------------------------------------------------------------------------
  # Net cost index scenario (this includes the vaccination and treatment costs)
  Total.Net.cost <- Net.total.vacc.cost + Net.total.treatment.cost + Net.total.prod.costs
  
  # Net cost baseline (In this case, baseline vaccination costs were zero)
  Total.Net.cost_baseline <- Net.total.treatment.cost_baseline + Net.total.prod.costs_baseline
  
  
  # Increment cost
  Increment.cost <-  Total.Net.cost - Total.Net.cost_baseline
  
  # return(list(Total.Net.cost=Total.Net.cost,
  #             Total.Net.cost_baseline = Total.Net.cost_baseline,
  #             Net.total.DALY = Net.total.DALY
  # 
  return(list(Total.Net.cost=Total.Net.cost,
              Total.Net.cost_baseline = Total.Net.cost_baseline,
              Net.total.DALY = Net.total.DALY,
              Net.total.DALY_baseline = Net.total.DALY_baseline,
              Net.total.vacc.cost = Net.total.vacc.cost,
              Net.total.vacc.cost.baseline = Net.total.vacc.cost.baseline,
              Net.total.treatment.cost = Net.total.treatment.cost,
              Net.total.treatment.cost_baseline = Net.total.treatment.cost_baseline,
              Net.total.deaths = Net.total.deaths,
              Net.total.deaths_baseline = Net.total.deaths_baseline,
              Net.total.cases = Net.total.cases,
              Net.total.cases_baseline =Net.total.cases_baseline,
              Net.total.YLL = Net.total.YLL,
              Net.total.YLL_baseline = Net.total.YLL_baseline,
              Net.total.YLD = Net.total.YLD,
              Net.total.YLD_baseline = Net.total.YLD_baseline,
              Net.total.prod.costs = Net.total.prod.costs,
              Net.total.prod.costs_baseline = Net.total.prod.costs_baseline

  ))
  gc()
  
} # end function
