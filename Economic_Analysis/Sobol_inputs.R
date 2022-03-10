#--------------------------------------------------------------------------------------------
# Sampling cost and disease inputs
#-------------------------------------------------------------------------------------------

# Draws a Sobol Samples
Create.Sobol.inputs <- function(scenario=scenario,baseline.scenario=baseline.scenario,
                                discount.cost=T, discount.daly=T){
  
  SobolSamples <- sobol(n = 10, dim = 12, scrambling = 3)
  
  colnames(SobolSamples) <-c(
    "unit.treat.cost.asymp",   # 2inpatient treatment cost for asymptomatics
    "unit.treat.cost.mild",    # inpatient treatment cost for mild cases
    "unit.treat.cost.severe",  #outpatient treatment cost for severe cases
    "unit.treat.cost.critical",  #outpatient treatment cost for critical cases
    "vaccine.delivery.unitcost_30", # vaccine delivery cost per unit at 30% vaccination
    "vaccine.delivery.unitcost_50", # vaccine delivery cost per unit at 50% vaccination
    "vaccine.delivery.unitcost_70", # vaccine delivery cost per unit at 70% vaccination
    "DW.critical.disease", # disability weight of the critical disease
    "DW.severe.disease",  # disability weight of the severe disease
    "DW.mild.disease",    # disability weight of the mild disease#
    "length.hosp.severe", #length of hospital stay for severe cases
    "length.ICU.critical" #length of ICU stay for critical cases
  )
  
  #save(SobolSamples,file="SobolSamples")
  
  SobolSamples[,"unit.treat.cost.asymp"] <- qgamma(SobolSamples[,"unit.treat.cost.asymp"], shape=105.1285,scale=0.0009398019)
  SobolSamples[,"unit.treat.cost.mild"] <- qgamma(SobolSamples[,"unit.treat.cost.mild"], shape=105.1285,scale=0.0009398019)
  SobolSamples[,"unit.treat.cost.severe"] <- qgamma(SobolSamples[,"unit.treat.cost.severe"], shape=105.1124,scale=1.184827)
  SobolSamples[,"unit.treat.cost.critical"] <- qgamma(SobolSamples[,"unit.treat.cost.critical"], shape=105.1319,scale=5.702456)
  SobolSamples[,"vaccine.delivery.unitcost_30"] <- qgamma(SobolSamples[,"vaccine.delivery.unitcost_30"], shape=10055.05,scale=0.0006076546)
  SobolSamples[,"vaccine.delivery.unitcost_50"] <- qgamma(SobolSamples[,"vaccine.delivery.unitcost_50"], shape=13679.26,scale=0.00030411)
  SobolSamples[,"vaccine.delivery.unitcost_70"] <- qgamma(SobolSamples[,"vaccine.delivery.unitcost_70"], shape=9221.409,scale=0.0004229289)
  SobolSamples[,"length.hosp.severe"] <- qgamma(SobolSamples[,"length.hosp.severe"], shape=53.273947475158401,scale=0.0003599898843749761)
  SobolSamples[,"length.ICU.critical"] <- qgamma(SobolSamples[,"length.ICU.critical"], shape=53.273947475158401,scale=0.0003599898843749761)
  SobolSamples[,"DW.critical.disease"] <- qbeta(SobolSamples[,"length.ICU.critical"], 102.8285,54.16159)# beta(102.8285,54.16159)
  SobolSamples[,"DW.severe.disease"] <- qbeta(SobolSamples[,"DW.severe.disease"],  21.59726,140.7882)  # beta(21.59726,140.7882)
  SobolSamples[,"DW.mild.disease"] <- qbeta(SobolSamples[,"DW.mild.disease"], 21.29268,396.2109)  # beta(21.29268,396.2109)
  
  
  # GDP deflators
  # Source: International Monetary Fund, World Economic Outlook Database, October 2009
  #         https://www.imf.org/external/pubs/ft/weo/2016/01/weodata/weorept.aspx?pr.x=43&pr.y=6&sy=2000&ey=2021&scsm=1&ssd=1&sort=country&ds=.&br=1&c=664&s=NGDP_D&grp=0&a=
  GDPdef <- read.csv(file.path("R_models", "Kenya_GDPdeflator_Data.csv"),header=T,check.names=F)
  # Vaccine delivery costs
  vaccine.delivery.unitcost_30 <- SobolSamples[,"vaccine.delivery.unitcost_30"]
  vaccine.delivery.unitcost_50 <- SobolSamples[,"vaccine.delivery.unitcost_50"] 
  vaccine.delivery.unitcost_70 <- SobolSamples[,"vaccine.delivery.unitcost_70"]
  #obtain the specific vaccine delivery costs depending on the chosen scenario and baseline.scenario
  if (scenario == 2 && baseline.scenario == 1){
    vaccine.delivery.unitcost <- vaccine.delivery.unitcost_30
    vaccine.delivery.unitcost.baseline <- replicate(nrow(SobolSamples),0)
  } else if (scenario == 3 && baseline.scenario == 2){
    vaccine.delivery.unitcost <- vaccine.delivery.unitcost_50
    vaccine.delivery.unitcost.baseline <- vaccine.delivery.unitcost_30
  } else if (scenario == 4 && baseline.scenario == 3){
    vaccine.delivery.unitcost <- vaccine.delivery.unitcost_70
    vaccine.delivery.unitcost.baseline <- vaccine.delivery.unitcost_50
  } else if (scenario == 5 && baseline.scenario == 1){
    vaccine.delivery.unitcost <- vaccine.delivery.unitcost_30
    vaccine.delivery.unitcost.baseline <- replicate(nrow(SobolSamples),0)
  } else if (scenario == 6 && baseline.scenario == 5){
    vaccine.delivery.unitcost <- vaccine.delivery.unitcost_50
    vaccine.delivery.unitcost.baseline <- vaccine.delivery.unitcost_30
  } else if (scenario == 7 && baseline.scenario == 6){
    vaccine.delivery.unitcost <- vaccine.delivery.unitcost_70
    vaccine.delivery.unitcost.baseline <- vaccine.delivery.unitcost_50
  }
  
  # treatment costs
  unit.treat.cost.asymp <- SobolSamples[,"unit.treat.cost.asymp"]*GDPdef[,"2021"]/GDPdef[,"2020"]
  unit.treat.cost.mild <- SobolSamples[,"unit.treat.cost.mild"]*GDPdef[,"2021"]/GDPdef[,"2020"]
  unit.treat.cost.severe <- SobolSamples[,"unit.treat.cost.severe"]*GDPdef[,"2021"]/GDPdef[,"2020"]
  unit.treat.cost.critical <- SobolSamples[,"unit.treat.cost.critical"]*GDPdef[,"2021"]/GDPdef[,"2020"]
  
  #length of stay in hospital
  length.hosp.severe <-  SobolSamples[,"length.hosp.severe"]#The IQR range is 0.008219-0.02466 days which is 5/365(3-9 days)
  length.ICU.critical <- SobolSamples[,"length.ICU.critical"]#The IQR range is 0.01096-0.03014 days which is 7/365(4-11 days)
  
  # You can define fixed variable here
  create.inputs<-function(x, dscnt.cost=discount.cost, dscnt.daly=discount.daly){
    # Vaccine price per dose
    vaccine.cost.per.dose = 8.67  # economic procurement cost including 10% vaccine wastage rate
    supplies.cost.per.dose = 0.08 # economic procurement cost including 5% syringe wastage rate
    
    # Disability weights
    DW.critical.disease <-  SobolSamples[,"DW.critical.disease"]
    DW.severe.disease <- SobolSamples[,"DW.severe.disease"]
    DW.mild.disease   <- SobolSamples[,"DW.mild.disease"]
    DW.asymptomatic <- 0
    
    # Duration of disease in Years 
    dur.critical.covid <- 20/365
    dur.severe.covid <-   12/365
    dur.mild.covid <- 7/365
    dur.asymptomatic.covid <- 7/365
    # Discount rates
    if (dscnt.cost) {cost.discount <- 0.03} else {cost.discount <- 0}
    if (dscnt.daly) {effect.discount <- 0.03} else {effect.discount <- 0}
    return(list(
      scenario = scenario,
      baseline.scenario = baseline.scenario,
      vaccine.cost.per.dose = vaccine.cost.per.dose,
      supplies.cost.per.dose = supplies.cost.per.dose,
      vaccine.delivery.unitcost=vaccine.delivery.unitcost,
      vaccine.delivery.unitcost.baseline=vaccine.delivery.unitcost.baseline,
      unit.treat.cost.asymp=unit.treat.cost.asymp,
      unit.treat.cost.mild=unit.treat.cost.mild,
      unit.treat.cost.severe=unit.treat.cost.severe,
      unit.treat.cost.critical=unit.treat.cost.critical,
      cost.discount=cost.discount,
      effect.discount=effect.discount,
      DW.critical.disease =DW.critical.disease,
      DW.severe.disease  =DW.severe.disease,
      DW.mild.disease    = DW.mild.disease,
      DW.asymptomatic =DW.asymptomatic,
      dur.critical.covid = dur.critical.covid,
      dur.severe.covid =  dur.severe.covid,
      dur.mild.covid = dur.mild.covid,
      dur.asymptomatic.covid = dur.asymptomatic.covid,
      length.hosp.severe = length.hosp.severe,
      length.ICU.critical = length.ICU.critical
    ))
  }
  SobolSamples.list <- lapply(seq_len(nrow(SobolSamples)), function(i) SobolSamples[i,])  # convert each sample into a list
  inputs<- lapply(SobolSamples.list,FUN=create.inputs)
  return(inputs)
}


