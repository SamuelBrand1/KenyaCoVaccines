#if (!require("pacman")) install.packages("pacman")
#pkgs = c("randtoolbox",
#        "parallel",
#       "foreach",
#      "doParallel", "doSNOW")
#pacman::p_load(pkgs, character.only = T)
library(randtoolbox)
library(parallel)
library(foreach)
library(doParallel)
library(doSNOW)

rm(list = ls()) #remove previous variables before beginning the comparative analysis
ptm <- proc.time() #start clock
#set working directory
setwd(".../KenyaCoVaccinesPrivate")

#Add source files
source("Economic_Analysis/Sobol_inputs.R")# reads in the Sobol pobability sampling function
source("Economic_Analysis/Read_scenario_files.R") # reads in the clinical output files  
source("Economic_Analysis/CEA_function_calc.R") #does the cost effective analysis 

#-------------------------------------------------------------------------------------------
#read vaccine costs inputs from all the Sobol Samples
#inputs <- Create.Sobol.inputs(scenario=scenario,baseline.scenario=baseline.scenario,
 #                               discount.cost=T, discount.daly=T)
inputs <- Create.Sobol.inputs(7,6,0,1)
for (i in 1:length(inputs)) {
  inputs[[i]]$vaccine.delivery.unitcost <-inputs[[i]]$vaccine.delivery.unitcost[i]
  inputs[[i]]$vaccine.delivery.unitcost.baseline <-inputs[[i]]$vaccine.delivery.unitcost.baseline[i]
  inputs[[i]]$unit.treat.cost.asymp<-inputs[[i]]$unit.treat.cost.asymp[i]
  inputs[[i]]$unit.treat.cost.mild <-inputs[[i]]$unit.treat.cost.mild[i]
  inputs[[i]]$unit.treat.cost.severe <-inputs[[i]]$unit.treat.cost.severe[i]
  inputs[[i]]$unit.treat.cost.critical <-inputs[[i]]$unit.treat.cost.critical[i]
  inputs[[i]]$DW.critical.disease <-inputs[[i]]$DW.critical.disease[i]
  inputs[[i]]$DW.mild.disease <-inputs[[i]]$DW.mild.disease[i]
  inputs[[i]]$DW.severe.disease <-inputs[[i]]$DW.severe.disease[i]
  inputs[[i]]$length.hosp.severe <-inputs[[i]]$length.hosp.severe[i]
  inputs[[i]]$length.ICU.critical <-inputs[[i]]$length.ICU.critical[i]
}
#----------------------------------------------------------------------
#read in the clinical inputs for the comparative scenario (index, baseline) for all the MCMC samples  
#Case.files <- Read.Case.files(scenario=scenario,baseline.scenario=baseline.scenario)
Case.files <- Read.Case.files(7,6)
#----------------------------------------------------------------------
#run the CEA function for all MCMC*SOBOL outputs
cl<-makeCluster((parallel::detectCores()-1))
registerDoSNOW(cl)  # register cluster to enable parallel computing
clusterExport(cl,list = c("CEA"), envir = environment())
LEN =length(Case.files)
outputs= lapply(1:LEN, function(x){
  parLapply(cl,inputs, Case.files[[x]], fun=CEA)
} )
CEA_output <-matrix(unlist(outputs),ncol=18,byrow=TRUE)
colnames(CEA_output) <- c("Total.Cost","Total.Cost.baseline","Total.Daly",
                          "Total.DALY_baseline","Total.vacc.cost",
                          "Total.vacc.cost.baseline","Total.treatment.cost",
                          "Total.treatment.cost_baseline","Total.deaths",
                          "Total.deaths.baseline","Total.cases","Total.cases_baseline",
                          "Total.YLL","Total.YLL_baseline","Total.YLD",
                          "Total.YLD.baseline","Total.prod.costs",
                          "Total.prod.costs.baseline")

stopCluster(cl)   # stop the clusters
proc.time() - ptm # Stop the clock

#----------------------------------------------------------------------
#save the output of the specific compartive analysis
filename<- paste("scenario_",inputs[[1]]$scenario,"baseline_scenario_",inputs[[1]]$baseline.scenario,sep="")
save(CEA_output, file=paste("model_outputs/societal_persp/",paste("CEA_MODEL_OUTPUT10SOBOL5under",filename,sep="_"),sep="",".rda"))



# `   