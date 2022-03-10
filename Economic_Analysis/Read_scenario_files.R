Read.Case.files <- function(scenario=scenario,baseline.scenario=baseline.scenario)
  {
  #---------------------------------------------------------------
  #load vaccination files-[changed this so as to have number of doses as a data frame]
  #----------------------------------------------------------------
  load(file = "generation_of_data_for_econ_analysis/vaccinations_per_day.rda")
  vaccination_rates <- vaccination_rates_kenya[-c(1:556,1162:1230),] #Leave 1.5years time frame (31Aug21-27Mar23)
  num.doses.procured.30 <- vaccination_rates[,-c(1,4:13)]
  num.doses.procured.30$tot_doses_30_perc <- num.doses.procured.30$number_first_dose_30_perc + num.doses.procured.30$number_second_dose_30_perc
  num.doses.procured_2 <- subset(num.doses.procured.30, select = tot_doses_30_perc )
  num.doses.procured.50 <- vaccination_rates[,-c(1:3,6:13)]
  num.doses.procured.50$tot_doses_50_perc <- num.doses.procured.50$number_first_dose_50_perc + num.doses.procured.50$number_second_dose_50_perc
  num.doses.procured_3 <- subset(num.doses.procured.50, select = tot_doses_50_perc )
  num.doses.procured.70 <- vaccination_rates[,-c(1:5,8:13)]
  num.doses.procured.70$tot_doses_70_perc <- num.doses.procured.70$number_first_dose_70_perc + num.doses.procured.70$number_second_dose_70_perc
  num.doses.procured_4 <- subset(num.doses.procured.70, select = tot_doses_70_perc )
  num.doses.procured.30rapid <- vaccination_rates[,-c(1:7,10:13)]
  num.doses.procured.30rapid$tot_doses_30_perc_rapid <- num.doses.procured.30rapid$number_first_dose_30_perc_rapid + num.doses.procured.30rapid$number_second_dose_30_perc_rapid
  num.doses.procured_5 <- subset(num.doses.procured.30rapid, select = tot_doses_30_perc_rapid )
  num.doses.procured.50rapid <- vaccination_rates[,-c(1:9,12:13)]
  num.doses.procured.50rapid$tot_doses_50_perc_rapid <- num.doses.procured.50rapid$number_first_dose_50_perc_rapid + num.doses.procured.50rapid$number_second_dose_50_perc_rapid
  num.doses.procured_6 <- subset(num.doses.procured.50rapid, select = tot_doses_50_perc_rapid )
  num.doses.procured.70rapid <- vaccination_rates[,-c(1:11)]
  num.doses.procured.70rapid$tot_doses_70_perc_rapid <- num.doses.procured.70rapid$number_first_dose_70_perc_rapid + num.doses.procured.70rapid$number_second_dose_70_perc_rapid
  num.doses.procured_7 <- subset(num.doses.procured.70rapid, select = tot_doses_70_perc_rapid )
  
  
  #----------------------------------------------------------------------  
  #load in and save the 2 comparative scenarios to be analysed 
  #the severe, critical and deaths are multiplied by 5 to assume
  #an under reporting of 5 in hoprital cases and in deaths
  #----------------------------------------------------------------------
  LEN_MCMC =2000 # number of MCMC samples used to forecast clinical outcomes per scenario
  #-----------------------------------------------------------------
  if (scenario == 2 && baseline.scenario == 1){
  #load scenario 1
  load(file = "generation_of_data_for_econ_analysis/immune_escape_scenario_1/kenya_inc_A_no_vac_imm_esc.rda")
  load(file = "generation_of_data_for_econ_analysis/immune_escape_scenario_1/kenya_inc_M_no_vac_imm_esc.rda")
  load(file = "generation_of_data_for_econ_analysis/immune_escape_scenario_1/kenya_inc_sev_no_vac_imm_esc.rda")
  load(file = "generation_of_data_for_econ_analysis/immune_escape_scenario_1/kenya_inc_crit_no_vac_imm_esc.rda")
  load(file = "generation_of_data_for_econ_analysis/immune_escape_scenario_1/kenya_inc_deaths_no_vac_imm_esc.rda")
  asymp_array_1 = kenya_inc_A[-c(1:273,848:912),,] # Leaves data for 1.1 years-31Aug21-27Mar23
  mild_array_1 = kenya_inc_M[-c(1:273,848:912),,] # Leaves data for 1.1 years-31Aug21-27Mar23
  severe_array_1 = kenya_inc_sev[-c(1:273,848:912),,]*5 # Leaves data for 1.1 years-31Aug21-27Mar23 and adjusts for under-reporting of hospitalizations by a factor of 1
  critical_array_1 = kenya_inc_crit[-c(1:273,848:912),,]*5 # Leaves data for 1.1 years-31Aug21-27Mar23 and adjusts for under-reporting of hospitalizations by a factor of 1
  deaths_array_1 = kenya_inc_deaths[-c(1:273,848:912),,]*5 # Leaves data for 1.1 years-31Aug21-27Mar23 and adjusts for under-reporting of deaths by a factor of 1
  doses_1 = 0
  rm(kenya_inc_A,kenya_inc_M,kenya_inc_sev,kenya_inc_crit,kenya_inc_deaths)  #delete previous array with same var name
  
  #load scenario 2
  load(file = "generation_of_data_for_econ_analysis/immune_escape_scenario_2/kenya_inc_A_30_perc_imm_esc.rda")
  load(file = "generation_of_data_for_econ_analysis/immune_escape_scenario_2/kenya_inc_M_30_perc_imm_esc.rda")
  load(file = "generation_of_data_for_econ_analysis/immune_escape_scenario_2/kenya_inc_sev_30_perc_imm_esc.rda")
  load(file = "generation_of_data_for_econ_analysis/immune_escape_scenario_2/kenya_inc_crit_30_perc_imm_esc.rda")
  load(file = "generation_of_data_for_econ_analysis/immune_escape_scenario_2/kenya_inc_deaths_30_perc_imm_esc.rda")
  asymp_array_2 = kenya_inc_A[-c(1:273,848:912),,] # Leaves data for 1.1 years-31Aug21-27Mar23
  mild_array_2 = kenya_inc_M[-c(1:273,848:912),,] # Leaves data for 1.1 years-31Aug21-27Mar23
  severe_array_2 = kenya_inc_sev[-c(1:273,848:912),,]*5# Leaves data for 1.1 years-31Aug21-27Mar23 and adjust for under-reporting of hospitalizations by a factor of 1
  critical_array_2 = kenya_inc_crit[-c(1:273,848:912),,]*5 # Leaves data for 1.1 years-31Aug21-27Mar23 and adjust for under-reporting of hospitalizations by a factor of 1
  deaths_array_2 = kenya_inc_deaths[-c(1:273,848:912),,]*5 # Leaves data for 1.1 years-31Aug21-27Mar23 and adjusts for under-reporting of deaths by a factor of 1 
  doses_2 = sum(num.doses.procured_2)
  rm(kenya_inc_A,kenya_inc_M,kenya_inc_sev,kenya_inc_crit,kenya_inc_deaths)  #delete previous array with same var name
  
  #save scenario(2) and baseline scenario(1) files in a looped list
  case.files_21 = list()
  case.files = list()
  for (i in 1:LEN_MCMC) {
      asymp_cases =   asymp_array_2[,,i] #1 
      mild_cases =   mild_array_2[,,i] #2
      severe_cases=   severe_array_2[,,i] #3
      critical_cases =   critical_array_2[,,i] #4
      deaths = deaths_array_2[,,i] #5
      num.doses.procured = doses_2 #6
      asymp_cases.baseline =   asymp_array_1[,,i] #7
      mild_cases.baseline =   mild_array_1[,,i] #8
      severe_cases.baseline =   severe_array_1[,,i] #9
      critical_cases.baseline =   critical_array_1[,,i] #10
      deaths.baseline= deaths_array_1[,,i] #11
      num.doses.procured.baseline = doses_1 #12
      case.files_21[[i]] = list(asymp_cases,mild_cases,severe_cases,critical_cases,
                               deaths,num.doses.procured,asymp_cases.baseline,
                               mild_cases.baseline,severe_cases.baseline,
                               critical_cases.baseline,deaths.baseline,
                               num.doses.procured.baseline)
      names(case.files_21[[i]]) <- c("asymp_cases", "mild_cases", "severe_cases",
                                    "critical_cases","deaths","num.doses.procured",
                                    "asymp_cases.baseline","mild_cases.baseline",
                                    "severe_cases.baseline","critical_cases.baseline",
                                    "deaths.baseline","num.doses.procured.baseline")
  }
  case.files <- case.files_21 
  }
#-------------------------------------------------------------------------
  else if (scenario == 3 && baseline.scenario == 2){
    #load scenario 2 files
    load(file = "generation_of_data_for_econ_analysis/immune_escape_scenario_2/kenya_inc_A_30_perc_imm_esc.rda")
    load(file = "generation_of_data_for_econ_analysis/immune_escape_scenario_2/kenya_inc_M_30_perc_imm_esc.rda")
    load(file = "generation_of_data_for_econ_analysis/immune_escape_scenario_2/kenya_inc_sev_30_perc_imm_esc.rda")
    load(file = "generation_of_data_for_econ_analysis/immune_escape_scenario_2/kenya_inc_crit_30_perc_imm_esc.rda")
    load(file = "generation_of_data_for_econ_analysis/immune_escape_scenario_2/kenya_inc_deaths_30_perc_imm_esc.rda")
    asymp_array_2 = kenya_inc_A[-c(1:273,848:912),,] # Leaves data for 1.1 years-31Aug21-27Mar23
    mild_array_2 = kenya_inc_M[-c(1:273,848:912),,] # Leaves data for 1.1 years-31Aug21-27Mar23
    severe_array_2 = kenya_inc_sev[-c(1:273,848:912),,]*5# Leaves data for 1.1 years-31Aug21-27Mar23 and adjust for under-reporting of hospitalizations by a factor of 1
    critical_array_2 = kenya_inc_crit[-c(1:273,848:912),,]*5 # Leaves data for 1.1 years-31Aug21-27Mar23 and adjust for under-reporting of hospitalizations by a factor of 1
    deaths_array_2 = kenya_inc_deaths[-c(1:273,848:912),,]*5 # Leaves data for 1.1 years-31Aug21-27Mar23 and adjusts for under-reporting of deaths by a factor of 1 
    doses_2 = sum(num.doses.procured_2)
    rm(kenya_inc_A,kenya_inc_M,kenya_inc_sev,kenya_inc_crit,kenya_inc_deaths)  #delete previous array with same var name
    
    #load scenario 3
    load(file = "generation_of_data_for_econ_analysis/immune_escape_scenario_3/kenya_inc_A_50_perc_imm_esc.rda")
    load(file = "generation_of_data_for_econ_analysis/immune_escape_scenario_3/kenya_inc_M_50_perc_imm_esc.rda")
    load(file = "generation_of_data_for_econ_analysis/immune_escape_scenario_3/kenya_inc_sev_50_perc_imm_esc.rda")
    load(file = "generation_of_data_for_econ_analysis/immune_escape_scenario_3/kenya_inc_crit_50_perc_imm_esc.rda")
    load(file = "generation_of_data_for_econ_analysis/immune_escape_scenario_3/kenya_inc_deaths_50_perc_imm_esc.rda")
    asymp_array_3 = kenya_inc_A[-c(1:273,848:912),,] # Leaves data for 1.1 years-31Aug21-27Mar23
    mild_array_3 = kenya_inc_M[-c(1:273,848:912),,] # Leaves data for 1.1 years-31Aug21-27Mar23
    severe_array_3 = kenya_inc_sev[-c(1:273,848:912),,]*5 # Leaves data for 1.5 years-31Aug21-27Mar23 and adjust for under-reporting of hospitalizations by a factor of 5
    critical_array_3 = kenya_inc_crit[-c(1:273,848:912),,]*5 # Leaves data for 1.5 years-31Aug21-27Mar23 and adjust for under-reporting of hospitalizations by a factor of 5
    deaths_array_3 = kenya_inc_deaths[-c(1:273,848:912),,]*5 # Leaves data for 1.5 years-31Aug21-27Mar23 and adjust for under-reporting of deaths by a factor of 5
    doses_3 = sum(num.doses.procured_3)
    rm(kenya_inc_A,kenya_inc_M,kenya_inc_sev,kenya_inc_crit,kenya_inc_deaths)  #delete previous array with same var name
    
  #save scenario(3) and baseline scenario(2) files in a looped list
    case.files_32 = list()
    case.files = list()
    for (i in 1:LEN_MCMC) {
        asymp_cases =   asymp_array_3[,,i] #1 
        mild_cases =   mild_array_3[,,i] #2
        severe_cases=   severe_array_3[,,i] #3
        critical_cases =   critical_array_3[,,i] #4
        deaths = deaths_array_3[,,i] #5
        num.doses.procured = doses_3 #6
        asymp_cases.baseline =   asymp_array_2[,,i] #7
        mild_cases.baseline =   mild_array_2[,,i] #8
        severe_cases.baseline =   severe_array_2[,,i] #9
        critical_cases.baseline =   critical_array_2[,,i] #10
        deaths.baseline= deaths_array_2[,,i] #11
        num.doses.procured.baseline = doses_2 #12
        case.files_32[[i]] = list(asymp_cases,mild_cases,severe_cases,critical_cases,
                                 deaths,num.doses.procured,asymp_cases.baseline,
                                 mild_cases.baseline,severe_cases.baseline,
                                 critical_cases.baseline,deaths.baseline,
                                 num.doses.procured.baseline)
        names(case.files_32[[i]]) <- c("asymp_cases", "mild_cases", "severe_cases",
                                      "critical_cases","deaths","num.doses.procured",
                                      "asymp_cases.baseline","mild_cases.baseline",
                                      "severe_cases.baseline","critical_cases.baseline",
                                      "deaths.baseline","num.doses.procured.baseline")
      }
    case.files <- case.files_32 
  }
#---------------------------------------------------------------------- else if (scenario == 4 && baseline.scenario == 3){
  else if (scenario == 4 && baseline.scenario == 3){
    #load scenario 3
    load(file = "generation_of_data_for_econ_analysis/immune_escape_scenario_3/kenya_inc_A_50_perc_imm_esc.rda")
    load(file = "generation_of_data_for_econ_analysis/immune_escape_scenario_3/kenya_inc_M_50_perc_imm_esc.rda")
    load(file = "generation_of_data_for_econ_analysis/immune_escape_scenario_3/kenya_inc_sev_50_perc_imm_esc.rda")
    load(file = "generation_of_data_for_econ_analysis/immune_escape_scenario_3/kenya_inc_crit_50_perc_imm_esc.rda")
    load(file = "generation_of_data_for_econ_analysis/immune_escape_scenario_3/kenya_inc_deaths_50_perc_imm_esc.rda")
    asymp_array_3 = kenya_inc_A[-c(1:273,848:912),,] # Leaves data for 1.1 years-31Aug21-27Mar23
    mild_array_3 = kenya_inc_M[-c(1:273,848:912),,] # Leaves data for 1.1 years-31Aug21-27Mar23
    severe_array_3 = kenya_inc_sev[-c(1:273,848:912),,]*5 # Leaves data for 1.5 years-31Aug21-27Mar23 and adjust for under-reporting of hospitalizations by a factor of 5
    critical_array_3 = kenya_inc_crit[-c(1:273,848:912),,]*5 # Leaves data for 1.5 years-31Aug21-27Mar23 and adjust for under-reporting of hospitalizations by a factor of 5
    deaths_array_3 = kenya_inc_deaths[-c(1:273,848:912),,]*5 # Leaves data for 1.5 years-31Aug21-27Mar23 and adjust for under-reporting of deaths by a factor of 5
    doses_3 = sum(num.doses.procured_3)
    rm(kenya_inc_A,kenya_inc_M,kenya_inc_sev,kenya_inc_crit,kenya_inc_deaths)  #delete previous array with same var name
    
    #load scenario 4
    load(file = "generation_of_data_for_econ_analysis/immune_escape_scenario_4/kenya_inc_A_70_perc_imm_esc.rda")
    load(file = "generation_of_data_for_econ_analysis/immune_escape_scenario_4/kenya_inc_M_70_perc_imm_esc.rda")
    load(file = "generation_of_data_for_econ_analysis/immune_escape_scenario_4/kenya_inc_sev_70_perc_imm_esc.rda")
    load(file = "generation_of_data_for_econ_analysis/immune_escape_scenario_4/kenya_inc_crit_70_perc_imm_esc.rda")
    load(file = "generation_of_data_for_econ_analysis/immune_escape_scenario_4/kenya_inc_deaths_70_perc_imm_esc.rda")
    asymp_array_4 = kenya_inc_A[-c(1:273,848:912),,] # Leaves data for 1.5 years-31Aug21-27Mar23
    mild_array_4 = kenya_inc_M[-c(1:273,848:912),,] # Leaves data for 1.5 years-31Aug21-27Mar23
    severe_array_4 = kenya_inc_sev[-c(1:273,848:912),,]*5 # Leaves data for 1.5 years-31Aug21-27Mar23 and adjust for under-reporting of hospitalizations by a factor of 5
    critical_array_4 = kenya_inc_crit[-c(1:273,848:912),,]*5 # Leaves data for 1.5 years-31Aug21-27Mar23 and adjust for under-reporting of hospitalizations by a factor of 5
    deaths_array_4 = kenya_inc_deaths[-c(1:273,848:912),,]*5 # Leaves data for 1.5 years-31Aug21-27Mar23 and adjust for under-reporting of deaths by a factor of 5
    doses_4 = sum(num.doses.procured_4)
    rm(kenya_inc_A,kenya_inc_M,kenya_inc_sev,kenya_inc_crit,kenya_inc_deaths)  #delete previous array with same var name
    
    #save scenario(4) and baseline scenario(3) files in a looped list
  case.files_43 = list()
  case.files = list()
  for (i in 1:LEN_MCMC) {
      asymp_cases =   asymp_array_4[,,i] #1 
      mild_cases =   mild_array_4[,,i] #2
      severe_cases=   severe_array_4[,,i] #3
      critical_cases =   critical_array_4[,,i] #4
      deaths = deaths_array_4[,,i] #5
      num.doses.procured = doses_4 #6
      asymp_cases.baseline =   asymp_array_3[,,i] #7
      mild_cases.baseline =   mild_array_3[,,i] #8
      severe_cases.baseline =   severe_array_3[,,i] #9
      critical_cases.baseline =   critical_array_3[,,i] #10
      deaths.baseline= deaths_array_3[,,i] #11
      num.doses.procured.baseline = doses_3 #12
      case.files_43[[i]] = list(asymp_cases,mild_cases,severe_cases,critical_cases,
                               deaths,num.doses.procured,asymp_cases.baseline,
                               mild_cases.baseline,severe_cases.baseline,
                               critical_cases.baseline,deaths.baseline,
                               num.doses.procured.baseline)
      names(case.files_43[[i]]) <- c("asymp_cases", "mild_cases", "severe_cases",
                                    "critical_cases","deaths","num.doses.procured",
                                    "asymp_cases.baseline","mild_cases.baseline",
                                    "severe_cases.baseline","critical_cases.baseline",
                                    "deaths.baseline","num.doses.procured.baseline")
  }
  case.files <- case.files_43 
  } 
  #-------------------------------------------------------------------------
  else if (scenario == 5 && baseline.scenario == 1){
    #load scenario 1
    load(file = "generation_of_data_for_econ_analysis/immune_escape_scenario_1/kenya_inc_A_no_vac_imm_esc.rda")
    load(file = "generation_of_data_for_econ_analysis/immune_escape_scenario_1/kenya_inc_M_no_vac_imm_esc.rda")
    load(file = "generation_of_data_for_econ_analysis/immune_escape_scenario_1/kenya_inc_sev_no_vac_imm_esc.rda")
    load(file = "generation_of_data_for_econ_analysis/immune_escape_scenario_1/kenya_inc_crit_no_vac_imm_esc.rda")
    load(file = "generation_of_data_for_econ_analysis/immune_escape_scenario_1/kenya_inc_deaths_no_vac_imm_esc.rda")
    asymp_array_1 = kenya_inc_A[-c(1:273,848:912),,] # Leaves data for 1.1 years-31Aug21-27Mar23
    mild_array_1 = kenya_inc_M[-c(1:273,848:912),,] # Leaves data for 1.1 years-31Aug21-27Mar23
    severe_array_1 = kenya_inc_sev[-c(1:273,848:912),,]*5 # Leaves data for 1.1 years-31Aug21-27Mar23 and adjusts for under-reporting of hospitalizations by a factor of 1
    critical_array_1 = kenya_inc_crit[-c(1:273,848:912),,]*5 # Leaves data for 1.1 years-31Aug21-27Mar23 and adjusts for under-reporting of hospitalizations by a factor of 1
    deaths_array_1 = kenya_inc_deaths[-c(1:273,848:912),,]*5 # Leaves data for 1.1 years-31Aug21-27Mar23 and adjusts for under-reporting of deaths by a factor of 1
    doses_1 = 0
    rm(kenya_inc_A,kenya_inc_M,kenya_inc_sev,kenya_inc_crit,kenya_inc_deaths)  #delete previous array with same var name
    
    #load scenario 5
    load(file = "generation_of_data_for_econ_analysis/immune_escape_scenario_5/kenya_inc_A_30_perc_rapid_imm_esc.rda")
    load(file = "generation_of_data_for_econ_analysis/immune_escape_scenario_5/kenya_inc_M_30_perc_rapid_imm_esc.rda")
    load(file = "generation_of_data_for_econ_analysis/immune_escape_scenario_5/kenya_inc_sev_30_perc_rapid_imm_esc.rda")
    load(file = "generation_of_data_for_econ_analysis/immune_escape_scenario_5/kenya_inc_crit_30_perc_rapid_imm_esc.rda")
    load(file = "generation_of_data_for_econ_analysis/immune_escape_scenario_5/kenya_inc_deaths_30_perc_rapid_imm_esc.rda")
    asymp_array_5 = kenya_inc_A[-c(1:273,848:912),,] # Leaves data for 1.5 years-31Aug21-27Mar23
    mild_array_5 = kenya_inc_M[-c(1:273,848:912),,] # Leaves data for 1.5 years-31Aug21-27Mar23
    severe_array_5 = kenya_inc_sev[-c(1:273,848:912),,]*5 # Leaves data for 1.5 years-31Aug21-27Mar23 and adjust for under-reporting of hospitalizations by a factor of 5
    critical_array_5 = kenya_inc_crit[-c(1:273,848:912),,]*5 # Leaves data for 1.5 years-31Aug21-27Mar23 and adjust for under-reporting of hospitalizations by a factor of 5
    deaths_array_5 = kenya_inc_deaths[-c(1:273,848:912),,]*5 # Leaves data for 1.5 years-31Aug21-27Mar23 and adjust for under-reporting of deaths by a factor of 5
    doses_5 = sum(num.doses.procured_5)
    rm(kenya_inc_A,kenya_inc_M,kenya_inc_sev,kenya_inc_crit,kenya_inc_deaths)  #delete previous array with same var name
    
    #save scenario(5) and baseline scenario(1) files in a looped list
    case.files_51 = list()
    case.files = list()
    for (i in 1:LEN_MCMC) {
      asymp_cases =   asymp_array_5[,,i] #1 
      mild_cases =   mild_array_5[,,i] #2
      severe_cases=   severe_array_5[,,i] #3
      critical_cases =   critical_array_5[,,i] #4
      deaths = deaths_array_5[,,i] #5
      num.doses.procured = doses_5 #6
      asymp_cases.baseline =   asymp_array_1[,,i] #7
      mild_cases.baseline =   mild_array_1[,,i] #8
      severe_cases.baseline =   severe_array_1[,,i] #9
      critical_cases.baseline =   critical_array_1[,,i] #10
      deaths.baseline= deaths_array_1[,,i] #11
      num.doses.procured.baseline = doses_1 #12
      case.files_51[[i]] = list(asymp_cases,mild_cases,severe_cases,critical_cases,
                               deaths,num.doses.procured,asymp_cases.baseline,
                               mild_cases.baseline,severe_cases.baseline,
                               critical_cases.baseline,deaths.baseline,
                               num.doses.procured.baseline)
      names(case.files_51[[i]]) <- c("asymp_cases", "mild_cases", "severe_cases",
                                    "critical_cases","deaths","num.doses.procured",
                                    "asymp_cases.baseline","mild_cases.baseline",
                                    "severe_cases.baseline","critical_cases.baseline",
                                    "deaths.baseline","num.doses.procured.baseline")
  }
  case.files <- case.files_51 
  } 
  #-------------------------------------------------------------------------
  else if (scenario == 6 && baseline.scenario == 5){
  #load scenario 5
  load(file = "generation_of_data_for_econ_analysis/immune_escape_scenario_5/kenya_inc_A_30_perc_rapid_imm_esc.rda")
  load(file = "generation_of_data_for_econ_analysis/immune_escape_scenario_5/kenya_inc_M_30_perc_rapid_imm_esc.rda")
  load(file = "generation_of_data_for_econ_analysis/immune_escape_scenario_5/kenya_inc_sev_30_perc_rapid_imm_esc.rda")
  load(file = "generation_of_data_for_econ_analysis/immune_escape_scenario_5/kenya_inc_crit_30_perc_rapid_imm_esc.rda")
  load(file = "generation_of_data_for_econ_analysis/immune_escape_scenario_5/kenya_inc_deaths_30_perc_rapid_imm_esc.rda")
  asymp_array_5 = kenya_inc_A[-c(1:273,848:912),,] # Leaves data for 1.5 years-31Aug21-27Mar23
  mild_array_5 = kenya_inc_M[-c(1:273,848:912),,] # Leaves data for 1.5 years-31Aug21-27Mar23
  severe_array_5 = kenya_inc_sev[-c(1:273,848:912),,]*5 # Leaves data for 1.5 years-31Aug21-27Mar23 and adjust for under-reporting of hospitalizations by a factor of 5
  critical_array_5 = kenya_inc_crit[-c(1:273,848:912),,]*5 # Leaves data for 1.5 years-31Aug21-27Mar23 and adjust for under-reporting of hospitalizations by a factor of 5
  deaths_array_5 = kenya_inc_deaths[-c(1:273,848:912),,]*5 # Leaves data for 1.5 years-31Aug21-27Mar23 and adjust for under-reporting of deaths by a factor of 5
  doses_5 = sum(num.doses.procured_5)
  rm(kenya_inc_A,kenya_inc_M,kenya_inc_sev,kenya_inc_crit,kenya_inc_deaths)  #delete previous array with same var name
  
  #load scenario 6
  #rm(kenya_inc_A,kenya_inc_M,kenya_inc_sev,kenya_inc_crit,kenya_inc_deaths)  #delete previous array with same var name
  load(file = "generation_of_data_for_econ_analysis/immune_escape_scenario_6/kenya_inc_A_50_perc_rapid_imm_esc.rda")
  load(file = "generation_of_data_for_econ_analysis/immune_escape_scenario_6/kenya_inc_M_50_perc_rapid_imm_esc.rda")
  load(file = "generation_of_data_for_econ_analysis/immune_escape_scenario_6/kenya_inc_sev_50_perc_rapid_imm_esc.rda")
  load(file = "generation_of_data_for_econ_analysis/immune_escape_scenario_6/kenya_inc_crit_50_perc_rapid_imm_esc.rda")
  load(file = "generation_of_data_for_econ_analysis/immune_escape_scenario_6/kenya_inc_deaths_50_perc_rapid_imm_esc.rda")
  asymp_array_6 = kenya_inc_A[-c(1:273,848:912),,] # Leaves data for 1.5 years-31Aug21-27Mar23
  mild_array_6 = kenya_inc_M[-c(1:273,848:912),,] # Leaves data for 1.5 years-31Aug21-27Mar23
  severe_array_6 = kenya_inc_sev[-c(1:273,848:912),,]*5# Leaves data for 1.5 years-31Aug21-27Mar23 and adjust for under-reporting of hospitalizations by a factor of 5
  critical_array_6 = kenya_inc_crit[-c(1:273,848:912),,]*5 # Leaves data for 1.5 years-31Aug21-27Mar23 and adjust for under-reporting of hospitalizations by a factor of 5
  deaths_array_6 = kenya_inc_deaths[-c(1:273,848:912),,]*5 # Leaves data for 1.5 years-31Aug21-27Mar23 and adjust for under-reporting of deaths by a factor of 5
  doses_6 = sum(num.doses.procured_6)
  rm(kenya_inc_A,kenya_inc_M,kenya_inc_sev,kenya_inc_crit,kenya_inc_deaths)  #delete previous array with same var name
  
  #save scenario(6) and baseline scenario(5) files in a looped list
  case.files_65 = list()
  case.files = list()
  for (i in 1:LEN_MCMC) {
      asymp_cases =   asymp_array_6[,,i] #1 
      mild_cases =   mild_array_6[,,i] #2
      severe_cases=   severe_array_6[,,i] #3
      critical_cases =   critical_array_6[,,i] #4
      deaths = deaths_array_6[,,i] #5
      num.doses.procured = doses_6 #6
      asymp_cases.baseline =   asymp_array_5[,,i] #7
      mild_cases.baseline =   mild_array_5[,,i] #8
      severe_cases.baseline =   severe_array_5[,,i] #9
      critical_cases.baseline =   critical_array_5[,,i] #10
      deaths.baseline= deaths_array_5[,,i] #11
      num.doses.procured.baseline = doses_5 #12
      case.files_65[[i]] = list(asymp_cases,mild_cases,severe_cases,critical_cases,
                               deaths,num.doses.procured,asymp_cases.baseline,
                               mild_cases.baseline,severe_cases.baseline,
                               critical_cases.baseline,deaths.baseline,
                               num.doses.procured.baseline)
      names(case.files_65[[i]]) <- c("asymp_cases", "mild_cases", "severe_cases",
                                    "critical_cases","deaths","num.doses.procured",
                                    "asymp_cases.baseline","mild_cases.baseline",
                                    "severe_cases.baseline","critical_cases.baseline",
                                    "deaths.baseline","num.doses.procured.baseline")
    }
  case.files <- case.files_65 
  } 
  #-------------------------------------------------------------------------
  else if (scenario == 7 && baseline.scenario == 6){
  #load scenario 6
  #rm(kenya_inc_A,kenya_inc_M,kenya_inc_sev,kenya_inc_crit,kenya_inc_deaths)  #delete previous array with same var name
  load(file = "generation_of_data_for_econ_analysis/immune_escape_scenario_6/kenya_inc_A_50_perc_rapid_imm_esc.rda")
  load(file = "generation_of_data_for_econ_analysis/immune_escape_scenario_6/kenya_inc_M_50_perc_rapid_imm_esc.rda")
  load(file = "generation_of_data_for_econ_analysis/immune_escape_scenario_6/kenya_inc_sev_50_perc_rapid_imm_esc.rda")
  load(file = "generation_of_data_for_econ_analysis/immune_escape_scenario_6/kenya_inc_crit_50_perc_rapid_imm_esc.rda")
  load(file = "generation_of_data_for_econ_analysis/immune_escape_scenario_6/kenya_inc_deaths_50_perc_rapid_imm_esc.rda")
  asymp_array_6 = kenya_inc_A[-c(1:273,848:912),,] # Leaves data for 1.5 years-31Aug21-27Mar23
  mild_array_6 = kenya_inc_M[-c(1:273,848:912),,] # Leaves data for 1.5 years-31Aug21-27Mar23
  severe_array_6 = kenya_inc_sev[-c(1:273,848:912),,]*5# Leaves data for 1.5 years-31Aug21-27Mar23 and adjust for under-reporting of hospitalizations by a factor of 5
  critical_array_6 = kenya_inc_crit[-c(1:273,848:912),,]*5 # Leaves data for 1.5 years-31Aug21-27Mar23 and adjust for under-reporting of hospitalizations by a factor of 5
  deaths_array_6 = kenya_inc_deaths[-c(1:273,848:912),,]*5 # Leaves data for 1.5 years-31Aug21-27Mar23 and adjust for under-reporting of deaths by a factor of 5
  doses_6 = sum(num.doses.procured_6)
  rm(kenya_inc_A,kenya_inc_M,kenya_inc_sev,kenya_inc_crit,kenya_inc_deaths)  #delete previous array with same var name
  
  #load scenario 7
  load(file = "generation_of_data_for_econ_analysis/immune_escape_scenario_7/kenya_inc_A_70_perc_rapid_imm_esc.rda")
  load(file = "generation_of_data_for_econ_analysis/immune_escape_scenario_7/kenya_inc_M_70_perc_rapid_imm_esc.rda")
  load(file = "generation_of_data_for_econ_analysis/immune_escape_scenario_7/kenya_inc_sev_70_perc_rapid_imm_esc.rda")
  load(file = "generation_of_data_for_econ_analysis/immune_escape_scenario_7/kenya_inc_crit_70_perc_rapid_imm_esc.rda")
  load(file = "generation_of_data_for_econ_analysis/immune_escape_scenario_7/kenya_inc_deaths_70_perc_rapid_imm_esc.rda")
  load(file = "generation_of_data_for_econ_analysis/vaccinations_per_day.rda")
  asymp_array_7 = kenya_inc_A[-c(1:273,848:912),,] # Leaves data for 1.5 years-31Aug21-27Mar23
  mild_array_7 = kenya_inc_M[-c(1:273,848:912),,] # Leaves data for 1.5 years-31Aug21-27Mar23
  severe_array_7 = kenya_inc_sev[-c(1:273,848:912),,]*5# Leaves data for 1.5 years-31Aug21-27Mar23 and adjust for under-reporting of hospitalizations by a factor of 5
  critical_array_7 = kenya_inc_crit[-c(1:273,848:912),,]*5 # Leaves data for 1.5 years-31Aug21-27Mar23 and adjust for under-reporting of hospitalizations by a factor of 5
  deaths_array_7 = kenya_inc_deaths[-c(1:273,848:912),,]*5 # Leaves data for 1.5 years-31Aug21-27Mar23 and adjust for under-reporting of deaths by a factor of 5
  doses_7 = sum(num.doses.procured_7)
  rm(kenya_inc_A,kenya_inc_M,kenya_inc_sev,kenya_inc_crit,kenya_inc_deaths)  #delete previous array with same var name
  
  #save scenario(7) and baseline scenario(6) files in a looped list
  case.files_76 = list()
  case.files = list()
  for (i in 1:LEN_MCMC) {
      asymp_cases =   asymp_array_7[,,i] #1 
      mild_cases =   mild_array_7[,,i] #2
      severe_cases=   severe_array_7[,,i] #3
      critical_cases =   critical_array_7[,,i] #4
      deaths = deaths_array_7[,,i] #5
      num.doses.procured = doses_7 #6
      asymp_cases.baseline =   asymp_array_6[,,i] #7
      mild_cases.baseline =   mild_array_6[,,i] #8
      severe_cases.baseline =   severe_array_6[,,i] #9
      critical_cases.baseline =   critical_array_6[,,i] #10
      deaths.baseline= deaths_array_6[,,i] #11
      num.doses.procured.baseline = doses_6 #12
      case.files_76[[i]] = list(asymp_cases,mild_cases,severe_cases,critical_cases,
                               deaths,num.doses.procured,asymp_cases.baseline,
                               mild_cases.baseline,severe_cases.baseline,
                               critical_cases.baseline,deaths.baseline,
                               num.doses.procured.baseline)
      names(case.files_76[[i]]) <- c("asymp_cases", "mild_cases", "severe_cases",
                                    "critical_cases","deaths","num.doses.procured",
                                    "asymp_cases.baseline","mild_cases.baseline",
                                    "severe_cases.baseline","critical_cases.baseline",
                                    "deaths.baseline","num.doses.procured.baseline")
    }
  case.files <- case.files_76
}
return(case.files)
}