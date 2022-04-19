## Gather data back from RDA files 
R"""
load(file = "generation_of_data_for_econ_analysis/waning_imm_scenario_1/kenya_inc_sev_no_vac_w_vac.rda")
"""
kenya_inc_sev_sc_1 = @rget kenya_inc_sev
kenya_inc_sev_sc_1 = sum(kenya_inc_sev_sc_1, dims = 2)[:, 1, :]
kenya_inc_sev = 0
R"""
kenya_inc_sev = 0 #Free memory
"""
R"""
load(file = "generation_of_data_for_econ_analysis/waning_imm_scenario_2/kenya_inc_sev_30_perc_w_vac.rda")
"""
kenya_inc_sev_sc_2 = @rget kenya_inc_sev
kenya_inc_sev_sc_2 = sum(kenya_inc_sev_sc_2, dims = 2)[:, 1, :]
kenya_inc_sev = 0
R"""
kenya_inc_sev = 0 #Free memory
"""
R"""
load(file = "generation_of_data_for_econ_analysis/waning_imm_scenario_3/kenya_inc_sev_50_perc_w_vac.rda")
"""
kenya_inc_sev_sc_3 = @rget kenya_inc_sev
kenya_inc_sev_sc_3 = sum(kenya_inc_sev_sc_3, dims = 2)[:, 1, :]
kenya_inc_sev = 0
R"""
kenya_inc_sev = 0 #Free memory
"""
R"""
load(file = "generation_of_data_for_econ_analysis/waning_imm_scenario_4/kenya_inc_sev_70_perc_w_vac.rda")
"""
kenya_inc_sev_sc_4 = @rget kenya_inc_sev
kenya_inc_sev_sc_4 = sum(kenya_inc_sev_sc_4, dims = 2)[:, 1, :]
kenya_inc_sev = 0
R"""
kenya_inc_sev = 0 #Free memory
"""
R"""
load(file = "generation_of_data_for_econ_analysis/waning_imm_scenario_5/kenya_inc_sev_30_perc_rapid_w_vac.rda")
"""
kenya_inc_sev_sc_5 = @rget kenya_inc_sev
kenya_inc_sev_sc_5 = sum(kenya_inc_sev_sc_5, dims = 2)[:, 1, :]
kenya_inc_sev = 0
R"""
kenya_inc_sev = 0 #Free memory
"""
R"""
load(file = "generation_of_data_for_econ_analysis/waning_imm_scenario_6/kenya_inc_sev_50_perc_rapid_w_vac.rda")
"""
kenya_inc_sev_sc_6 = @rget kenya_inc_sev
kenya_inc_sev_sc_6 = sum(kenya_inc_sev_sc_6, dims = 2)[:, 1, :]
kenya_inc_sev = 0
R"""
kenya_inc_sev = 0 #Free memory
"""
R"""
load(file = "generation_of_data_for_econ_analysis/waning_imm_scenario_7/kenya_inc_sev_70_perc_rapid_w_vac.rda")
"""
kenya_inc_sev_sc_7 = @rget kenya_inc_sev
kenya_inc_sev_sc_7 = sum(kenya_inc_sev_sc_7, dims = 2)[:, 1, :]
kenya_inc_sev = 0
R"""
kenya_inc_sev = 0 #Free memory
"""