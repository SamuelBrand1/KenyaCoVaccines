##Create callbacks for vaccine function based on 1st Dec 2020 start for sim

function start_vac1_30perc!(integrator)
    integrator.p[29] = fully_vaccinated_30perc[end] / (t_end - t_0)
end
function start_vac1_50perc!(integrator)
    integrator.p[29] = fully_vaccinated_50perc[end] / (t_end - t_0)
end
function start_vac1_70perc!(integrator)
    integrator.p[29] = fully_vaccinated_70perc[end] / (t_end - t_0)
end
function start_vac1_30perc_rapid!(integrator)
    integrator.p[29] = fully_vaccinated_30perc[end] / (t_end_rapid - t_0)
end
function start_vac1_50perc_rapid!(integrator)
    integrator.p[29] = fully_vaccinated_50perc[end] / (t_end_rapid - t_0)
end
function start_vac1_70perc_rapid!(integrator)
    integrator.p[29] = fully_vaccinated_70perc[end] / (t_end_rapid - t_0)
end

function start_vac2_30perc!(integrator)
    integrator.p[30] = fully_vaccinated_30perc[end] / (t_end - t_0)
end
function start_vac2_50perc!(integrator)
    integrator.p[30] = fully_vaccinated_50perc[end] / (t_end - t_0)
end
function start_vac2_70perc!(integrator)
    integrator.p[30] = fully_vaccinated_70perc[end] / (t_end - t_0)
end
function start_vac2_30perc_rapid!(integrator)
    integrator.p[30] = fully_vaccinated_30perc[end] / (t_end_rapid - t_0)
end
function start_vac2_50perc_rapid!(integrator)
    integrator.p[30] = fully_vaccinated_50perc[end] / (t_end_rapid - t_0)
end
function start_vac2_70perc_rapid!(integrator)
    integrator.p[30] = fully_vaccinated_70perc[end] / (t_end_rapid - t_0)
end

function stop_vac1!(integrator)
    integrator.p[29] = 0.0
end
function stop_vac2!(integrator)
    integrator.p[30] = 0.0
end

vac_cbs = [CallbackSet(),
    CallbackSet(
        PresetTimeCallback([(Date(2021, 9, 27) - Date(2020, 12, 1)).value - 56 + 14], start_vac1_30perc!, save_positions = (false, false)),
        PresetTimeCallback([(Date(2021, 9, 27) - Date(2020, 12, 1)).value + 14], start_vac2_30perc!, save_positions = (false, false)),
        PresetTimeCallback([(Date(2021, 9, 27) + Month(18) - Date(2020, 12, 1)).value - 56 + 14], stop_vac1!, save_positions = (false, false)),
        PresetTimeCallback([(Date(2021, 9, 27) + Month(18) - Date(2020, 12, 1)).value + 14], stop_vac2!, save_positions = (false, false))
    ),
    CallbackSet(
        PresetTimeCallback([(Date(2021, 9, 27) - Date(2020, 12, 1)).value - 56 + 14], start_vac1_50perc!, save_positions = (false, false)),
        PresetTimeCallback([(Date(2021, 9, 27) - Date(2020, 12, 1)).value + 14], start_vac2_50perc!, save_positions = (false, false)),
        PresetTimeCallback([(Date(2021, 9, 27) + Month(18) - Date(2020, 12, 1)).value - 56 + 14], stop_vac1!, save_positions = (false, false)),
        PresetTimeCallback([(Date(2021, 9, 27) + Month(18) - Date(2020, 12, 1)).value + 14], stop_vac2!, save_positions = (false, false))
    ),
    CallbackSet(
        PresetTimeCallback([(Date(2021, 9, 27) - Date(2020, 12, 1)).value - 56 + 14], start_vac1_70perc!, save_positions = (false, false)),
        PresetTimeCallback([(Date(2021, 9, 27) - Date(2020, 12, 1)).value + 14], start_vac2_70perc!, save_positions = (false, false)),
        PresetTimeCallback([(Date(2021, 9, 27) + Month(18) - Date(2020, 12, 1)).value - 56 + 14], stop_vac1!, save_positions = (false, false)),
        PresetTimeCallback([(Date(2021, 9, 27) + Month(18) - Date(2020, 12, 1)).value + 14], stop_vac2!, save_positions = (false, false))
    ),
    CallbackSet(
        PresetTimeCallback([(Date(2021, 9, 27) - Date(2020, 12, 1)).value - 56 + 14], start_vac1_30perc_rapid!, save_positions = (false, false)),
        PresetTimeCallback([(Date(2021, 9, 27) - Date(2020, 12, 1)).value + 14], start_vac2_30perc_rapid!, save_positions = (false, false)),
        PresetTimeCallback([(Date(2021, 9, 27) + Month(18) - Date(2020, 12, 1)).value - 56 + 14], stop_vac1!, save_positions = (false, false)),
        PresetTimeCallback([(Date(2021, 9, 27) + Month(18) - Date(2020, 12, 1)).value + 14], stop_vac2!, save_positions = (false, false))
    ),
    CallbackSet(
        PresetTimeCallback([(Date(2021, 9, 27) - Date(2020, 12, 1)).value - 56 + 14], start_vac1_50perc_rapid!, save_positions = (false, false)),
        PresetTimeCallback([(Date(2021, 9, 27) - Date(2020, 12, 1)).value + 14], start_vac2_50perc_rapid!, save_positions = (false, false)),
        PresetTimeCallback([(Date(2021, 9, 27) + Month(18) - Date(2020, 12, 1)).value - 56 + 14], stop_vac1!, save_positions = (false, false)),
        PresetTimeCallback([(Date(2021, 9, 27) + Month(18) - Date(2020, 12, 1)).value + 14], stop_vac2!, save_positions = (false, false))
    ),
    CallbackSet(
        PresetTimeCallback([(Date(2021, 9, 27) - Date(2020, 12, 1)).value - 56 + 14], start_vac1_70perc_rapid!, save_positions = (false, false)),
        PresetTimeCallback([(Date(2021, 9, 27) - Date(2020, 12, 1)).value + 14], start_vac2_70perc_rapid!, save_positions = (false, false)),
        PresetTimeCallback([(Date(2021, 9, 27) + Month(18) - Date(2020, 12, 1)).value - 56 + 14], stop_vac1!, save_positions = (false, false)),
        PresetTimeCallback([(Date(2021, 9, 27) + Month(18) - Date(2020, 12, 1)).value + 14], stop_vac2!, save_positions = (false, false))
    )
]