cost_func_pars = function(x, pars) {
  # simulate outbreak

  R0_obs = pars[["outbk_R0"]]
  pop_obs = pars[["popsize"]]
  threshold_obs = pars[["threshold"]]
  params$R0 = pars[["outbk_R0"]]
  size_obs = pars[["outbk_size"]]
  dur_obs = pars[["outbk_dur"]]
  timepeak_obs = pars[["outbk_time_to_peak"]]

  outbk = simulate_outbreak(params = params,
                            saveat = 7,
                            stat_window = 7,# weekly statistics!
                            popsize = pop_obs,
                            outbreak_threshold = threshold_obs,
                            prop_eff_pop=x)


  size_sim = outbk$outbreak_char[["outbreak_size"]]
  dur_sim = outbk$outbreak_char[["outbreak_dur"]]
  timepeak_sim = outbk$outbreak_char[["time_to_peak"]]
  # Least absolute relative errors (LARE)
  # https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3762514/
  abs((size_obs - size_sim)/size_obs) +
    abs((size_obs - size_sim)/size_sim)+
    abs((dur_obs - dur_sim)/dur_obs) +
    abs((dur_obs - dur_sim)/dur_sim)+
    abs((timepeak_obs - timepeak_sim)/timepeak_obs) +
    abs((timepeak_obs - timepeak_sim)/timepeak_sim)
}
