# run the model for the given parameter set and extract the characteristics
# of the resulting outbreak
simulate_outbreak = function(params=NULL,
                        model=seiarw_2ag_erlang,
                        popsize=1e5,
                        inc_ref_pop=1e3,
                        prop_s=1,
                        prop_i=0.001,
                        prop_eff_pop=1,
                        output_var=c("CI.1","CI.2"),
                        saveat=1,
                        stat_window=7,
                        outbreak_threshold=1){

  if (is.null(params)) {
    params = initialize_params()
  }

  effpop = popsize * prop_eff_pop

  params$init$R = (1-prop_s) * effpop
  params$init$I = prop_i * effpop
  params$init$S = effpop - params$init$I - params$init$R

  out <- run_model(model, params,
                   saveat = saveat,
                   output_var = output_var)

  out$ci = out$CI.1 + out$CI.2
  char = extract_outbreak_char(out=out,
                               stat_window=stat_window,
                               outbreak_threshold=outbreak_threshold,
                               popsize=popsize,
                               inc_ref_pop=inc_ref_pop)

  return(list(timeseries = out,
              outbreak_char = char))
}
  # out$time_scaled = out$time / output_time # output_time (usu. a week)
  #
  # incid_scaled = diff(out$ci) # incidence per output_time
  # over = which(incid_scaled > outbreak_threshold) # index for the week of which the incid is over the threshold
  # start = NA
  # end = NA
  # if(length(over) > 1){
  #   start = out$time_scaled[min(over)]
  #   end = out$time_scaled[max(over)]
  # }
  # peak_index = which.max(incid_scaled)
  # peak_time = out$time_scaled[peak_index]
  # peak_inc_per_1e3 = incid_scaled[peak_index] / popsize * inc_ref_pop
  # time_to_peak = (peak_time - start)
  # outbreak_dur = (end - start)
  # outbreak_size = tail(out$ci, 1)
  # attack_rate_per_1e3 = outbreak_size / popsize * inc_ref_pop
  # prop_peak_inc = incid_scaled[peak_index] / outbreak_size
  #
  # return(list(timeseries = out,
  #             outbreak_char = c(outbreak_dur = outbreak_dur,
  #               time_to_peak = time_to_peak,
  #               peak_inc_per_1e3 = peak_inc_per_1e3,
  #               outbreak_size = outbreak_size,
  #               attack_rate_per_1e3 = attack_rate_per_1e3,
  #               prop_peak_inc = prop_peak_inc)))

# }
