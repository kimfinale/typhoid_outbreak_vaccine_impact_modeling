# run the model for the given parameter set and extract the characteristics
# of the resulting outbreak
# simulate outbreak function created and this will be deleted?
run_outbreak = function(params,
                        model=seiarw_2ag_erlang,
                        popsize=1e5,
                        inc_ref_pop=1e3,
                        prop_s=1,
                        prop_i=0.001,
                        prop_eff_pop=1,
                        output_var=c("CI.1","CI.2"),
                        output_time=7,
                        threshold=1){

  popsize_scaled = popsize * prop_eff_pop

  params$init$R = (1-prop_s) * popsize_scaled
  params$init$I = prop_i * popsize_scaled
  params$init$S = popsize_scaled - params$init$I - params$init$R

  out <- run_model(model, params,
                   output_time=output_time,
                   output_var=output_var)

  df = data.frame(ci = out$CI.1 + out$CI.2)
  df$time = out$time / output_time # make it a week

  weekly = diff(df$ci)
  over = which(weekly > threshold) # index for the week of which the incid is over the threshold
  start = df$time[min(over)]
  end = df$time[max(over)]
  peak = df$time[which.max(weekly)]
  peak_inc_per_1e3 = weekly[peak] / popsize * inc_ref_pop
  time_to_peak = (peak - start)
  outbreak_dur = (end - start)
  outbreak_size = tail(df$ci, 1)
  attack_rate_per_1e3 = outbreak_size / popsize * inc_ref_pop

  return(c(outbreak_dur=outbreak_dur,
           time_to_peak=time_to_peak,
           peak_inc_per_1e3=peak_inc_per_1e3,
           outbreak_size=outbreak_size,
           attack_rate_per_1e3=attack_rate_per_1e3))
}
