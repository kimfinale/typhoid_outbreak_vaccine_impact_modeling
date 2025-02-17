extract_outbreak_char <- function(out,
                                  stat_window=7,
                                  outbreak_threshold=1,
                                  popsize=1e5,
                                  inc_ref_pop=1e3){

  out$time_scaled = out$time / stat_window # stat_window (usu. a week)
  out = out[out$time %% stat_window == 0, ]
  incid_scaled = diff(out$ci) # incidence per stat_window
  over = which(incid_scaled > outbreak_threshold) # index for the week of which the incid is over the threshold
  start = NA
  end = NA
  if(length(over) > 1){
    start = out$time_scaled[min(over)]
    end = out$time_scaled[max(over)]
  }
  peak_index = which.max(incid_scaled)
  peak_time = out$time_scaled[peak_index]
  peak_inc_per_1e3 = incid_scaled[peak_index] / popsize * inc_ref_pop
  time_to_peak = (peak_time - start)
  outbreak_dur = (end - start)
  outbreak_size = tail(out$ci, 1)
  attack_rate_per_1e3 = outbreak_size / popsize * inc_ref_pop
  prop_peak_inc = incid_scaled[peak_index] / outbreak_size

  return(outbreak_char = c(outbreak_dur = outbreak_dur,
                    time_to_peak = time_to_peak,
                    peak_inc_per_1e3 = peak_inc_per_1e3,
                    outbreak_size = outbreak_size,
                    attack_rate_per_1e3 = attack_rate_per_1e3,
                    prop_peak_inc = prop_peak_inc))
}
