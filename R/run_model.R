# run the model and extract output for the assigned time and variable
run_model <- function(model, params, saveat=1, output_var=NULL){
  out <- model(params)
  time_filter <- seq(1, by=ceiling(saveat/params$tau),
                     length.out=(ceiling(params$ndays/saveat)))

  if(is.null(output_var)) {
    output_var <- params$measure_var
  }
  out <- out[time_filter, output_var, drop=FALSE]
  out$time <- seq(0, length.out=nrow(out), by=saveat)

  return(out)
}
