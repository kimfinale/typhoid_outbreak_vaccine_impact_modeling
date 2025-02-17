#' Incidence
#' Runs the model to extract difference of the state variables across
#' the unit_time in days
#'
#' @param pars # parameter values to be varied
#' @param pars_baseline # This represents a complete set of parameters needed to run the model including initial values for the state variables. Some of the parameter values may be updated based on the pars
#' @param model # model (e.g., seirw)
#' @param output_time #
#' @param output_state
#'
#' @return
#' @export
#'
#' @examples
#'
incidence <- function(pars,
                      pars_baseline=NULL,
                      model=NULL){


  if (is.null(pars_baseline)) {
    # PARAMETERS is a global variable holding baseline parameters
    stop("Baseline parameters must be provided")
  }

  params = pars_baseline
  # updating happens through the use of names
  # names(pars) <- c("day1", "prop_immune", "R0", "prop_report")

  if (!is.null(pars)) params = update_params(pars, params)

  if (is.null(model)) {
    out <- params$model(params)
  }
  else {
    out <- model(params)
  }

  # extract by time_filter to match incidence across the
  # time_filter <- seq(1, by=round(output_time/params[["tau"]]),
  #                    length.out=(round(params[["ndays"]]/output_time)+1))

  # since the incidence during the initial params[["Day1"]] are not reported,
  # the time steps to be extracted are selected from the end. The time steps are
  # selected every round(output_time/params[["tau"]]) intervals. Then the
  # list of steps are then reversed such that values are extracted
  # in the correct order
  # time_filter <- rev(seq(nrow(out), by= - round(output_time/params[["tau"]]),
  #                    length.out=(round(params[["sobs_length"]]/output_time)+1)))
  # incidence is more useful when it returns the full length (i.e., ndays) of
  # the output. After that, the length could be adjusted
  time_filter <-
    rev(seq(nrow(out),
            by= - round(params[["output_days"]]/params[["tau"]]),
            length.out=ceiling(params[["ndays"]]/params[["output_days"]])))

  output_state = params[["output_state"]]
  ids = grep(output_state, names(out))
  output_state = names(out)[ids]
  # out <- out[time_filter, cols, drop=FALSE]
  out <- out[time_filter, output_state, drop=FALSE]
  # if measured variables is more than one
  df <- data.frame(matrix(NA, nrow=nrow(out)-1, ncol=length(output_state)))
  names(df) <- output_state
  for (i in 1:length(output_state)) {
    df[,i] <- diff(out[, output_state[i]])
  }

  return(list(inc=df, param=params))
}
