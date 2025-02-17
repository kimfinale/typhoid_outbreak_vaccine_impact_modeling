# extract start and end dates
ext_outbreak_st <- function (outbreak_id) { # st = space and time
  loc <- gsub("[0-9+]|-", "", outbreak_id)
  dates = gsub("[a-zA-Z]|::", "", outbreak_id)
  str = strsplit(dates, "-")
  start_d <- as.Date(paste0(str[[1]][2], "-", str[[1]][3], "-", str[[1]][4]))
  end_d <- as.Date(paste0(str[[1]][5], "-", str[[1]][6], "-", str[[1]][7]))

  return(list(location=loc, start_date=start_d, end_date=end_d))

}
