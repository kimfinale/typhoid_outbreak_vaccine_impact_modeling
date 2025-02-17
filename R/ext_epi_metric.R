ext_epi_metric <- function(inc) {
  list(total_case = sum(inc$inc),
       child_case = sum(inc$inc[,1]),
       adult_case = sum(inc$inc[,2]),
       max_IR = max(rowSums(inc$inc)) / inc$param$population,
       epi_dur_day = sum(rowSums(inc$inc) > inc$param$case_threshold))
}
