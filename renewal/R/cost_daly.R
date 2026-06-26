# =============================================================================
# Cost / DALY integration.
#
# REUSES the manuscript's CEA functions verbatim (compute_yld, compute_yll,
# add_cea_results in R/2-functions.R). This module only (a) sets the globals
# those functions expect and (b) builds the country lookups they need
# (coi_per_patient, vaccine delivery cost, life expectancy, GDP), pointed at the
# files actually present in data/ — whose names differ from R/1-data.R's
# hardcoded ones. Logic mirrors R/1-data.R (CPI inflation, severity-weighted COI).
#
# If any required input file is missing, setup_cost_env() returns NULL and the
# caller skips the cost arm (the transmission outputs are unaffected).
# =============================================================================

# Files needed for the cost/DALY arm (beyond parameters.csv / country_costs.xlsx).
cost_files <- function() list(
  costs = "data/country_costs.xlsx",
  cpi   = "data/imf-dm-export-20241028.xls",
  life  = "data/wpp_life_expectancy_20241022.csv",
  gdp   = "data/GDP_WorldBank.xls"
)

cost_inputs_available <- function() all(vapply(cost_files(), file.exists, logical(1)))

# Build lookups + set globals used by add_cea_results(); returns lookups or NULL.
setup_cost_env <- function(cfg) {
  if (!cost_inputs_available()) {
    warning("Cost/DALY input file(s) missing; skipping cost arm.")
    return(NULL)
  }
  suppressMessages({ library(readxl); library(dplyr); library(tidyr) })
  f <- cost_files()
  ge <- globalenv()

  # --- parameters (mirror R/1-data.R global setup) ---------------------------
  parms <- as.data.frame(read.csv("data/parameters.csv", stringsAsFactors = FALSE))
  parms$Parameter <- trimws(parms$Parameter)   # CSV has trailing spaces on some names
  parms$Value <- suppressWarnings(as.numeric(parms$Value))
  parms <- parms[parms$Disease == "Typhoid", ]
  getp <- function(p) as.numeric(parms[parms$Parameter == p, "Value"][1])
  assign("parms", parms, envir = ge)
  assign("vacc_cost_per_dose", getp("Vaccine_Cost"), envir = ge)
  vacc_equipment_cost <- getp("Vaccine_Equipment_Cost")
  assign("dose_regimen", 1, envir = ge)
  avg_age_inf <- getp("Mean_Age_Infection"); assign("avg_age_inf", avg_age_inf, envir = ge)
  pr_mod <- getp("Prop_Moderate"); pr_sev <- getp("Prop_Severe")
  pr_gi  <- getp("Prop_GI_Bleed"); pr_ip  <- getp("Prop_Intestinal_Perf")
  pr_tot <- pr_mod + pr_sev + pr_gi + pr_ip

  # --- CPI inflation factor (mirror R/1-data.R:49-72) ------------------------
  cpi <- as.data.frame(read_xls(f$cpi))
  region_col <- names(cpi)[1]
  yrs <- as.character(2017:2025)
  wrow <- which(cpi[[region_col]] == "World")            # original deflator: World
  if (length(wrow) == 0) {                               # provided file may lack it
    avail <- cpi[[region_col]][!is.na(cpi[[region_col]]) & cpi[[region_col]] != ""]
    if (length(avail) == 0) stop("No usable CPI series in ", basename(f$cpi))
    warning(sprintf("No 'World' CPI row in %s; using '%s' as inflation proxy.",
                    basename(f$cpi), avail[1]))
    wrow <- which(cpi[[region_col]] == avail[1])
  }
  xi <- as.numeric(cpi[wrow[1], yrs]) / 100 + 1
  xi[is.na(xi)] <- 1                                      # missing year -> no inflation
  cumulative_inflation_factor  <- prod(xi)               # 2017..2025
  cumulative_inflation_factor2 <- prod(tail(xi, 3))      # 2023..2025
  assign("vacc_equipment_cost_adj", vacc_equipment_cost * cumulative_inflation_factor, envir = ge)

  # --- COI per patient by country (mirror R/1-data.R:66-85) ------------------
  costs <- as.data.frame(read_xlsx(f$costs))
  special <- c("DEU", "CHN", "USA", "MEX", "FJI")     # use 3-yr factor for these
  costs$Value_New <- ifelse(!(costs$ISO3 %in% special),
                            costs$Value * cumulative_inflation_factor,
                            costs$Value * cumulative_inflation_factor2)
  costs2 <- costs %>% distinct(ISO3, Estimate, Value_New) %>%
    pivot_wider(names_from = Estimate, values_from = Value_New)
  costs2$coi_per_patient <-
    (pr_mod / pr_tot * costs2$Cost_Outpatient) +
    (pr_sev / pr_tot * costs2$Cost_Inpatient) +
    (pr_gi  / pr_tot * costs2$Cost_Inpatient) +
    (pr_ip  / pr_tot * costs2$Cost_Inpatient)

  # --- Life expectancy at mean age of infection (mirror R/1-data.R:14,87) ----
  le <- read.csv(f$life, check.names = FALSE, stringsAsFactors = FALSE)
  le <- le[le$Variant == "Estimates", ]
  life_exp <- data.frame(
    country = le[["ISO3 Alpha-code"]],
    year    = suppressWarnings(as.integer(le[["Year"]])),
    life_exp = suppressWarnings(as.numeric(le[[as.character(avg_age_inf)]])),
    stringsAsFactors = FALSE)
  life_exp <- life_exp[!is.na(life_exp$country) & life_exp$country != "", ]

  # --- GDP per capita long (mirror R/1-data.R:97-102) ------------------------
  gdp <- as.data.frame(read_xls(f$gdp))
  names(gdp)[names(gdp) == "Country Code"] <- "country"
  yr_cols <- intersect(as.character(2001:2023), names(gdp))
  gdp_long <- gdp %>% select(c("country", all_of(yr_cols))) %>%
    pivot_longer(cols = all_of(yr_cols), names_to = "year", values_to = "gdp")
  gdp_long$year <- as.integer(gdp_long$year)
  gdp_long$gdp  <- suppressWarnings(as.numeric(gdp_long$gdp))

  list(costs2 = as.data.frame(costs2), life_exp = life_exp,
       gdp_long = as.data.frame(gdp_long), avg_age_inf = avg_age_inf)
}

# Attach cost inputs to the scenario df, then run the reused add_cea_results().
# Year is clamped to GDP coverage [2001, 2023] for the lookup only.
add_cost_daly <- function(df, lookups) {
  df$country <- trimws(df$country)               # defensive: drop any stray spaces
  yr <- pmin(pmax(df$year, 2001L), 2023L)
  m  <- match(df$country, trimws(lookups$costs2$ISO3))
  df$coi_per_patient       <- lookups$costs2$coi_per_patient[m]
  df$vaccine_delivery_cost <- lookups$costs2$Vaccine_Delivery[m]
  li <- match(paste(df$country, df$year), paste(lookups$life_exp$country, lookups$life_exp$year))
  df$life_exp <- as.numeric(lookups$life_exp$life_exp[li])
  gi <- match(paste(df$country, yr), paste(lookups$gdp_long$country, lookups$gdp_long$year))
  df$gdp <- as.numeric(lookups$gdp_long$gdp[gi])
  add_cea_results(df)   # <-- reused verbatim from R/2-functions.R
}
