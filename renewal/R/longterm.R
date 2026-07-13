# =============================================================================
# Long-term (beyond-outbreak) TCV impact.
#
# A reactive campaign vaccinates a cohort of Npop; beyond averting the index
# outbreak, that cohort is protected against the region's ONGOING typhoid burden
# for a waning window. Two hazards contribute (both in TRUE-typhoid units):
#   (1) endemic background incidence  h_endemic = (inc_confirmed/1e5) / Se_BC
#   (2) recurrent future outbreaks    h_recur   = recur_per_year * mean_outbreak_true / Npop
# Protection wanes exponentially with mean duration tau; integrating to a horizon H,
# each vaccinee contributes  VE_eff * tau*(1-exp(-H/tau))  protected person-years,
# with VE_eff = coverage * VE0 * (1 + indirect).
#
# NB units: background incidence estimates (GBD, SEFI, surveillance) are CONFIRMED
# or modelled TRUE cases, so we divide CONFIRMED incidence by Se_BC to reach true
# typhoid, matching the outbreak-averted quantity (which is already true after the
# PPV pi-scaling). Set se_bc=1 if the incidence source is already true-typhoid.
# =============================================================================

# Protected person-years per vaccinee over horizon H (exponential waning, mean tau).
protected_py <- function(tau_years, horizon_years) tau_years * (1 - exp(-horizon_years / tau_years))

longterm_averted <- function(Npop, endemic_inc_conf_per100k = 0,
                             recur_per_year = 0, mean_outbreak_true = 0,
                             VE0 = 0.80, coverage = 0.85, indirect = 0,
                             tau_years = 15, horizon_years = 7, se_bc = 0.62) {
  h_endemic <- (endemic_inc_conf_per100k / 1e5) / se_bc          # true, per person-year
  h_recur   <- if (Npop > 0) recur_per_year * mean_outbreak_true / Npop else 0
  hazard    <- h_endemic + h_recur
  VE_eff    <- min(coverage * VE0 * (1 + indirect), 1)
  py        <- protected_py(tau_years, horizon_years)
  averted   <- Npop * hazard * VE_eff * py
  list(averted = averted,
       averted_endemic = Npop * h_endemic * VE_eff * py,
       averted_recur   = Npop * h_recur   * VE_eff * py,
       VE_eff = VE_eff, protected_py = py, hazard = hazard)
}
