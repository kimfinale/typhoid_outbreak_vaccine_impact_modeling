# =============================================================================
# Does covariate-assisted PPV assignment for UNANCHORED outbreaks change the
# burden? Compares two assignment rules at production ORI timing:
#   (a) "pooled"     - unanchored outbreaks get the fitted population PPV
#   (b) "strictness" - unanchored outbreaks get logit(pi_pop) + gamma*z(strictness),
#                      EXCEPT those whose strictness comes from the severity/epi-link
#                      gate beyond anchor support (gate > max anchor gate), which are
#                      held at the pooled value (unsupported extrapolation).
# Anchored outbreaks keep their own fitted PPV in both rules.
#   Rscript revision/run_ppv_assignment_comparison.R
# =============================================================================
setwd(Sys.getenv("RENEWAL_ROOT", "."))
source("revision/allocation_scenarios.R")
dir.create("revision/tables", showWarnings=FALSE, recursive=TRUE)
outdir <- "revision/main_text_figures"; dir.create(outdir, showWarnings=FALSE, recursive=TRUE)

cfg   <- yaml::read_yaml("renewal/config.yml")
prep  <- prep_outbreaks(cfg); w <- gi_from_config(cfg)
coverage <- cfg$vaccine$coverage_base; psi_T <- 0.8*cfg$vaccine$psi_mean
delta <- round((cfg$timing$campaign_duration_days+cfg$timing$immunity_onset_days)/cfg$step_days)
T_EFF <- cfg$timing$delay_base_weeks + delta                    # production timing

pic  <- read.csv("latent_class_ppv/tables/final_pi_community.csv", check.names=FALSE, stringsAsFactors=FALSE)
pars <- read.csv("latent_class_ppv/tables/final_parameters.csv", check.names=FALSE)
PI_POP <- pars$med[pars$param=="pi_population_median"]
strict <- read.csv("revision/case_definition_strictness.csv", check.names=FALSE, stringsAsFactors=FALSE)
gam <- read.csv("revision/tables/strictness_sensitivity/strictness_summary.csv", check.names=FALSE)$gamma_med[1]

anch <- strict[strict$set=="anchor", ]
S_MEAN <- mean(anch$strictness); S_SD <- sd(anch$strictness); GATE_MAX <- max(anch$gate)
cat(sprintf("gamma=%.2f  pop PPV=%.3f  anchor strictness mean=%.2f sd=%.2f  max anchor gate=%d\n",
            gam, PI_POP, S_MEAN, S_SD, GATE_MAX))

ppv_strict <- function(ob) {
  r <- strict[strict$study==ob, ]
  if (nrow(r)==0) return(list(pi=PI_POP, note="no coding -> pooled"))
  if (r$gate[1] > GATE_MAX)                      # gate beyond anchor support
    return(list(pi=PI_POP, note=sprintf("gate=%d unsupported -> held pooled", r$gate[1])))
  list(pi=plogis(qlogis(PI_POP) + gam*(r$strictness[1]-S_MEAN)/S_SD),
       note=sprintf("strictness=%d", r$strictness[1]))
}

PLAUSIBLE <- list(
  "Lewis 2005"=c("Additive","Multiplicative"), "Davis 2018"=c("Reporting-shock","Additive"),
  "N'Cho 2019"=c("Additive","Time-varying PPV","Reporting-shock"),
  "Muyembe-Tamfum 2009"=c("Multiplicative","Additive"),
  "Neil 2012"=c("Renewal-consistent","Multiplicative","Additive"),
  "Qamar 2018"=c("Renewal-consistent","Multiplicative"),
  "Yousafzai 2019"=c("Renewal-consistent","Multiplicative"),
  "Muti 2014"=c("Renewal-consistent","Multiplicative","Time-varying PPV"),
  "Polonsky 2014, Dzivaresekwa"=c("Renewal-consistent","Multiplicative","Time-varying PPV"),
  "Polonsky 2014_Kuwadzana"=c("Renewal-consistent","Multiplicative","Time-varying PPV"),
  "Ali 2017"=c("Time-varying PPV","Additive","Renewal-consistent"),
  "Aye 2004"=c("Additive","Multiplicative"),
  "Kabwama 2017"=c("Additive","Multiplicative","Renewal-consistent"))

rows <- list()
for (ob in names(prep$series)) {
  S <- prep$series[[ob]]; Tn <- length(S)
  anchored <- !is.na(pic$pi_med[match(ob, pic$study)])
  plaus <- PLAUSIBLE[[ob]]; if (is.null(plaus)) plaus <- names(ALLOCS)
  for (mode in c("pooled","strictness")) {
    if (anchored) { pi_o <- pic$pi_med[match(ob, pic$study)]; note <- "anchored" }
    else if (mode=="pooled") { pi_o <- PI_POP; note <- "pooled" }
    else { z <- ppv_strict(ob); pi_o <- z$pi; note <- z$note }
    for (nm in plaus) {
      I <- if (nm=="Time-varying PPV") ALLOCS[[nm]](S, pi_o, 2) else ALLOCS[[nm]](S, pi_o)
      Bkg <- pmax(S-I, 0)
      vc <- vaccinate(I, Bkg, w, T_EFF, coverage, psi_T, delta, cfg$timing$c_shape)
      rows[[length(rows)+1L]] <- data.frame(outbreak=ob, anchored=anchored, mode=mode, note=note,
        pi=pi_o, alloc=nm, susp_total=sum(S), true_total=sum(I),
        averted=sum(I)-sum(vc$Iv), stringsAsFactors=FALSE)
    }
  }
}
res <- do.call(rbind, rows)
write.csv(res, "revision/tables/ppv_assignment_comparison.csv", row.names=FALSE)

per <- do.call(rbind, lapply(split(res, list(res$outbreak,res$mode), drop=TRUE), function(d) data.frame(
  outbreak=d$outbreak[1], anchored=d$anchored[1], mode=d$mode[1], note=d$note[1], pi=d$pi[1],
  true_total=d$true_total[1], averted_lo=min(d$averted), averted_hi=max(d$averted), stringsAsFactors=FALSE)))
wide <- merge(per[per$mode=="pooled", c("outbreak","anchored","pi","true_total","averted_lo","averted_hi")],
              per[per$mode=="strictness", c("outbreak","note","pi","true_total","averted_lo","averted_hi")],
              by="outbreak", suffixes=c("_pool","_str"))
wide$pi_ratio <- wide$pi_str/wide$pi_pool
wide <- wide[order(wide$pi_ratio), ]

cat("\n=== unanchored PPV assignment: pooled vs strictness-conditioned ===\n")
u <- wide[!wide$anchored, ]
print(data.frame(outbreak=u$outbreak, note=u$note,
                 PPV=sprintf("%.3f -> %.3f (%.2fx)", u$pi_pool, u$pi_str, u$pi_ratio),
                 averted=sprintf("%.0f-%.0f -> %.0f-%.0f", u$averted_lo_pool,u$averted_hi_pool,
                                 u$averted_lo_str,u$averted_hi_str)), row.names=FALSE)

agg <- function(sfx) c(true=sum(wide[[paste0("true_total_",sfx)]]),
                       lo=sum(wide[[paste0("averted_lo_",sfx)]]), hi=sum(wide[[paste0("averted_hi_",sfx)]]))
A <- agg("pool"); B <- agg("str"); SUSP <- sum(unique(res[,c("outbreak","susp_total")])$susp_total)
cat(sprintf("\n=== AGGREGATE (production timing, week %d) ===\n", T_EFF))
cat(sprintf("suspected cases: %s\n", format(round(SUSP), big.mark=",")))
cat(sprintf("POOLED     : true typhoid %s | averted %s-%s | true red %.1f-%.1f%% | visible %.1f-%.1f%%\n",
    format(round(A['true']),big.mark=","), format(round(A['lo']),big.mark=","), format(round(A['hi']),big.mark=","),
    100*A['lo']/A['true'], 100*A['hi']/A['true'], 100*A['lo']/SUSP, 100*A['hi']/SUSP))
cat(sprintf("STRICTNESS : true typhoid %s | averted %s-%s | true red %.1f-%.1f%% | visible %.1f-%.1f%%\n",
    format(round(B['true']),big.mark=","), format(round(B['lo']),big.mark=","), format(round(B['hi']),big.mark=","),
    100*B['lo']/B['true'], 100*B['hi']/B['true'], 100*B['lo']/SUSP, 100*B['hi']/SUSP))
cat(sprintf("\nAggregate averted shift: %+.0f%% (lo) / %+.0f%% (hi);  true-typhoid denominator shift %+.0f%%\n",
    100*(B['lo']-A['lo'])/A['lo'], 100*(B['hi']-A['hi'])/A['hi'], 100*(B['true']-A['true'])/A['true']))
