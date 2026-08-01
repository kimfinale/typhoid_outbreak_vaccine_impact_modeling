# =============================================================================
# Build the outbreak epidemiologic + transmission-classification summary table.
#   Rscript revision/make_outbreak_epi_summary.R
# Sources: source_decomposition/tables/tab_transmission_mode.csv (theta, mode,
# evidence), renewal/tables/tab_resolution.csv (resolution), and primary-source
# verifications for the 6 monthly series (Lutterloh, Imanishi, Walters, Hendriksen).
# Transmission-route indicators: 2=primary/strong, 1=contributory, 0=not evidenced.
# =============================================================================
setwd(Sys.getenv("RENEWAL_ROOT", "."))
OUT <- "revision"; dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

d <- data.frame(stringsAsFactors = FALSE,
  outbreak = c("Lewis 2005","Aye 2004","Muyembe-Tamfum 2009","Ali 2017","Kabwama 2017",
               "Lutterloh 2012","Hendriksen 2015","Davis 2018","N'Cho 2019",
               "Neil 2012","Walters 2014 (Kasese)","Walters 2014 (Bundibugyo)",
               "Polonsky 2014 (Dzivaresekwa)","Polonsky 2014 (Kuwadzana)","Muti 2014",
               "Imanishi 2014 (Dzivaresekwa)","Imanishi 2014 (Kuwadzana)",
               "Qamar 2018","Yousafzai 2019"),
  location = c("Kathmandu","Madaya, Mandalay","Kinshasa","Kikwit","Kampala",
               "Neno/Tsangano (border)","Lusaka","Harare","Harare",
               "Kasese","Kasese","Bundibugyo","Harare (Dzivaresekwa)","Harare (Kuwadzana)",
               "Harare (Dzivaresekwa)","Harare (Dzivaresekwa)","Harare (Kuwadzana)",
               "Hyderabad","Hyderabad"),
  country = c("Nepal","Myanmar","DR Congo","DR Congo","Uganda","Malawi/Mozambique","Zambia",
              "Zimbabwe","Zimbabwe","Uganda","Uganda","Uganda","Zimbabwe","Zimbabwe","Zimbabwe",
              "Zimbabwe","Zimbabwe","Pakistan","Pakistan"),
  period = c("2002","2000","2004","2017","Jan-Jun 2015","Mar-Nov 2009","2010-2012",
             "2016-2017","2017-2018","2008-2009","2009-2011","2011","2012","2012","2011",
             "2011-2012","2011-2012","2016-2017","2016-2018"),
  suspected = c("~5,900","49","peritonitis series","~1,500","10,230","303","~2,040",
                "~880","~580","709 (+ IP-heavy)","709","333","916","1,860","1,078",
                "4,181 (Harare total)","4,181 (Harare total)","~586","101 (confirmed)"),
  bc_confirmed = c("culture-confirmed","3/22","11/16 (peritonitis)","culture-confirmed","56/364",
                   "44/108","94 isolates (pooled)","borehole-linked","74/286 (mixed spec.)",
                   "19/54","7/74","15/72","62 pooled/3,795","(pooled)","24 pooled/1,078",
                   "52 pooled/4,181","(pooled)","70/170 (XDR)","101/101 (confirmed)"),
  vehicle = c("Single municipal water supply","Myaung River water + contact","Presumed municipal water",
              "Camp water + overcrowding","Vended/spring water (recurrent)","Spring water",
              "Water/sanitation (endemic)","Boreholes","Boreholes/shallow wells","No vehicle implicated",
              "Piped water (taps/rivers)","Gravity-flow-scheme piped water","Borehole water","Borehole water",
              "Wells","Boreholes/shallow wells","Boreholes/shallow wells","Sewage-contaminated water",
              "Sewage-contaminated water"),
  # route indicators: point/common-source ; environmental (waterborne, long-cycle) ; person-to-person
  pt_source = c(2,1,1,2,1,1,0,1,1,0,0,2,0,0,0,1,1,0,0),
  env_water = c(2,2,1,2,2,2,2,2,2,1,2,2,2,2,2,2,2,2,2),
  p2p       = c(0,2,0,1,1,1,2,0,0,2,1,1,2,2,2,0,0,2,2),
  evidence = c("S. Typhi from tap PFGE-matched to cases; ended on chlorination (point-source)",
               "Water aOR 12.5 AND person-to-person aOR 22; small township",
               "Untested; severe peritonitis case series (different estimand)",
               "Military-camp water then dissemination to general population",
               "Recurrent contamination sustained transmission over months; adults 15-59",
               "E. coli in springs; single PFGE clone (35/42); prolonged ~8 months",
               "Single H58 clone 83%; multi-peak 33-month endemic/propagated",
               "13/18 boreholes fecal-coliform+; declined after repair (common-source)",
               "Boreholes/wells amid water shortage; WASH-terminated",
               "No vehicle; multi-strain PFGE; prolonged ~19 mo; propagated",
               "E. coli taps/rivers; NO clustering by water source; continuous multi-year",
               "Single GFS implicated; clustering; E. coli intake->taps (common-source)",
               "Borehole; 'long-cycle transmission predominated'; poor common-source fit (r2 0.07)",
               "Borehole; epicentre shifted here; long-cycle",
               "Wells E. coli (8/8) AND home-contact aOR 8.3 (propagated)",
               "Case-control implicated water; boreholes E. coli; no secondary analysis",
               "Same (Harare 2011-12 epidemic)",
               "Sewage S. Typhi DNA; 2-weekly growing peaks -> strongly propagated (r2 0)",
               "Cyclical curve -> propagated XDR clone"),
  theta_cs = c(0.80,0.50,0.65,0.50,0.60,0.55,0.45,0.725,0.80,0.20,0.20,0.70,0.50,0.50,0.50,0.50,0.50,0.50,0.50),
  mode = c("A common-source","C mixed","A common-source","C mixed","C mixed","C mixed","C mixed",
           "A common-source","A common-source","B propagated","B propagated (=Neil)","A common-source",
           "C mixed (propagated-lean)","C mixed","C mixed","A/C water","A/C water",
           "C mixed (propagated, data theta 0.05)","C mixed (propagated, data theta 0.38)"),
  confidence = c("high","high","low","med","med","med","med(est.)","med","med","med",
                 "med(est.)","med(est.)","med","med","high","med(est.)","med(est.)","high","med"),
  resolution = c("daily","weekly","daily","weekly","daily","fortnightly","monthly","daily","daily",
                 "weekly","monthly","monthly","weekly","weekly","daily","monthly","monthly","weekly","weekly"),
  model = c("Additive renewal-with-source (source-dominant)","Additive renewal-with-source",
            "exclude (complications)","Additive renewal-with-source (+spatial)",
            "Additive renewal-with-source","Static only","Static only (+ODE)",
            "Additive renewal-with-source (+ODE)","Additive renewal-with-source (+ODE)",
            "Additive renewal-with-source (propagated-dominant)","Static (redundant)",
            "Static only","Additive renewal-with-source (propagated-lean)",
            "Additive renewal-with-source","Additive renewal-with-source (redundant)",
            "Static (redundant)","Static (redundant)",
            "Additive renewal-with-source (propagated-lean)",
            "Additive renewal-with-source (redundant)"),
  indep_group = c("independent","independent","independent","independent","independent","independent",
                  "independent","Harare-recurrent (season 1)","Harare-recurrent (season 2)",
                  "KASESE-UGANDA (canonical)","KASESE-UGANDA (=Neil)","KASESE-UGANDA (linked spread)",
                  "HARARE 2011-12 (canonical)","HARARE 2011-12","HARARE 2011-12",
                  "HARARE 2011-12","HARARE 2011-12","HYDERABAD-XDR (canonical)","HYDERABAD-XDR"))

write.csv(d, file.path(OUT, "outbreak_epi_summary.csv"), row.names = FALSE)

## ---- render a markdown table ------------------------------------------------
sym <- function(x) c("–","◐","●")[x + 1]     # 0 - , 1 half, 2 full
md <- c("| Outbreak | Location, country | Period | Suspected | Vehicle / source | Pt-src | Env-water | P2P | theta (mode) | Res. | Model | Independence |",
        "|---|---|---|---|---|:--:|:--:|:--:|---|---|---|---|")
for (i in seq_len(nrow(d))) md <- c(md, sprintf(
  "| %s | %s, %s | %s | %s | %s | %s | %s | %s | %.2f (%s) | %s | %s | %s |",
  d$outbreak[i], d$location[i], d$country[i], d$period[i], d$suspected[i], d$vehicle[i],
  sym(d$pt_source[i]), sym(d$env_water[i]), sym(d$p2p[i]), d$theta_cs[i], d$mode[i],
  d$resolution[i], d$model[i], d$indep_group[i]))
writeLines(c("# Outbreak epidemiologic & transmission-classification summary",
             "", "Route indicators: ● primary/strong  ◐ contributory  – not evidenced.",
             "theta = common-source fraction (source_decomposition); 1-theta = propagated.",
             "Model: additive renewal-with-source for resolution-eligible series; theta allocates source and propagated incidence inside the recursion. Monthly/fortnightly -> static only.",
             "", md, "",
             "**Independent epidemics: ~12** (from 19 series). Groups collapse: Kasese-Uganda (Neil+Walters x2), Harare 2011-12 (Polonsky x2 + Muti + Imanishi x2), Hyderabad-XDR (Qamar+Yousafzai); Muyembe excluded (complications)."),
           file.path(OUT, "outbreak_epi_summary.md"))
cat("Wrote", file.path(OUT, "outbreak_epi_summary.csv"), "and .md\n")
cat(sprintf("%d series; independent epidemics ~12\n", nrow(d)))
