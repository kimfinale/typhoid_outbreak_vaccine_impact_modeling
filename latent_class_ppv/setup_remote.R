# =============================================================================
# One-time remote-machine setup for the manuscript pipeline.
# RUN THIS AFTER installing (OS-level, cannot be done from R):
#   1. R >= 4.5   2. RTools45 (Windows; matching R 4.5)   3. (optional) Quarto CLI
# It installs cmdstanr + cmdstan + the R packages the pipeline needs, and verifies
# the toolchain with a trivial compile+sample.
#   Rscript latent_class_ppv/setup_remote.R
# =============================================================================
cat("R version:", R.version.string, "\n")
if (getRversion() < "4.4") stop("R is too old (need >= 4.4, ideally 4.5.x). Update R + RTools first.")

repos <- c("https://mc-stan.org/r-packages/", getOption("repos"))
if (is.na(repos[2]) || repos[2] == "@CRAN@") repos[2] <- "https://cloud.r-project.org"

pkgs <- c("cmdstanr","posterior","ggplot2","dplyr","tidyr","yaml","scales","knitr","rmarkdown","quarto")
for (p in pkgs) if (!requireNamespace(p, quietly = TRUE)) {
  cat("installing", p, "...\n"); install.packages(p, repos = repos)
}

suppressMessages(library(cmdstanr))
cat("\n-- C++ toolchain check (needs RTools45 on Windows) --\n")
tc <- tryCatch(check_cmdstan_toolchain(fix = FALSE), error = function(e) conditionMessage(e))
print(tc)

have_cmdstan <- !is.null(tryCatch(cmdstan_path(), error = function(e) NULL))
if (!have_cmdstan) { cat("\ninstalling cmdstan (a few minutes)...\n"); install_cmdstan(cores = 2) }
cat("cmdstan version:", tryCatch(cmdstan_version(), error = function(e) "NOT INSTALLED"), "\n")

cat("\n-- trivial compile + sample test --\n")
ok <- tryCatch({
  f <- write_stan_file("parameters { real y; } model { y ~ normal(0, 1); }")
  m <- cmdstan_model(f)
  fit <- m$sample(chains = 1, iter_warmup = 100, iter_sampling = 100, refresh = 0, show_messages = FALSE)
  nrow(fit$draws()) > 0
}, error = function(e) { cat("FAILED:", conditionMessage(e), "\n"); FALSE })
cat(if (ok) "TOOLCHAIN OK — ready to run the pipeline.\n" else
    "TOOLCHAIN NOT READY — check RTools/compiler (see message above).\n")

cat("\nquarto CLI:", ifelse(nzchar(Sys.which("quarto")), Sys.which("quarto"),
    "NOT FOUND (install Quarto CLI for PDF/docx; else render falls back to rmarkdown/.md)"), "\n")
