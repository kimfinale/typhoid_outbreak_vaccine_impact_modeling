# final epidemic size
# the simplest SIR model
# R = 1 - exp(-R0*R) where R is the final epidemic size (or R(\infty) for the SIR model)
final_epidemic_size <- function(R0 = 2) {
  y = function(x) x - 1 + exp(-R0*x)
  final_size <- uniroot(y, interval=c(1e-6,1-1e-6))$root

  return(final_size)

}

# print parameter values
print_params <- function(params){
  n <- names(params)
  for(i in seq_along(n)){
    cat(paste0(n[i], "=", params[n[i]]), ", ")
  }
}

expit <- function(x) 1/(1+exp(-x))
logit <- function(x) {
  if(x>=1 | x<=0){
    stop("The function can take the numbers between 0 and 1 with both exclusive")
  }
  log(x/(1-x))
}

get_output_days <- function(x) {
  unitdays = gsub("^([0-9]+).*", "\\1", x$date_range)
  l = length(unique(unitdays))
  if (l != 1){
    stop("Output days have more than one kind")
  }
  return(as.double(unique(unitdays)))
}
#
my_discrete_colors <-
  c("dodgerblue2", "#E31A1C", "green4", "#6A3D9A",
    "#FF7F00","black", "gold1", "skyblue2", "palegreen2", "#FDBF6F",
    "gray70", "maroon", "orchid1", "darkturquoise", "darkorange4", "brown")


# tstamp
# create the time stamp yearmonthday by default and hour, minute, and second can be added
tstamp <- function(year=TRUE, month=TRUE, day=TRUE,
                   hour=FALSE, minute=FALSE, second=FALSE) {
  stamp1 <- c()
  stamp2 <- c()
  if (year & !month & !day) {
    stamp <- format(Sys.time(), "%Y")
  } else if (year & month & !day) {
    stamp1 <- format(Sys.time(), "%Y%m")
  } else if (year & month & day) {
    stamp1 <- format(Sys.time(), "%Y%m%d")
  } else if (!year & month & day) {
    stamp1 <- format(Sys.time(), "%m%d")
  } else if (year & !month & day) {
    stamp1 <- format(Sys.time(), "%Y%d")
  } else if (!year & month & !day) {
    stamp1 <- format(Sys.time(), "%m")
  } else if (!year & !month & day) {
    stamp1 <- format(Sys.time(), "%d")
  } else{ stamp1 <- "You'd better select parameters well."}

  if (hour & !minute & !second) {
    stamp2 <- format(Sys.time(), "%H")
  } else if (hour & minute & !second) {
    stamp2 <- format(Sys.time(), "%H%M")
  } else if (hour & minute & second) {
    stamp2 <- format(Sys.time(), "%H%M%S")
  } else if (!hour & minute & !second) {
    stamp2 <- format(Sys.time(), "%M")
  } else if (!hour & !minute & second) {
    stamp2 <- format(Sys.time(), "%S")
  } else if (!hour & minute & second) {
    stamp2 <- format(Sys.time(), "%M%S")
  } else{}

  if (!is.null(stamp2)) {
    stamp1 <- paste0(stamp1, "T", stamp2)
  }
  return (stamp1)
}

print_params <- function(params){
  n <- names(params)
  str <- paste0(n, "=", params[n])
  print(str)
}

# calculates alpha and beta parameters for the Beta distribution
calc_beta_params <- function(mu, var) {
  alpha <- ((1 - mu) / var - 1 / mu) * mu ^ 2
  beta <- alpha * (1 / mu - 1)
  return(params = list(alpha = alpha, beta = beta))
}


# calculates alpha and beta parameters for the Beta distribution
calc_lognorm_params <- function(mean, sd){
  v <- sd*sd
  m <- mean
  phi <- sqrt(v + m*m);
  mu <- log(m*m/phi);                # mean of log(Y)
  sigma <- sqrt(log(phi*phi/(m*m))); # std dev of log(Y)

  return(list(mu = mu, sigma = sigma))
}

calc_sd_lognorm <- function(mu, p, alpha=0.025, max=10) {
  eq <- function(x) exp(mu + x * qnorm(alpha)) - p
  uniroot(eq, interval=c(0,max))$root
}

erf <- function(x) 2 * pnorm(x * sqrt(2)) - 1
quan_func <- function(x) exp(mu + sqrt(2*sigma*sigma) / erf(2*p-1))
# MLE grid search -  1D ---------------------------------------------------------

## copied from https://stackoverflow.com/questions/29067916/error-function-erfz
## if you want the so-called 'error function'
erf <- function(x) 2 * pnorm(x * sqrt(2)) - 1
## (see Abramowitz and Stegun 29.2.29)
## and the so-called 'complementary error function'
erfc <- function(x) 2 * pnorm(x * sqrt(2), lower = FALSE)
## and the inverses
erfinv <- function (x) qnorm((1 + x)/2)/sqrt(2)
erfcinv <- function (x) qnorm(x/2, lower = FALSE)/sqrt(2)

MLE_check <- function(p_name = "local_rep_prop", theta_tab,nn=1e3){

  # theta_tab <- seq(0.001,0.01,0.001)
  store_lik <- NULL

  for(ii in 1:length(theta_tab)){

    theta[[p_name]] <- theta_tab[ii]

    # Run SMC and output likelihoods
    output_smc <- smc_model(theta,
                            nn=1e3 # number of particles
    )
    store_lik <- rbind(store_lik,c(theta_tab[ii],output_smc$lik))

  }

  colnames(store_lik) <- c("param","lik")
  store_lik <- as_tibble(store_lik)

}

# MLE grid search -  2D ---------------------------------------------------------

MLE_check_2D <- function(p1_name = "local_rep_prop", p2_name = "confirmed_prop",
                         theta_tab1, theta_tab2,nn=1e3,
                         filename = 1){

  # p1_name = "local_rep_prop"; p2_name = "confirmed_prop"; theta_tab1 = seq(0.01,0.05,0.01); theta_tab2 = seq(0.3,1,0.1)

  store_lik <- NULL

  for(ii in 1:length(theta_tab1)){

    for(jj in 1:length(theta_tab2)){

      theta[[p1_name]] <- theta_tab1[ii]
      theta[[p2_name]] <- theta_tab2[jj]

      # Run SMC and output likelihooda
      output_smc <- smc_model(theta,
                              nn=1e3 # number of particles
      )
      store_lik <- rbind(store_lik,c(theta_tab1[ii],theta_tab2[jj],output_smc$lik))

    }
  }

  colnames(store_lik) <- c("param1","param2","lik")
  store_lik <- as_tibble(store_lik)

  write_csv(store_lik,paste0("outputs/param_search_",filename,".csv"))

}

# MLE grid search -  2D ---------------------------------------------------------

MLE_check_3D <- function(p1_name = "local_rep_prop",
                         p2_name = "confirmed_prop",
                         p3_name = "betavol",
                         theta_tab1,
                         theta_tab2,
                         theta_tab3,
                         nn=1e3,
                         filename = 1){

  # p1_name = "local_rep_prop"; p2_name = "confirmed_prop"; p3_name = "betavol"; theta_tab1 = seq(0.01,0.05,0.02); theta_tab2 = seq(0.6,1,0.2); theta_tab3 = seq(0.1,0.3,0.1)

  store_lik <- NULL

  out_fit <- foreach(ii = 1:length(theta_tab1)) %dopar% {
    #for(ii in 1:length(theta_tab1)){
    for(jj in 1:length(theta_tab2)){
      for(kk in 1:length(theta_tab3)){
        theta[[p1_name]] <- theta_tab1[ii]
        theta[[p2_name]] <- theta_tab2[jj]
        theta[[p3_name]] <- theta_tab3[kk]
        # Run SMC and output likelihooda
        output_smc <- smc_model(theta,
                                nn=1e3 # number of particles
        )
        store_lik <- rbind(store_lik,c(theta_tab1[ii],theta_tab2[jj],theta_tab3[kk],output_smc$lik))
      }
    }
    store_lik
  }
  # Collate results
  store_lik <- NULL
  for(ii in 1:length(theta_tab1)){

    store_lik <- rbind(store_lik,out_fit[[ii]])

  }

  colnames(store_lik) <- c("param1","param2","param3","lik")
  store_lik <- as_tibble(store_lik)

  write_csv(store_lik,paste0("outputs/param_search_",filename,".csv"))

}

# Compute acceptance probability ------------------------------------------



figure_size <- data.frame(journal=c("Nature","Elsevier","Lancet"),
                          single=c(89,90,75),
                          double=c(183,190,154),
                          unit2=c("mm","mm","mm"))




seir_jl <- function(u, p, t){

  S <- u[1];
  E <- u[2];
  # the number of compartments to model the duration of infectiousness
  I <- u[3];
  R <- u[4];
  C <- u[5];

  # population size
  pop <- S+E+I+R

  epsilon <- p[1] # 1/latent period
  gamma <- p[2] # 1/duration of infectiousness
  beta <- p[3] # transmission rate
  omega <- p[4] # 1/omega = duration of natural immunity
  # force of infection
  foi <- beta*I/pop

  muEI <- epsilon
  muIR <- gamma
  muRS <- omega

  # differential equations
  dS <- - foi*S + muRS*R
  dE <- foi*S - muEI*E
  dI <- muEI*E - muIR*I
  dR <- muIR*I - muRS*R

  dC <- muEI*E

  return(c(dS,dE,dI,dR,dC))
}

run_seir_jl <- function(model=NULL,
                        epsilon=1/4,
                        gamma=1/7,
                        omega=1/(4*365),
                        R0=2,
                        tend=100,
                        saveat=1,
                        nms=c("t","S","E","I","R","C")){
  # library(diffeqr)
  de <- diffeqr::diffeq_setup()
  # simulation
  # R0 <- 2 # basic reproduction number
  # epsilon <- 1/4 # 1/epsilon = incubation period
  # gamma <- 1/7 # 1/gamma = duration of infectiousness
  beta <- R0*gamma # instantaneous transmission rate
  # omega <- 1/(4*365) # natural immunity waning rate
  # parameters
  params <- c(epsilon=epsilon, gamma=gamma, beta=beta, omega=omega)
  u0 <- c(0.99, 0, 0.01, 0, 0)
  tend <- 100 #
  tspan <- c(0.0, tend)

  prob <- de$ODEProblem(model, u0, tspan, params)
  sol <- de$solve(prob, de$Tsit5(), saveat=saveat)
  mat <- sapply(sol$u, identity)
  udf <- as.data.frame(t(mat))
  out <- cbind(data.frame(t=sol$t), udf)
  names(out) <- nms

  return(out)
}



# ggplot2 themes


function (base_size = 11, base_family = "", base_line_size = base_size/22,
          base_rect_size = base_size/22)
{
  theme_grey(base_size = base_size, base_family = base_family,
             base_line_size = base_line_size, base_rect_size = base_rect_size) %+replace%
    theme(panel.background = element_rect(fill = "white", colour = NA),
          panel.border = element_rect(fill = NA, colour = "grey20"),
          panel.grid = element_line(colour = "grey92"),
          panel.grid.minor = element_line(linewidth = rel(0.5)),
          strip.background = element_rect(fill = "grey85",
                                          colour = "grey20"),
          legend.key = element_rect(fill = "white",
                                                                                        colour = NA), complete = TRUE)
}

theme_pub <- function(base_size=12, base_family="sans") {
  library(grid)
  library(ggthemes)
  (theme_foundation(base_size=base_size, base_family=base_family)
    + theme(plot.title=element_text(face="bold",
                                    size=rel(1.2),
                                    hjust=0.5,
                                    margin=margin(0,0,20,0)),
            text = element_text(),
            panel.background = element_rect(colour = NA),
            plot.background = element_rect(colour = NA),
            panel.border = element_rect(colour = NA),
            axis.title = element_text(face="plain", size=rel(1)),
            axis.title.y = element_text(angle=90, vjust=2),
            axis.title.x = element_text(vjust=-0.2),
            axis.text = element_text(),
            axis.line.x = element_line(colour="black"),
            axis.line.y = element_line(colour="black"),
            axis.ticks = element_line(),
            panel.grid.major = element_line(colour="#f0f0f0"),
            panel.grid.minor = element_blank(),
            legend.key = element_rect(colour = NA),
            legend.position = "bottom",
            legend.direction = "horizontal",
            legend.box = "vetical",
            legend.key.size= unit(3,"mm"),
            legend.margin = unit(0, "cm"),
            legend.title = element_text(face="plain"),
            plot.margin=unit(c(3,3,3,3),"mm"),
            strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
            strip.text = element_text(face="plain")
    ))
}

library(ggplot2)
library(grid)

# define consistent ggplot theme to apply to all figures
theme_ms <- function(base_size=12, base_family="Helvetica") {
  library(grid)
  (theme_bw(base_size = base_size, base_family = base_family)+
      theme(text=element_text(color="black"),
            axis.title=element_text(face="bold", size = rel(1.3)),
            axis.text=element_text(size = rel(1), color = "black"),
            legend.title=element_text(face="bold"),
            legend.text=element_text(face="bold"),
            legend.background=element_rect(fill="transparent"),
            legend.key.size = unit(0.8, 'lines'),
            panel.border=element_rect(color="black",size=1),
            panel.grid=element_blank()
      ))
}

theme_Publication <- function(base_size=14, base_family="sans") {
  library(grid)
  library(ggthemes)
  (theme_foundation(base_size=base_size, base_family=base_family)
    + theme(plot.title = element_text(face = "bold",
                                      size = rel(1.2), hjust = 0.5, margin = margin(0,0,20,0)),
            text = element_text(),
            panel.background = element_rect(colour = NA),
            plot.background = element_rect(colour = NA),
            panel.border = element_rect(colour = NA),
            axis.title = element_text(face = "bold",size = rel(1)),
            axis.title.y = element_text(angle=90,vjust =2),
            axis.title.x = element_text(vjust = -0.2),
            axis.text = element_text(),
            axis.line.x = element_line(colour="black"),
            axis.line.y = element_line(colour="black"),
            axis.ticks = element_line(),
            panel.grid.major = element_line(colour="#f0f0f0"),
            panel.grid.minor = element_blank(),
            legend.key = element_rect(colour = NA),
            legend.position = "bottom",
            legend.direction = "horizontal",
            legend.box = "vetical",
            legend.key.size= unit(0.5, "cm"),
            #legend.margin = unit(0, "cm"),
            legend.title = element_text(face="italic"),
            plot.margin=unit(c(10,5,5,5),"mm"),
            strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
            strip.text = element_text(face="bold")
    ))

}

scale_fill_Publication <- function(...){
  library(scales)
  discrete_scale("fill","Publication",manual_pal(values = c("#386cb0","#f87f01","#7fc97f","#ef3b2c","#feca01","#a6cee3","#fb9a99","#984ea3","#8C591D")), ...)

}

scale_colour_Publication <- function(...){
  library(scales)
  discrete_scale("colour","Publication",manual_pal(values = c("#386cb0","#f87f01","#7fc97f","#ef3b2c","#feca01","#a6cee3","#fb9a99","#984ea3","#8C591D")), ...)

}


### Dark theme for ggplot plots

theme_dark_grey <- function(base_size=14, base_family="sans") {
  library(grid)
  library(ggthemes)
  (theme_foundation(base_size=base_size, base_family=base_family)
    + theme(plot.title = element_text(face = "bold", colour = '#ffffb3',
                                      size = rel(1.2), hjust = 0.5, margin = margin(0,0,20,0)),
            text = element_text(),
            panel.background = element_rect(colour = NA, fill = 'grey20'),
            plot.background = element_rect(colour = NA, fill = '#262626'),
            panel.border = element_rect(colour = NA),
            axis.title = element_text(face = "bold",size = rel(1), colour = 'white'),
            axis.title.y = element_text(angle=90,vjust =2),
            axis.title.x = element_text(vjust = -0.2),
            axis.text = element_text(colour = 'grey70'),
            axis.line.x = element_line(colour="grey70"),
            axis.line.y = element_line(colour="grey70"),
            axis.ticks = element_line(colour="grey70"),
            panel.grid.major = element_line(colour="#262626"),
            panel.grid.minor = element_blank(),
            legend.background = element_rect(fill ='#262626'),
            legend.text = element_text(color = 'white'),
            legend.key = element_rect(colour = NA, fill = '#262626'),
            legend.position = "bottom",
            legend.direction = "horizontal",
            legend.box = "vetical",
            legend.key.size= unit(0.5, "cm"),
            #legend.margin = unit(0, "cm"),
            legend.title = element_text(face="italic", colour = 'white'),
            plot.margin=unit(c(10,5,5,5),"mm"),
            strip.background=element_rect(colour="#2D3A4C",fill="#2D3A4C"),
            strip.text = element_text(face="bold", colour = 'white')
    ))
}

scale_fill_Publication_dark <- function(...){
  library(scales)
  discrete_scale("fill","Publication",manual_pal(values = c("#fbb4ae","#b3cde3","#ccebc5","#decbe4","#fed9a6","#ffffcc","#e5d8bd","#fddaec","#f2f2f2")), ...)

}

scale_colour_Publication_dark <- function(...){
  library(scales)
  discrete_scale("colour","Publication",manual_pal(values = c("#fbb4ae","#b3cde3","#ccebc5","#decbe4","#fed9a6","#ffffcc","#e5d8bd","#fddaec","#f2f2f2")), ...)

}


# theme_transparent <- function(base_size=14, base_family="sans") {
#    library(grid)
#    library(ggthemes)
#    (theme_foundation(base_size=base_size, base_family=base_family)
#       + theme(plot.title = element_text(face = "bold", colour = '#ffffb3',
#                                         size = rel(1.2), hjust = 0.5),
#               text = element_text(),
#               panel.background = element_rect(colour = NA, fill = 'transparent'),
#               plot.background = element_rect(colour = NA, fill = 'transparent'),
#               panel.border = element_rect(colour = NA),
#               axis.title = element_text(face = "bold",size = rel(1), colour = 'white'),
#               axis.title.y = element_text(angle=90,vjust =2),
#               axis.title.x = element_text(vjust = -0.2),
#               axis.text = element_text(colour = 'grey70'),
#               axis.line.x = element_line(colour="grey70"),
#               axis.line.y = element_line(colour="grey70"),
#               axis.ticks = element_line(colour="grey70"),
#               panel.grid.major = element_line(colour="#262626"),
#               panel.grid.minor = element_blank(),
#               legend.background = element_rect(fill = 'transparent'),
#               legend.text = element_text(color = 'white'),
#               legend.key = element_rect(colour = NA, fill = 'grey20'),
#               legend.position = "bottom",
#               legend.direction = "horizontal",
#               legend.box = "vetical",
#               legend.key.size= unit(0.5, "cm"),
#               #legend.margin = unit(0, "cm"),
#               legend.title = element_text(face="italic", colour = 'white'),
#               plot.margin=unit(c(10,5,5,5),"mm"),
#               strip.background=element_rect(colour="#2D3A4C",fill="#2D3A4C"),
#               strip.text = element_text(face="bold", colour = 'white')
#       ))
# }

theme_2D_contour <- function(){
  theme(panel.background=element_blank(),
      plot.background=element_blank(),
      # panel.border=element_blank(),
      panel.border = element_rect(fill = NA, colour = "grey20"),
      # axis.line.x=element_line(colour="black"),
      # axis.line.y=element_line(colour="black"),
      axis.text=element_text(size=11),
      axis.title=element_text(size=11),
      legend.position = "bottom",
      legend.direction = "horizontal",
      legend.box = "vetical",
      legend.key.size = unit(1.2,"mm"),
      legend.text = element_text(size=10),
      legend.key.height = unit(3,"mm"),
      legend.key.width = unit(2,"mm"),
      legend.margin = margin(-2, 0, 0, 0, unit="mm"),
      legend.title = element_text(face="plain",size=10))
}

theme_dark_blue <- function(base_size=14, base_family="sans") {
  library(grid)
  library(ggthemes)
  (theme_foundation(base_size=base_size, base_family=base_family)
    + theme(plot.title = element_text(face = "bold", colour = '#ffffb3',
                                      size = rel(1.2), hjust = 0.5, margin = margin(0,0,20,0)),
            text = element_text(),
            panel.background = element_rect(colour = NA, fill = '#282C33'),
            plot.background = element_rect(colour = NA, fill = '#282C33'),
            panel.border = element_rect(colour = NA),
            axis.title = element_text(face = "bold",size = rel(1), colour = 'white'),
            axis.title.y = element_text(angle=90,vjust =2),
            axis.title.x = element_text(vjust = -0.2),
            axis.text = element_text(colour = 'grey70'),
            axis.line.x = element_line(colour="grey70"),
            axis.line.y = element_line(colour="grey70"),
            axis.ticks = element_line(colour="grey70"),
            panel.grid.major = element_line(colour="#343840"),
            panel.grid.minor = element_blank(),
            legend.background = element_rect(fill ='#282C33'),
            legend.text = element_text(color = 'white'),
            legend.key = element_rect(colour = NA, fill = '#282C33'),
            legend.position = "bottom",
            legend.direction = "horizontal",
            legend.box = "vetical",
            legend.key.size= unit(0.5, "cm"),
            #legend.margin = unit(0, "cm"),
            legend.title = element_text(face="italic", colour = 'white'),
            plot.margin=unit(c(10,5,5,5),"mm"),
            strip.background=element_rect(colour="#2D3A4C",fill="#2D3A4C"),
            strip.text = element_text(face="bold", colour = 'white')
    ))
}


