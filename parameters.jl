# initial values
using LabelledArrays
using NamedTupleTools
using DifferentialEquations

function adjust_pop_dist_seiars!(p)
  # adjust the population distribution based on population, i0, and s0
  pop = p.population * p.n0;
  prop_S = p.s0;
  S = pop * prop_S; # can represent a proportion of the population who is at risk of infection
  I = S * p.i0; # parameterized in a way that a fraction of S (i.e, pop at risk) is I
  prop_u5 = p.prop_children; #
  p.u0.S_1 = (S - I) * prop_u5; # note that S - I is the number of susceptible
  p.u0.S_2 = (S - I) * (1 - prop_u5);
  # I and R states were modelled with two compartments. We simplify and equally divide the pop into two compartments
  p.u0.R_1 = pop * (1-prop_S) * prop_u5;
  p.u0.R_2 = pop * (1-prop_S) * (1 - prop_u5);
  
  latent_pd = 1/p.epsilon;
  infect_pd = 1/p.gamma;
  fE = latent_pd/(latent_pd + infect_pd); # approximate fraction of E, i.e., fE for every I

  p.u0.E_1 = I * fE * prop_u5;
  p.u0.E_2 = I * fE * (1 - prop_u5);

  p.u0.I_1 = I * (1-fE) * (1-p.fA) * prop_u5;
  p.u0.I_2 = I * (1-fE) * (1-p.fA) * (1 - prop_u5);
  
  p.u0.A_1 = I * (1-fE) * p.fA * prop_u5;
  p.u0.A_2 = I * (1-fE) * p.fA * (1 - prop_u5);
  return p    
end  

function update_LArray!(la1, la2)
  for k in keys(la2)
    la1[k] = la2[k];
  end
  return la1
end

# variables updated to include vaccinees with one or two doses who are infected
function initialize_u0_seiars(p=nothing)
  u0_1 = LVector(S_1=0.99, E_1=0.0, I_1=0.005, A_1=0.0, R_1=0.0, CE_1=0.0, CI_1=0.0, RS_1=0);
  u0_2 = LVector(S_2=0.99, E_2=0.0, I_2=0.005, A_2=0.0, R_2=0.0, CE_2=0.0, CI_2=0.0, RS_2=0);
  u0 = [u0_1*0.3; u0_2*0.7];
  if !isnothing(p)
    u0 = update_LArray!(u0, p)
  end
  return u0
end

function initialize_params_seiars(p = nothing)
  param_vector = (
    u0 = initialize_u0_seiars(),
    s0 = 1.0, 
    i0 = 0.01,
    R0 = 3.0,
    n0 = 1.0, # effective proportion of population (population at risk = population * n0)
    population = 1.0,
    prop_detection = 1.0, #
    # baseline population size is critical for the outbreak size
    # and the currently is the the population size of the admin in which the
    # outbreak was reported. Actual population is likely to be a fraction
    # prop_eff_pop = 1 #
    # time window over which the number of cases is tracked (i.e., weekly reported)
    # obs_length refers the duration during which data are available in days
    # this is updated for each data set
    tend = 100.0,
    report_freq = 1.0,
    obs_length = 365.0, # 20 weeks
    # this refers to the total simulation days and has to be larger than or equal
    # to the obs_length
    # tau = 0.01, # time step size for numerical integration
    epsilon = 1/4, # mean latent period = 1/epsilon
    gamma = 1/14, # mean infectious period = 1/gamma
    fA = 0.5, # fraction of asymptomatic state
    bA = 0.05, # relative infectiousness of asymptomatic state
    kappa = 575.0, # excretion rate cells per person per day
    prop_children = 0.193, #
  # exponent used to model sub-exponential growth (0 < expon < 1) (foi=I^{expon}*S/N)
    expon = 1.0, # exponential growth when expo = 1
    xi = 1/21, # mean decay rate of Vibrio cholerae
    K = 10000.0, # half-infective bacteria dose (10,000 cells/ml)
    R0W = 0.0, #
    sigma = 1/(4*365), # 1/sigma = mean duration of natural immunity
    # 1/sigma_vacc_1d = mean duration of OCV-induced immunity (1st dose)
    sigma_vacc_1d = 1/(2*365),
    # 1/sigma_vacc_2d = mean duration of OCV-induced immunity (2nd dose)
    sigma_vacc_2d = 1/(4*365),
    # vaccination-campaign associated
    vacc_1d_cov = 0.8,
    vacc_2d_cov = 0.8,
    campaign_1d_start = 1_000.0,
    campaign_2d_start = 1_000.0,
    campaign_1d_dur = 7.0, # duration of vaccination campaign in days
    campaign_2d_dur = 7.0, # delay_until_2nd_campaign = 14 # 30 days of delay between the 1st and the 2nd campaign
    vacc_1d_eff_1 = 0.3, # younger age group
    vacc_1d_eff_2 = 0.6, # 
    vacc_2d_eff_1 = 0.3, # younger age group
    vacc_2d_eff_2 = 0.6, # 
    vacc_1d_immunity_rate = 1/14, # vaccine-induced immunity develops after 14d on average after the first dose
    vacc_2d_immunity_rate = 1/7, # vaccine-induced immunity develops after 7d on average after the second dose
    # case threshold over which intervention will be implemented
    vacc_1d_times = LVector(start=0.0, stop=0.0),
    vacc_2d_times = LVector(start=0.0, stop=0.0),
    vacc_rates = LVector(vacc_1d_rate=0.0, vacc_2d_rate=0.0), # array instead of tuple is used such that vacc_*d_rate can be modified externally
    case_threshold = 10.0,
    alpha = 0.0, # proportional reduction in R0
    case_track_window = 7.0, 
    callback = nothing,
    ode = nothing,
    prob_ode = nothing,
    ode_solver = AutoTsit5(Rosenbrock23())
  )
  if !isnothing(p)
    param_vector = merge(param_vector, p)
  end
  return param_vector
end

