using Printf
function seiars_2ag!(du, u, p, t)
  # initial values
  prop_detection = p.prop_detection;
  epsilon = p.epsilon; # 1 / latent period
  kappa = p.kappa; # excretion rate
  xi = p.xi; # decay rate
  K = p.K; # half-infective dose (eg, 10,000 doses/ml)
  gamma = p.gamma; # 1 / recovery period
  sigma = p.sigma; # 1 / duration of natural immunity
  sigma_O = p.sigma_vacc_1d; # 1 / duration of OCV-induced immunity (one dose)
  sigma_T = p.sigma_vacc_2d; # 1 / duration of OCV-induced immunity (two doses)
  # 1 / rate from pre-symptomatic (P) state to infectious (I) states
  fA = p.fA; # fraction asymptomatic
  bA = p.bA; # relative infectivity of A to I
  R0 = p.R0;
  expon = p.expon; # exponent (0< & <1) to control the exponential-ness of the FOI
 
  β = R0 / ((bA*fA + (1-fA))/gamma);
  #  βW = β;
  # derived variables 
  N_1 = u.S_1 + u.E_1 + u.I_1 + u.A_1 + u.R_1;
  N_2 = u.S_2 + u.E_2 + u.I_2 + u.A_2 + u.R_2;
  N = N_1 + N_2;

  isum = u.I_1 + u.I_2;
  asum = u.A_1 + u.A_2;
  # @printf("expon = %.6f\n", expon);  
  foi = β * (isum + bA * asum) / N;
  # if isum > 0 && asum > 0
  #   foi = β * (isum^expon + bA*asum^expon) / N;   
  # end
  # @show p.vacc_rate
  # @printf("t = %.4f, vacc rate = %.4f\n", t, p.vacc_rate.first);

  S_1toE_1 = u.S_1 * foi;
  RS_1toE_1 = u.RS_1 * foi;
  E_1toI_1 = u.E_1 * (1-fA) * epsilon;
  E_1toA_1 = u.E_1 * fA * epsilon;
  I_1toR_1 = u.I_1 * gamma;
  A_1toR_1 = u.A_1 * gamma;
  R_1toRS_1 = u.R_1 * sigma;
  
  # Differential equations
  du.S_1 = - S_1toE_1;
  du.E_1 = + S_1toE_1 + RS_1toE_1 - E_1toI_1 - E_1toA_1;
  du.I_1 = + E_1toI_1 - I_1toR_1;
  du.A_1 = + E_1toA_1 - A_1toR_1;
  du.R_1 = + I_1toR_1 + A_1toR_1 - R_1toRS_1;
  du.RS_1 = + R_1toRS_1 - RS_1toE_1;
      
  du.CE_1 = + S_1toE_1;
  du.CI_1 = + prop_detection * E_1toI_1;
  
  S_2toE_2 = u.S_2 * foi;
  RS_2toE_2 = u.RS_2 * foi;
  E_2toI_2 = u.E_2 * (1-fA) * epsilon;
  E_2toA_2 = u.E_2 * fA * epsilon;
  I_2toR_2 = u.I_2 * gamma;
  A_2toR_2 = u.A_2 * gamma;
  R_2toRS_2 = u.R_2 * sigma;

  # Differential equations
  du.S_2 = - S_2toE_2;
  du.E_2 = + S_2toE_2 + RS_2toE_2 - E_2toI_2 - E_2toA_2;
  du.I_2 = + E_2toI_2 - I_2toR_2;
  du.A_2 = + E_2toA_2 - A_2toR_2;
  du.R_2 = + I_2toR_2 + A_2toR_2 - R_2toRS_2;
  du.RS_2 = + R_2toRS_2 - RS_2toE_2;
      
  du.CE_2 = + S_2toE_2;
  du.CI_2 = + prop_detection * E_2toI_2;
end
