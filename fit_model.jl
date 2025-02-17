# using Revise
using DifferentialEquations
using Plots
using CSV
using DataFrames
using Random
using LabelledArrays
using Distributions
using NamedTupleTools
using Dates
using Printf
using Optimization
using ForwardDiff
using OptimizationOptimJL
using OptimizationBBO
using JLD2
# load model function, parameter LabelledArrays
# path_juliafiles = "G:\\My Drive\\Projects\\VIMC\\VIMC 2.0\\CholeraOutbreakModel\\julia\\";
path_workspace = "G:\\My Drive\\Projects\\VIMC\\VIMC 2.0\\TyphoidOutbreakModel";
include(joinpath(path_workspace, "parameters.jl"));
include(joinpath(path_workspace, "seiars_2ag.jl"));
# include(joinpath(path_workspace, "seiarw_2ag_erlang_vacc_two_rounds.jl"));
include(joinpath(path_workspace, "utils.jl"));

# check the ODE model
parm = initialize_params_seiars((ode=seiars_2ag!, tend=200.0, ode_solver=Tsit5(), population=1e3, report_freq=7));
parm = adjust_pop_dist_seiars!(parm);
prob = ODEProblem(parm.ode, parm.u0, (0.0, float(parm.tend)), parm); 
sol = solve(prob, parm.ode_solver, saveat = Int(parm.report_freq)); 
plot(sol)
# @show sol
X = [sol[i].CI_1 + sol[i].CI_2 for i in 2:length(sol.t)] - [sol[i].CI_1 + sol[i].CI_2 for i in 1:(length(sol.t)-1)];
plot(X)               
dat = rand.(Poisson.(X));

inits = LVector(x0=[log(2)], upper=[log(20)], lower=[log(0.99)]);
parm = initialize_params_seiars((ode=seiars_2ag!, tend=200.0, ode_solver=Tsit5(), population=1e3, report_freq=7));
# params for fitting
parm = merge(parm, (loss=nll_1p, inits=inits, data=dat, de_params=LVector(nd=1, pop=500)));
prob_opt = OptimizationProblem(parm.loss, parm.inits.x0, parm, lb=parm.inits.lower, ub=parm.inits.upper);
sol = solve(prob_opt, BBO_adaptive_de_rand_1_bin_radiuslimited(), 
    NumDimensions=parm.de_params.nd, PopulationSize=parm.de_params.pop);
exp(sol.u[1]) 
parm.R0

