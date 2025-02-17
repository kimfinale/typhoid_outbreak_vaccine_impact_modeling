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
path_workspace = "G:\\My Drive\\Projects\\VIMC\\VIMC 2.0\\CholeraOutbreakModel\\julia";
include(joinpath(path_workspace, "parameters.jl"));
include(joinpath(path_workspace, "seiarw_2ag_erlang_vacc.jl"));
include(joinpath(path_workspace, "seiarw_2ag_erlang_vacc_two_rounds.jl"));
include(joinpath(path_workspace, "utils.jl"));

Random.seed!(42);
# fitted parameter
f = jldopen(string(path_workspace, "\\outputs\\fit_4p_20240207T2155.jld2"), "r");
fit_params = f["fit_params"];
p = fit_params[1];
# 
u0 = initialize_u0();
params = initialize_params();
# delete(nt, :u0)
# params = merge(params, p); that would be ideal, but I need to do as below
params = merge(params, (s0=p.s0, i0=p.i0, R0=p.R0, n0=p.n0, tend=p.tend, report_freq=p.report_freq,
TL=p.TL, data_id=p.data_id, fit_id=p.fit_id, data=p.data, location=p.location, population=p.population));
params = merge(params, (ode=seiarw_2ag_erlang_vacc_two_rounds!, u0=u0));

params = adjust_pop_dist!(params);
prob = ODEProblem(params.ode, params.u0, (0.0, float(params.tend)), params);
sol = solve(prob, params.ode_solver, saveat = params.report_freq);

df = DataFrame(sol);
rename!(df, [:t; collect(keys(params.u0))]);

ids_cumul = findall(x -> x ∈ [:CE_1, :CI_1, :CE_2, :CI_2], keys(params.u0))
ids = deleteat!(collect(2:73), ids_cumul);
# ks = deleteat!(keys(params.u0), ids)
params.population * params.n0
sum(df[1, ids])
round(params.population*params.n0, sigdigits=3) == round(sum(df[1, ids]), sigdigits=3)
round(sum(df[17, ids]), sigdigits=3) == round(sum(df[1, ids]), sigdigits=3)


# x = run_model(params);
# params = merge(params, (modeled = x,));

# plt = plot_incidence(params);
# plt

# params = adjust_pop_dist!(params);
# prob = ODEProblem(params.ode, params.u0, (0.0, float(params.tend)), params);
# sol = solve(prob, params.ode_solver, callback=params.callback, tstops=params.vacc_1d_times, saveat=params.report_freq);
# df = DataFrame(sol);
# rename!(df, [:t; collect(keys(params.u0))]);
# size(df)

# params = adjust_vacc_rate!(params);
# params = merge(params, (vacc_1d_eff_1=0.2, vacc_1d_eff_2=0.6, vacc_1d_cov=0.2,
#     vacc_1d_immunity_rate= 1/7.0,
#     campaign_1d_start=28.0, campaign_1d_dur=7.0));
# params = merge(params, (vacc_1d_times = LVector(start = params.campaign_1d_start, 
# stop = params.campaign_1d_start + params.campaign_1d_dur),));

# params = set_vaccination_callback(params);
# params = adjust_pop_dist!(params);
# prob_v = ODEProblem(params.ode, params.u0, (0.0, float(params.tend)), params);
# sol_v = solve(prob_v, params.ode_solver, callback=params.callback, tstops=params.vacc_1d_times, saveat=params.report_freq);
# df_v = DataFrame(sol_v);
# rename!(df_v, [:t; collect(keys(params.u0))]);
# sum(df_v[1, ids])
# round(sum(df_v[1, ids]), sigdigits=3) == round(sum(df[1, ids]), sigdigits=3)
# df_v[1:10, 1:22]
# df_v[end-10:end, 1:20]
# df_v[end-10:end, 21:40]
# df_v[end-10:end, 41:60]
# df_v[end-10:end, 61:73]
# x_v = run_model_vacc(params);
# params = merge(params, (modeled_vacc = x_v,));

# plot!(params.TL, params.modeled_vacc, color=:green);
# plt

params = merge(params, (vacc_1d_eff_1=0.2, vacc_1d_eff_2=0.6, vacc_1d_cov=0.8,
    vacc_1d_immunity_rate= 1/7.0, campaign_1d_start=28.0, campaign_1d_dur=7.0));
params = merge(params, (vacc_1d_times = LVector(start = params.campaign_1d_start, 
    stop = params.campaign_1d_start + params.campaign_1d_dur),));
params = merge(params, (report_freq=1.0,));
# prob_v = ODEProblem(params.ode, params.u0, (0.0, float(params.tend)), params);
# sol_v = solve(prob_v, params.ode_solver, callback=params.callback, tstops=params.vacc_1d_times, saveat=params.report_freq);
# df_v = DataFrame(sol_v);
# rename!(df_v, [:t; collect(keys(params.u0))]);
# sum(df_v[1, ids])
# round(sum(df_v[1, ids]), sigdigits=3) == round(sum(df[1, ids]), sigdigits=3)
# df_v[1:10, 1:22]
# df_v[end-10:end, 1:20]
# df_v[end-10:end, 21:40]
# df_v[end-10:end, 41:60]
# df_v[end-10:end, 61:73]

# params = merge(params, (ode_solver = AutoVern7(Rodas4P()),));
x = run_model(params);
params = merge(params, (modeled = x,));
x_v = run_model_vacc(params);
params = merge(params, (modeled_vacc = x_v,));
# params = merge(params, eval(Meta.parse(string("(inc_v", x, " = x_v,)"))));
# learn to use unit testing
plt = plot(x, color=:black);
plot!(x_v, color=:tomato)
plt

##---------------------------------------------------------------------------------------------
# vaccine impact simulation
vacc_start_weeks = 3.0:7.0;
vacc_start_color = cgrad(:matter, 5, categorical = true);
# vaccine coverage 0.2 - 0.9, vaccine start from 3 weeks to 
pathdir = string(path_workspace, "\\plots\\vacc\\");
vacc_fit_params = Vector{Any}(undef, length(fit_params));
# undefined elements replaced with 0
idx = [52,91,142,143,649]
for x in idx
    # fit_params[x] = 0; 
    vacc_fit_params[x] = 0;  
end 
fit_params[52] == 0

for i in 1:length(fit_params)
    @printf("running %d of %d\n", i, length(fit_params));
    p = fit_params[i]; # extract fitted parameter values
    if p == 0
        @printf("empty element\n"); 
        vacc_fit_params[i] = 0;       
    else
        u0 = initialize_u0();
        params = initialize_params();
        params = merge(params, (s0=p.s0, i0=p.i0, R0=p.R0, n0=p.n0, tend=p.tend, report_freq=p.report_freq,
            TL=p.TL, data_id=p.data_id, fit_id=p.fit_id, data=p.data, location=p.location, population=p.population));
        params = merge(params, (ode=seiarw_2ag_erlang_vacc_two_rounds!, u0=u0, report_freq=1.0)); # vaccine model, new 
        inc_novacc = run_model(params);
        params = merge(params, (inc_novacc = inc_novacc,));
        plt = plot(inc_novacc, label="no vacc", linewidth=2, linecolor=:black,
            title=string("ID: ", p.data_id, ", loc: ", p.location), xlabel = "Day", ylabel = "Incidence");
        for j in 1:length(vacc_start_weeks)
            vw = vacc_start_weeks[j] * 7.0;
            if (vw + params.campaign_1d_dur) < params.tend  
                params = merge(params, (vacc_1d_eff_1=0.2, vacc_1d_eff_2=0.6, vacc_1d_cov=0.8,
                    vacc_1d_immunity_rate=1/7.0, campaign_1d_start=vw, campaign_1d_dur=7.0));
                params = merge(params, (vacc_1d_times = LVector(start = params.campaign_1d_start, 
                stop = params.campaign_1d_start + params.campaign_1d_dur),));
                inc_v = run_model_vacc(params);
                vw_int = round(Int, vw);
                # C[1-9]T[1-9] represent the coverage proportion (*10) and start timing (weeks) for the vaccination campaign 
                params = merge(params, eval(Meta.parse(string("(inc_C8T", vw_int ÷ 7, " = inc_v,)"))));
                
                plot!(inc_v, color=vacc_start_color[j]);
            end
        end
        vacc_fit_params[i] = params;
        dt = Dates.now();
        tstamp = Dates.format(dt, dateformat"yyyymmdd\THHMM");
        png(plt, string(pathdir, "fig_", p.fit_id, "_id_", p.data_id, "_", tstamp, ".png"));
    end
end

using Test

jldsave(string(path_workspace, "\\outputs\\vacc_4p_", tstamp, ".jld2"); vacc_fit_params);
# # sanity check
# p = initialize_params();
# p = merge(p, (R0=2.0, ode = seiarw_2ag_erlang_vacc_two_rounds!, sigma=0.0, tend=1000.0));

# prob = ODEProblem(p.ode, p.u0, (0.0, float(p.tend)), p);
# sol = solve(prob, p.ode_solver, saveat = p.report_freq);
# ci = sol[end].CI_1 + sol[end].CI_2
# cr = sol[end].R1_1 + sol[end].R1_2
# cr/p.population
# final_epidemic_size(p.R0)

# paste("outbreak size (seirw_R) =", sum(tail(out[,grep("R", names(out))],1))/PARAMS$population)
# sum(head(out[,grep("S|^E|^I|R", names(out))],1))
# sum(tail(out[,grep("S|^E|^I|R", names(out))],1)) # ^E|^I to exclude CE and CI

# ## vaccination check
# params = initialize_params();
# params = merge(params, (tend = 60.0,));
# params = merge(params, (ode = seiarw_2ag_erlang!,));
# params = merge(params, (ode_solver = AutoTsit5(Rosenbrock23()),)); 
# params = merge(params, (campaign_start = 50.0,));
# prob_ode = ODEProblem(params.ode, params.u0, (0.0, float(params.tend)), params);
# params = merge(params, (prob_ode = prob_ode,));

# # params = adjust_vacc_rate!(params);
# params = merge(params, (campaign_start=38.0,));
# params = merge(params, (campaign_dur=18.0,));
# params = merge(params, (vaccination_times = LVector(start = params.campaign_start, 
# stop = params.campaign_start + params.campaign_dur),));
# # params = merge(params, (vaccination_times = vaccination_times,));

# params = set_vaccination_callback(params);
# # params = merge(params, (cbs = cbs,));
# prob = remake(params.prob_ode; u0=params.u0, p=params); 
# sol = solve(prob, params.ode_solver, callback=params.callback, 
#     tstops = params.vaccination_times, saveat=params.report_freq);

# condition(u,t,integrator) = t ∈ [0.0, 100.0]s
# condition(u,t,integrator) = t < vaccination_times[1] # condition(u,t,integrator) = t < params.vaccination_times[1] || t > params.vaccination_times[2]
# function affect!(integrator)
#     if integrator.t < vaccination_times[1] || integrator.t > vaccination_times[2] 
#         integrator.p.vacc_rate.first = 0.0;
#     else
#         integrator.p.vacc_rate.first = 0.1;
#     end
# end

# cb = PresetTimeCallback(vaccination_times, affect!);
# dosetimes = [4.0, 8.0]
# affect!(integrator) = integrator.u[1] += 10
# cb = PresetTimeCallback(dosetimes, affect!)

# cb = DiscreteCallback(condition, affect!);
# p = (α=1.0, β=2.0);
# function change1!(p)
#     p = merge(p, (β=0.2,));
# end
# change1!(p)
# print(p)

# l = LVector(α=1.0, β=2.0);
# function change2!(l)
#     l.β = 0.2; 
# end
# change2!(l)
# print(l)

# function change3(p)
#     p = merge(p, (β=0.2,));
#     return p
# end    
# p = change3(p);
# print(p)

 
# Keep the codes for Bayesian fitting
# # Bayesian parameter estimation using Turing package
# using Turing

# @model bayesian_fit(y) = begin  
#     p ~ Uniform(0.0, 10.0)
#     # params.R0 = exp(p[1]);
#     params.R0 = exp(p);
#     # @show params.R0;
#     tend = 10;
#     # type conversion required for autodiff (e.g., to use NUTS)
#     prob = remake(prob_ode; u0=convert.(eltype(p), prob_ode.u0), p = params)
#     # prob = remake(prob_ode; p = params)
#     # prob = ODEProblem(seiarw_2ag_erlang!, u0, (0.0, float(tend)), params);
#     # @show prob
#     sol = solve(prob, Tsit5(), saveat=1);
#     X = [sol[i].CI_1 for i in 2:round(Int, tend+1)] - [sol[i].CI_1 for i in 1:round(Int, tend)];

#     for i in 1:length(y)
#       y[i] ~ Poisson(X[i])
#     end
#   end;


# @time chain_nuts = sample(bayesian_fit(Y), NUTS(0.65), 10_000);
# histogram(chain_nuts)
# describe(chain_nuts)

# chaindf = DataFrame(chain_nuts)
# plot(chaindf.iteration, chaindf.p)

# using DynamicHMC
# dynamic_nuts = externalsampler(DynamicHMC.NUTS())
# @time chain_dnuts = sample(bayesian_fit(Y), dynamic_nuts, 10_000)

# keep the codes for initial modeling efforts
# # initial values for the state variables
# params = initialize_params((R0=3.5,));
# true_parms = (s0=0.9, i0=0.001, R0=2.9);
# params = merge(params, true_parms);
# params = merge(params, (population=1e5,));
# adjust_pop_dist!(params);
# # adjust the initial population

# prob_ode = ODEProblem(seiarw_2ag_erlang!, params.u0, (0.0, params.tend), params);
# prob = remake(prob_ode, p=params);
# params = merge(params, (prob_ode=prob_ode,));

# @time sol = solve(prob, Tsit5(), saveat=params.report_freq);
# @time sol = solve(prob, AutoTsit5(Rosenbrock23()), saveat=params.report_freq);

# # using BenchmarkTools
# # @btime sol = solve(prob, Tsit5(), saveat=params.report_freq);
# # @btime sol = solve(prob, AutoTsit5(Rosenbrock23()), saveat=params.report_freq);

# plot(sol)
# plot(sol.t, [sol[i].CI_1 for i in eachindex(sol.t)])
# plot(sol.t, [sol[i].I1_1 for i in eachindex(sol.t)])

# incid = [sol[i].CI_1 + sol[i].CI_2 for i in 2:length(sol.t)] - [sol[i].CI_1 + sol[i].CI_2 for i in 1:(length(sol.t)-1)];
# for i in eachindex(incid)
#     if incid[i] < 1e-6
#         incid[i] = 1e-6
#     end
# end
# # Assume y_t \distributedas Y_t
# y = rand.(Poisson.(incid));
# params = merge(params, (data=y,));


# x0 = [logit(0.2), logit(0.1), log(2.0)];
# collect(true_parms)

# # BlackBox optimization
# lower = [logit(0.01), logit(0.001), log(1.1)];
# upper = [logit(0.99), logit(0.9), log(20.0)];

# prob = OptimizationProblem(nll_3p, x0, params, lb=lower, ub=upper);
# @time sol = solve(prob, BBO_adaptive_de_rand_1_bin_radiuslimited());
# @time sol = solve(prob, BBO_adaptive_de_rand_1_bin_radiuslimited(), NumDimensions=3, PopulationSize=500);
# print_xhat(sol.minimizer)
# # @printf("s0 = %.4f, i0 = %.4f, R0 = %.4f\n", expit(sol[1]), expit(sol[2]), exp(sol[3]));
# printf_sol(sol);


# x0 = [logit(0.2), logit(0.1), log(2.0)];
# lower = [logit(0.01), logit(0.001), log(1.1)];
# upper = [logit(0.99), logit(0.9), log(20.0)];
# inits = LVector(x0=x0, lower=lower, upper=upper);

# params = merge(params, (inits=inits,));
