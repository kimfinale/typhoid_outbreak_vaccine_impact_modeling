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
path_workspace = "G:\\My Drive\\Projects\\VIMC\\VIMC 2.0\\TyphoidOutbreakModel\\";
include(joinpath(path_workspace, "parameters.jl"));
include(joinpath(path_workspace, "seiars_2ag.jl"));
# include(joinpath(path_workspace, "seiarw_2ag_erlang_vacc_two_rounds.jl"));
# include(joinpath(path_workspace, "utils.jl"));
include(joinpath(path_workspace, "utils_temp.jl"));

Random.seed!(42);

dat_ts = DataFrame(CSV.File(joinpath(path_workspace, "data", "outbreak_data_ts.csv"); header=1, delim=","));
dat = DataFrame(CSV.File(joinpath(path_workspace, "data", "outbreak_data_summary.csv"); header=1, delim=","));
# create ID variables to compare with other data set (e.g., time series data)
dat.ID = 1:size(dat,1); 
dat.ID_outbreak = string.(dat.location, "-", dat.start_date, "-", dat.end_date);
dat.attack_rate_naive = dat.total_suspected_cases ./ dat.population;
# data for the outbreaks affected by OCV
dat_ocv = DataFrame(CSV.File(joinpath(path_workspace, "data", "ocv_long_dataset.csv"); header=1, delim=","));
dat_noocv = dat[dat.ID_outbreak .∉ Ref(dat_ocv.ID_outbreak), :]; # n=824
# # include the outbreak occurred in admin 2 or smaller
# dat_noocv_23 = dat_noocv[dat_noocv.spatial_scale .∉ Ref(["admin1"]), :]; # n=766 %in%
# min_delay_vacc = 3; # weeks min = 12 d + 0 d + 4d + 0
# med_delay_vacc = 9; # weeks
# dat_noocv_23_shortdelay = dat_noocv_23[dat_noocv_23.duration .> min_delay_vacc, :]; #n=521
# ids = dat_noocv_23_shortdelay.ID_outbreak;
dat_fit = dat_noocv;
ids = dat_fit.ID_outbreak;
n_ids = length(ids);
# to store output
df = DataFrame(data_id=zeros(n_ids), fit_id=zeros(n_ids), s0=zeros(n_ids), i0=zeros(n_ids), R0=zeros(n_ids),
    n0=zeros(n_ids), objective=zeros(n_ids), retcode=zeros(n_ids)); 
fit_params = Vector{Any}(undef, n_ids);
# # initial conditions for fitting
# x0 = [logit(0.2), logit(0.1), log(2.0)];
# lower = [logit(0.01), logit(0.001), log(1.1)];
# upper = [logit(0.99), logit(0.9), log(20.0)];
# inits = LVector(x0=x0, lower=lower, upper=upper);

x0 = [logit(0.2), logit(0.1), log(2.0), logit(0.5)];
lower = [logit(0.1), logit(0.0001), log(1.1), logit(0.1)];
upper = [logit(0.99), logit(0.9), log(20.0), logit(0.99)];
inits = LVector(x0=x0, lower=lower, upper=upper);

frac_asymp = [0.2];
for fa in frac_asymp 
    @printf("running fA = %.2f\n", fa);    
    # plotdir = string(path_workspace, "\\plots\\4p_fA", round(Int, fa*100), "\\");
    plotdir = @sprintf("%s\\plots\\NB_size10_4p_fA%03d\\", path_workspace, round(Int, fa*100))
    ispath(plotdir) || mkdir(plotdir)
    outputfilename = @sprintf("%s\\outputs\\NB_size10_fit_4p_fA%03d_", path_workspace, round(Int, fa*100))
        
    for i in eachindex(ids)
        @printf("running %d of %d\n", i, length(ids));
        # @printf "running %d of %d\n" i length(ids);
        # params = initialize_params_novacc((loss=nll_4p, fA=0.5, ode=seiarw_2ag_erlang!,
        #     ode_solver = AutoTsit5(Rosenbrock23()), inits = inits));
        params = initialize_params((fA=fa, loss=nll_4p_NB, NB_size=10, ode=seiarw_2ag_erlang_vacc_two_rounds!, inits=inits));    
        # population size of the admin unit in which the outbreak was observed
        # id_pop = dat_noocv_23_shortdelay[dat_noocv_23_shortdelay.ID_outbreak .== ids[i], [:population, :ID]];
        id_pop = dat_fit[dat_fit.ID_outbreak .== ids[i], [:population, :ID]];
        # update the params.u0 based on the pop
        params = merge(params, (population = id_pop.population[1],)); # pop is Vector{Float64}
        params = merge(params, (data_id = id_pop.ID[1],)); # pop is Vector{Float64}
        params = merge(params, (fit_id = i,)); # pop is Vector{Float64}
        # update_params!(params, LVector(population = pop[1]));  
        d = dat_ts[dat_ts.ID_outbreak .== ids[i], [:sCh, :temporal_scale, :TL, :location]];
        params = merge(params, (TL = d.TL,));  
        params = merge(params, (location = d.location[1],));  
        # nobs = size(d, 1);
        if d.temporal_scale[1] == "weekly"
            params = merge(params, (report_freq = 7.0,));
        else
            params = merge(params, (report_freq = 1.0,));
        end
        y = round.(Int, d.sCh); # data
        params = merge(params, (data = y,));   
        params = merge(params, (tend = params.report_freq * length(y),));   
        # prob_ode = ODEProblem(params.ode, params.u0, (0.0, float(params.tend)), params);
        # params = merge(params, (prob_ode = prob_ode,));
        # fit the parameters
        try 
            sol = fit_model(params);
        # store results
            df[i, :data_id] = id_pop.ID[1];
            df[i, :fit_id] = i;
            xhat = extract_xhat(sol);
            params = merge(params, xhat);
            df[i, :s0] = xhat.s0;
            df[i, :i0] = xhat.i0;
            df[i, :R0] = xhat.R0;
            if length(sol.u) > 3
                df[i, :n0] = xhat.n0;
            end    
            df[i, :objective] = sol.objective;
            df[i, :retcode] = Int(sol.retcode); 
            params = merge(params, (sol_optim = sol,));
            x = run_model(params);
            params = merge(params, (modeled = x,));
            fit_params[i] = params;
            # save_plot(params, string(path_workspace, "\\plots\\"));
            save_plot(params, plotdir);
        catch e
            @printf("fitting for i = %d failed\n", i);
        end
    end
    
    dt = Dates.now();
    tstamp = Dates.format(dt, dateformat"yyyymmdd\THHMM");
    jldsave(string(outputfilename, tstamp, ".jld2"); fit_params); # rdata, rds
    jldsave(string(outputfilename, tstamp, ".jld2"); df);
    CSV.write(string(outputfilename, tstamp, ".jld2"), df);
end

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

# sol = solve(prob, Tsit5(), callback = cbs, tstops = tstop)
# using Plots;
# plot(sol);

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
