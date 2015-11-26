using Pathogen, DataFrames, Distributions, ProgressMeter

# Simulate
init_seq = create_seq(100, 0.25, 0.25, 0.25, 0.25)
init_var = rand(Uniform(0,25), (2,100))

pop = create_population(init_seq, init_var)

powerlaw = create_powerlaw(3., 4., 0.001)
latency = create_constantrate(Inf)
recovery = create_constantrate(1/7.)
substitution = jc69q([0.001])

ratearray = create_ratearray(pop, powerlaw, substitution)

while length(pop.timeline[1]) < 300
  onestep!(ratearray, pop, powerlaw, latency, recovery, substitution)
end

# Inference
# α, β: powerlaw exposure kernel parameters
# η: external pressure rate
# γ: recovery rate (1/mean infectious period)
# λ: JC69 transition/transversion rate

actual, obs = surveil(pop)

ilm_priors = SIR_priors(Gamma(3.),
                        Gamma(4.),
                        Uniform(0., 0.002),
                        Gamma(1/7))

mutation_priors = JC69_priors(Uniform(0., 0.002))

ilm_trace, mutation_trace = MCMC(200000, ilm_priors, mutation_priors, obs)

# # Tune the transition kernel's covariance matrix
# n = 300
# progressbar = Progress(n, 5, "Performing $n tuning MCMC stages...", 25)
# for i = 1:n
#   # Tune transition matrix
#   opt_cov = cov([ilm_trace.α ilm_trace.β ilm_trace.γ ilm_trace.η mutation_trace.λ])*(2.38^2)/5.
#
#   # Perform 1000 MCMC iterations
#   MCMC(1000, opt_cov, ilm_trace, mutation_trace, ilm_priors, mutation_priors, obs, false, false)
#   next!(progressbar)
# end
#
# opt_cov = cov([ilm_trace.α ilm_trace.β ilm_trace.γ ilm_trace.η mutation_trace.λ])*(2.38^2)/5.
# MCMC(100000, opt_cov, ilm_trace, ilm_priors, obs)

using Gadfly, DataFrames
cd("/Users/justin/Desktop/pathogen")

# Simulation/Maximum posteriori visualization
images = 500
max_tracelp=findfirst(ilm_trace.logposterior.==maximum(ilm_trace.logposterior))

for time = 1:images
  states, routes = plotdata(actual,
                            pop,
                            (time*maximum([maximum(ilm_trace.aug[max_tracelp]), maximum(obs)])/images))
  p1 = plot(layer(states, x="x", y="y", color="state", Geom.point),
            layer(routes, x="x", y="y", group="line", Geom.polygon),
            Theme(panel_opacity=1.,
                  panel_fill=colorant"white",
                  default_color=colorant"black",
                  background_color=colorant"white"))

  states, routes = plotdata(obs,
                            ilm_trace,
                            max_tracelp,
                            (time*maximum([maximum(ilm_trace.aug[max_tracelp]), pop.timeline[1][end]])/images))
  p2 = plot(layer(states, x="x", y="y", color="state", Geom.point),
            layer(routes, x="x", y="y", group="line", Geom.polygon),
            Theme(panel_opacity=1.,
                  panel_fill=colorant"white",
                  default_color=colorant"black",
                  background_color=colorant"white"))

  filenumber = time/images
  filenumber = prod(split("$filenumber", ".", limit=2))
  filenumber *= prod(fill("0", 5-length(filenumber)))
  draw(PNG("SIR_simulation_$filenumber.png", 15cm, 20cm), vstack(p1,p2))
end

# Assemble into animation
run(`convert -delay 10 -loop 0 -layers optimize SIR_simulation_*.png SIR3_animation_combined.gif`)

# Remove frames
for time = 1:images
  filenumber = time/images
  filenumber = prod(split("$filenumber", ".", limit=2))
  filenumber *= prod(fill("0", 5-length(filenumber)))
  rm("SIR_simulation_$filenumber.png")
end


# Inference visualization
# Joint trace plots (last 100k iterations)
plotdf = DataFrame(iteration = rep(1:200000, 5),
                   value = [ilm_trace.α[end-199999:end];
                            ilm_trace.β[end-199999:end];
                            ilm_trace.η[end-199999:end];
                            ilm_trace.γ[end-199999:end];
                            mutation_trace.λ[end-199999:end]],
                   parameter = [rep("alpha", 200000);
                                rep("beta", 200000);
                                rep("eta", 200000);
                                rep("gamma", 200000);
                                rep("lambda", 200000)])

draw(PNG("SIR3_traceplot.png", 20cm, 15cm),
     plot(plotdf,
          x="iteration",
          y="value",
          color="parameter",
          Geom.line,
          Theme(panel_opacity=1.,
                panel_fill=colorant"white",
                background_color=colorant"white")))


# logposterior plot (last 100k iterations)
draw(PNG("SIR3_logposterior.png", 20cm, 15cm),
     plot(x=1:200000,
          y=ilm_trace.logposterior[end-199999:end],
          Geom.line,
          Theme(panel_opacity=1.,
                panel_fill=colorant"white",
                background_color=colorant"white")))
