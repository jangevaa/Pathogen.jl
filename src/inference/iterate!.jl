function _initial_Σrp(mcmc::MCMC{S, M}) where {
                      S <: DiseaseStateSequence,
                      M <: ILM}
  rp = mcmc.risk_priors
  return Float64[i == j ? var(rp[i]) : 0.0 for i in 1:length(rp), j in 1:length(rp)]/(10 * length(rp))
end

function _initial_Σsm(mcmc::MCMC{S, M}) where {
                      S <: DiseaseStateSequence,
                      M <: PhyloILM}
  smp = mcmc.substitution_model_prior
  return Float64[i == j ? var(smp[i]) : 0.0 for i in 1:length(smp), j in 1:length(smp)]/(10 * length(smp))
end

function Base.convert(::Type{Array{Float64, 2}}, sm::Vector{NucleicAcidSubstitutionModel})
  prop_names = [propertynames(sm[1])...]
  return [getproperty(sm[i], θ) for i = 1:length(sm), θ in prop_names]
end

function iterate!(
  mc::MarkovChain{S, M},
  mcmc::MCMC{S, M},
  n::Int64,
  Σrp::Array{Float64, 2},
  σ::Float64;
  condition_on_network::Bool=false,
  event_batches::Int64=1,
  adapt_cov::Int64=100) where {
    S <: DiseaseStateSequence,
    M <: TNILM}
  if adapt_cov < 0
    @warn "`adapt_cov` argument indicates the increment in iterations in which the covariance matrix is updated and must be ≧ 0. Setting to 0 for non-adaptive Metropolis-Hastings."
    adapt_cov = 0
  end
  pmeter = Progress(n, "MCMC progress")
  use_adapted_Σrp = false
  local adapted_Σrp
  if (adapt_cov != 0 & mc.iterations > adapt_cov) && LinearAlgebra.isposdef(value(mc.Σrp))
    use_adapted_Σrp = true
    adapted_Σrp = OnlineStats.value(mc.Σrp) * 2.38^2 / length(mc.risk_parameters[1])
  end
  tr_cache = TransmissionRateCache(individuals(mcmc.population))
  for i = 1:n
    if use_adapted_Σrp
      update!(mc, mcmc, adapted_Σrp, σ, condition_on_network = condition_on_network, event_batches = event_batches, transmission_rate_cache=tr_cache)
    else
      update!(mc, mcmc, Σrp, σ, condition_on_network = condition_on_network, event_batches = event_batches, transmission_rate_cache=tr_cache)
    end
    next!(pmeter)
    @logmsg LogLevel(-5000) "MCMC progress" progress = i/n
    if adapt_cov > 0 && mod(i, adapt_cov) == 0
      OnlineStats.fit!(mc.Σrp, eachrow(convert(Array{Float64, 2}, mc.risk_parameters[end-adapt_cov:end])))
      if LinearAlgebra.isposdef(value(mc.Σrp))
        if !use_adapted_Σrp
          use_adapted_Σrp = true
          @debug "Now using adapted covariance matrix for adaptive MCMC on core $(Distributed.myid()) (updated every $adapt_cov iterations)"
        end
        adapted_Σrp = OnlineStats.value(mc.Σrp) * 2.38^2 / length(mc.risk_parameters[1])
        @debug "Covariance matrix updated for Adaptive Metropolis-Hastings MCMC" Covariance = adapted_Σrp
      end
    end
  end
  return mc
end

function iterate!(
  mc::MarkovChain{S, M},
  mcmc::MCMC{S, M},
  n::Int64,
  Σrp::Array{Float64, 2},
  Σsm::Array{Float64, 2},
  σ::Float64;
  event_batches::Int64=1,
  adapt_cov::Int64=100) where {
    S <: DiseaseStateSequence,
    M <: PhyloILM}
  if adapt_cov < 0
    @warn "`adapt_cov` argument indicates the increment in iterations in which the covariance matrix is updated and must be ≧ 0. Setting to 0 for non-adaptive Metropolis-Hastings."
    adapt_cov = 0
  end
  pmeter = Progress(n, "MCMC progress")
  use_adapted_Σrp = false
  use_adapted_Σsm = false
  local adapted_Σrp
  local adapted_Σsm

  if (adapt_cov != 0 & mc.iterations > adapt_cov)
    if LinearAlgebra.isposdef(value(mc.Σrp))
      use_adapted_Σrp = true
      adapted_Σrp = OnlineStats.value(mc.Σrp) * 2.38^2 / length(mc.risk_parameters[1])
    end
    if LinearAlgebra.isposdef(value(mc.Σsm))
      use_adapted_Σsm = true
      adapted_Σsm = OnlineStats.value(mc.Σsm) * 2.38^2 / length(propertynames(mcmc.substitution_model))
    end
  end
  tr_cache = TransmissionRateCache(individuals(mcmc.population))
  for i = 1:n
    if use_adapted_Σrp & use_adapted_Σsm
      update!(mc, mcmc, adapted_Σrp, adapted_Σsm, σ, event_batches = event_batches, transmission_rate_cache=tr_cache)
    elseif use_adapted_Σrp & !use_adapted_Σsm
      update!(mc, mcmc, adapted_Σrp, Σsm, σ, event_batches = event_batches, transmission_rate_cache=tr_cache)
    elseif !use_adapted_Σrp & use_adapted_Σsm
      update!(mc, mcmc, Σrp, adapted_Σsm, σ, event_batches = event_batches, transmission_rate_cache=tr_cache)
    elseif !use_adapted_Σrp & !use_adapted_Σsm
      update!(mc, mcmc, Σrp, Σsm, σ, event_batches = event_batches, transmission_rate_cache=tr_cache)
    end
    next!(pmeter)
    @logmsg LogLevel(-5000) "MCMC progress" progress = i/n
    if adapt_cov > 0 && mod(i, adapt_cov) == 0
      OnlineStats.fit!(mc.Σrp, eachrow(convert(Array{Float64, 2}, mc.risk_parameters[end-adapt_cov:end])))
      OnlineStats.fit!(mc.Σsm, eachrow(convert(Array{Float64, 2}, mc.substitution_model[end-adapt_cov:end])))
      if LinearAlgebra.isposdef(value(mc.Σrp))
        if !use_adapted_Σrp
          use_adapted_Σrp = true
          @debug "Now using adapted covariance matrix for adaptive MCMC for Risk Function Parameters" core = Distributed.myid() UpdateFreq = adapt_cov
        end
        adapted_Σrp = OnlineStats.value(mc.Σrp) * 2.38^2 / length(mc.risk_parameters[1])
        @debug "Covariance matrix updated for Adaptive Metropolis-Hastings MCMC for Risk Function Parameters" Σrp = adapted_Σrp
      end
      if LinearAlgebra.isposdef(value(mc.Σsm))
        if !use_adapted_Σsm
          use_adapted_Σsm = true
          @debug "Now using adapted covariance matrix for adaptive MCMC on core $(Distributed.myid()) (updated every $adapt_cov iterations)"
        end
        adapted_Σsm = OnlineStats.value(mc.Σsm) * 2.38^2 / length(propertynames(mcmc.substitution_model))
        @debug "Covariance matrix updated for Adaptive Metropolis-Hastings MCMC" Σsm = adapted_Σsm
      end
    end
  end
  return mc
end

function iterate!(mc::MarkovChain{S, M},
                  mcmc::MCMC{S, M},
                  n::Int64,
                  Σrp::Array{Float64, 2},
                  σ::Float64,
                  progress_channel::RemoteChannel;
                  condition_on_network::Bool=false,
                  event_batches::Int64=1,
                  adapt_cov::Int64=100) where {
                  S <: DiseaseStateSequence,
                  M <: TNILM}
  if adapt_cov < 0
    @warn "`adapt_cov` argument indicates the increment in iterations in which the covariance matrix is updated and must be ≧ 0. Setting to 0 for non-adaptive Metropolis-Hastings."
    adapt_cov = 0
  end
  use_adapted_Σrp = false
  local adapted_Σrp
  if (adapt_cov != 0 & mc.iterations > adapt_cov) && LinearAlgebra.isposdef(value(mc.Σrp))
    use_adapted_Σrp = true
    adapted_Σrp = OnlineStats.value(mc.Σrp) * 2.38^2 / length(mc.risk_parameters[1])
  end
  tr_cache = TransmissionRateCache(individuals(mcmc.population))
  for i = 1:n
    if use_adapted_Σrp
      update!(mc, mcmc, adapted_Σrp, σ, condition_on_network = condition_on_network, event_batches = event_batches, transmission_rate_cache=tr_cache)
    else
      update!(mc, mcmc, Σrp, σ, condition_on_network = condition_on_network, event_batches = event_batches, transmission_rate_cache=tr_cache)
    end
    put!(progress_channel, true)
    @logmsg LogLevel(-5000) "MCMC progress" progress = i/n
    if adapt_cov > 0 && mod(i, adapt_cov) == 0
      OnlineStats.fit!(mc.Σrp, eachrow(convert(Array{Float64, 2}, mc.risk_parameters[end-adapt_cov:end])))
      if LinearAlgebra.isposdef(value(mc.Σrp))
        if !use_adapted_Σrp
          use_adapted_Σrp = true
          @debug "Now using adapted covariance matrix for adaptive MCMC on core $(Distributed.myid()) (updated every $adapt_cov iterations)"
        end
        adapted_Σrp = OnlineStats.value(mc.Σrp) * 2.38^2 / length(mc.risk_parameters[1])
        @debug "Covariance matrix updated for Adaptive Metropolis-Hastings MCMC" Covariance = adapted_Σrp
      end
    end
  end
  return mc
end

function iterate!(
  mc::MarkovChain{S, M},
  mcmc::MCMC{S, M},
  n::Int64,
  Σrp::Array{Float64, 2},
  Σsm::Array{Float64, 2},
  σ::Float64,
  progress_channel::RemoteChannel;
  event_batches::Int64=1,
  adapt_cov::Int64=100) where {
    S <: DiseaseStateSequence,
    M <: PhyloILM}
  if adapt_cov < 0
    @warn "`adapt_cov` argument indicates the increment in iterations in which the covariance matrix is updated and must be ≧ 0. Setting to 0 for non-adaptive Metropolis-Hastings."
    adapt_cov = 0
  end
  use_adapted_Σrp = false
  use_adapted_Σsm = false
  local adapted_Σrp
  local adapted_Σsm

  if (adapt_cov != 0 & mc.iterations > adapt_cov)
    if LinearAlgebra.isposdef(value(mc.Σrp))
      use_adapted_Σrp = true
      adapted_Σrp = OnlineStats.value(mc.Σrp) * 2.38^2 / length(mc.risk_parameters[1])
    end
    if LinearAlgebra.isposdef(value(mc.Σsm))
      use_adapted_Σsm = true
      adapted_Σsm = OnlineStats.value(mc.Σsm) * 2.38^2 / length(propertynames(mcmc.substitution_model))
    end
  end
  tr_cache = TransmissionRateCache(individuals(mcmc.population))
  for i = 1:n
    if use_adapted_Σrp & use_adapted_Σsm
      update!(mc, mcmc, adapted_Σrp, adapted_Σsm, σ, event_batches = event_batches, transmission_rate_cache=tr_cache)
    elseif use_adapted_Σrp & !use_adapted_Σsm
      update!(mc, mcmc, adapted_Σrp, Σsm, σ, event_batches = event_batches, transmission_rate_cache=tr_cache)
    elseif !use_adapted_Σrp & use_adapted_Σsm
      update!(mc, mcmc, Σrp, adapted_Σsm, σ, event_batches = event_batches, transmission_rate_cache=tr_cache)
    elseif !use_adapted_Σrp & !use_adapted_Σsm
      update!(mc, mcmc, Σrp, Σsm, σ, event_batches = event_batches, transmission_rate_cache=tr_cache)
    end
    put!(progress_channel, true)
    @logmsg LogLevel(-5000) "MCMC progress" progress = i/n
    if adapt_cov > 0 && mod(i, adapt_cov) == 0
      OnlineStats.fit!(mc.Σrp, eachrow(convert(Array{Float64, 2}, mc.risk_parameters[end-adapt_cov:end])))
      OnlineStats.fit!(mc.Σsm, eachrow(convert(Array{Float64, 2}, mc.substitution_model[end-adapt_cov:end])))
      if LinearAlgebra.isposdef(value(mc.Σrp))
        if !use_adapted_Σrp
          use_adapted_Σrp = true
          @debug "Now using adapted covariance matrix for adaptive MCMC for Risk Function Parameters" core = Distributed.myid() UpdateFreq = adapt_cov
        end
        adapted_Σrp = OnlineStats.value(mc.Σrp) * 2.38^2 / length(mc.risk_parameters[1])
        @debug "Covariance matrix updated for Adaptive Metropolis-Hastings MCMC for Risk Function Parameters" Σrp = adapted_Σrp
      end
      if LinearAlgebra.isposdef(value(mc.Σsm))
        if !use_adapted_Σsm
          use_adapted_Σsm = true
          @debug "Now using adapted covariance matrix for adaptive MCMC on core $(Distributed.myid()) (updated every $adapt_cov iterations)"
        end
        adapted_Σsm = OnlineStats.value(mc.Σsm) * 2.38^2 / length(propertynames(mcmc.substitution_model))
        @debug "Covariance matrix updated for Adaptive Metropolis-Hastings MCMC" Σsm = adapted_Σsm
      end
    end
  end
  return mc
end

function iterate!(
  mcmc::MCMC{S, M},
  n::Int64,
  Σrp::Array{Float64, 2},
  σ::Float64;
  condition_on_network::Bool=false,
  event_batches::Int64=1,
  adapt_cov::Int64=100) where {
    S <: DiseaseStateSequence,
    M <: TNILM}
  if length(mcmc.markov_chains) == 1
    iterate!(mcmc.markov_chains[1], mcmc, n, Σrp, σ,
             condition_on_network = condition_on_network,
             event_batches = event_batches,
             adapt_cov = adapt_cov)
    return mcmc
  else
    if adapt_cov < 0
      @warn "`adapt_cov` argument indicates the increment in iterations in which the covariance matrix is updated and must be ≧ 0. Setting to 0 for non-adaptive Metropolis-Hastings."
      adapt_cov = 0
    end
    pmeter = Progress(n*length(mcmc.markov_chains), "MCMC progress")
    pchannel = RemoteChannel(()->Channel{Bool}(10), 1)
    mc_Futures = Future[]
    @debug "Starting MCMC..."
    for mc in mcmc.markov_chains
      push!(mc_Futures, @spawn iterate!(mc, mcmc, n, Σrp, σ, pchannel, condition_on_network = condition_on_network, event_batches = event_batches, adapt_cov = adapt_cov))
    end
    @debug "MCMC in progress..."
    while !all(isready.(mc_Futures))
      take!(pchannel) && next!(pmeter)
    end
    close(pchannel)
    finish!(pmeter)
    mcmc.markov_chains = [fetch(i) for i in mc_Futures]
    @debug "MCMC complete..."
    return mcmc
  end
end


function iterate!(
  mcmc::MCMC{S, M},
  n::Int64,
  Σrp::Array{Float64, 2},
  Σsm::Array{Float64, 2},
  σ::Float64;
  condition_on_network::Bool=false,
  event_batches::Int64=1,
  adapt_cov::Int64=100) where {
    S <: DiseaseStateSequence,
    M <: PhyloILM}
  if length(mcmc.markov_chains) == 1
    iterate!(mcmc.markov_chains[1], mcmc, n, Σrp, Σsm, σ,
             event_batches = event_batches,
             adapt_cov = adapt_cov)
    return mcmc
  else
    if adapt_cov < 0
      @warn "`adapt_cov` argument indicates the increment in iterations in which the covariance matrix is updated and must be ≧ 0. Setting to 0 for non-adaptive Metropolis-Hastings."
      adapt_cov = 0
    end
    pmeter = Progress(n*length(mcmc.markov_chains), "MCMC progress")
    pchannel = RemoteChannel(()->Channel{Bool}(10), 1)
    mc_Futures = Future[]
    @debug "Starting MCMC..."
    for mc in mcmc.markov_chains
      push!(mc_Futures, @spawn iterate!(mc, mcmc, n, Σrp, Σsm, σ, pchannel, event_batches = event_batches, adapt_cov = adapt_cov))
    end
    @debug "MCMC in progress..."
    while !all(isready.(mc_Futures))
      take!(pchannel) && next!(pmeter)
    end
    close(pchannel)
    finish!(pmeter)
    mcmc.markov_chains = [fetch(i) for i in mc_Futures]
    @debug "MCMC complete..."
    return mcmc
  end
end

function iterate!(
  mc::MarkovChain{S, M},
  mcmc::MCMC{S, M},
  n::Int64,
  σ::Float64;
  condition_on_network::Bool=false,
  event_batches::Int64=1,
  adapt_cov::Int64=100) where {
    S <: DiseaseStateSequence,
    M <: TNILM}
  return iterate!(
    mc,
    mcmc,
    n,
    _initial_Σrp(mcmc),
    σ,
    condition_on_network = condition_on_network,
    event_batches = event_batches,
    adapt_cov = adapt_cov)
end

function iterate!(
  mc::MarkovChain{S, M},
  mcmc::MCMC{S, M},
  n::Int64,
  σ::Float64;
  event_batches::Int64=1,
  adapt_cov::Int64=100,
  DA_couple::Bool=false) where {
    S <: DiseaseStateSequence,
    M <: PhyloILM}
  return iterate!(
    mc,
    mcmc,
    n,
    _initial_Σrp(mcmc),
    _initial_Σsm(mcmc),
    σ,
    event_batches = event_batches,
    adapt_cov = adapt_cov,
    DA_couple = DA_couple)
end

function iterate!(
  mcmc::MCMC{S, M},
  n::Int64,
  σ::Float64;
  condition_on_network::Bool=false,
  event_batches::Int64=1,
  adapt_cov::Int64=100) where {
    S <: DiseaseStateSequence,
    M <: TNILM}
  return iterate!(
    mcmc,
    n,
    _initial_Σrp(mcmc),
    σ,
    condition_on_network = condition_on_network,
    event_batches = event_batches,
    adapt_cov = adapt_cov)
end

function iterate!(
  mcmc::MCMC{S, M},
  n::Int64,
  σ::Float64;
  event_batches::Int64=1,
  adapt_cov::Int64=100,
  DA_couple::Bool=false) where {
    S <: DiseaseStateSequence,
    M <: PhyloILM}
  return iterate!(
    mcmc,
    n,
    _initial_Σrp(mcmc),
    _initial_Σsm(mcmc),
    σ,
    event_batches = event_batches,
    adapt_cov = adapt_cov)
end