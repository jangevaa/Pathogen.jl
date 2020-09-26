mutable struct MCMC{S <: DiseaseStateSequence, M <: ILM}
  event_observations::EventObservations{S, M}
  event_extents::EventExtents{S}
  population::Population
  starting_states::DiseaseStates
  risk_functions::RiskFunctions{S}
  risk_priors::RiskPriors{S}
  substitution_model::Union{Nothing, Type{<:NucleicAcidSubstitutionModel}}
  substitution_model_prior::Union{Nothing, Vector{UnivariateDistribution}}
  transmission_network_prior::Union{Nothing, TNPrior}
  markov_chains::Vector{MarkovChain{S, M}}

  function MCMC(
    obs::EventObservations{S, M},
    ee::EventExtents{S},
    pop::Population,
    states::DiseaseStates,
    rf::RiskFunctions{S},
    rp::RiskPriors{S};
    tnprior::Union{Nothing, TNPrior}=nothing) where {
    S <: DiseaseStateSequence,
    M <: TNILM}
    return new{S, M}(obs, ee, pop, states, rf, rp, nothing, nothing, tnprior, MarkovChain{S, M}[])
  end

  function MCMC(
    obs::EventObservations{S, M},
    ee::EventExtents{S},
    pop::Population,
    states::DiseaseStates,
    rf::RiskFunctions{S},
    rp::RiskPriors{S},
    sm::Type{N},
    smp::Vector{UnivariateDistribution};
    tnprior::Union{Nothing, TNPrior}=nothing) where {
    S <: DiseaseStateSequence,
    M <: PhyloILM,
    N <: NucleicAcidSubstitutionModel}
    return new{S, M}(obs, ee, pop, states, rf, rp, sm, smp, tnprior, MarkovChain{S, M}[])
  end
end

function MCMC(
  obs::EventObservations{S, M},
  ee::EventExtents{S},
  pop::Population,
  rf::RiskFunctions{S},
  rp::RiskPriors{S};
  tnprior::Union{Nothing, TNPrior}=nothing) where {
  S <: DiseaseStateSequence,
  M <: TNILM}
  return MCMC(obs, ee, pop, fill(State_S, individuals(pop)), rf, rp, tnprior=tnprior)
end

function MCMC(
  obs::EventObservations{S, M},
  ee::EventExtents{S},
  pop::Population,
  rf::RiskFunctions{S},
  rp::RiskPriors{S},
  sm::Type{N},
  smp::Vector{UnivariateDistribution};
  tnprior::Union{Nothing, TNPrior}=nothing) where {
  S <: DiseaseStateSequence,
  M <: PhyloILM,
  N <: NucleicAcidSubstitutionModel}
  return MCMC(obs, ee, pop, fill(State_S, individuals(pop)), rf, rp, sm, smp, tnprior=tnprior)
end

function Base.show(io::IO, x::MCMC{S, M}) where {
                   S <: DiseaseStateSequence,
                   M <: ILM}
  return print(io, "$S $M MCMC with $(length(x.markov_chains)) Markov chains")
end

function TNDistribution(x::MCMC; burnin::Int64=0, thin::Int64=1)
  return TNDistribution([TNDistribution(y, burnin=burnin, thin=thin) for y in x.markov_chains])
end