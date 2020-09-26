function generate(
  ::Type{Transmission},
  tr::TransmissionRates,
  tnd::Nothing,
  id::Int64)
  external_or_internal = Weights([tr.external[id];
                                  sum(tr.internal[:, id])])
  if sum(external_or_internal) == 0.0
    @error "All transmission rates = 0.0, No transmission generated" external_or_internal
    return NoTransmission()
  elseif sample([true; false], external_or_internal)
    @debug "Exogenous tranmission generated" external_or_internal
    return ExogenousTransmission(id)
  else
    source = sample(1:individuals(tr), Weights(tr.internal[:, id]))
    @debug "Endogenous transmission generated (source id = $source)" external_or_internal
    return EndogenousTransmission(id, source)
  end
end

function generate(
  ::Type{Transmission},
  tr::TransmissionRates,
  tnd::TNDistribution,
  id::Int64)
  internalweights = tr.internal[:, id] .* tnd.internal[:, id]
  external_or_internal = Weights([tr.external[id] * tnd.external[id];
                                 sum(internalweights)])
  if sum(external_or_internal) == 0.0
    @error "All transmission rates = 0.0, No transmission generated" external_or_internal
    return NoTransmission()
  elseif sample([true; false], external_or_internal)
    @debug "Exogenous tranmission generated" external_or_internal
    return ExogenousTransmission(id)
  else
    source = sample(1:individuals(tr), Weights(internalweights))
    @debug "Endogenous transmission generated (source id = $source)" external_or_internal
    return EndogenousTransmission(id, source)
  end
end

function generate(
  ::Type{TransmissionNetwork},
  starting_states::DiseaseStates,
  tr::TransmissionRates,
  tnd::Union{Nothing, TNDistribution},
  ids::Vector{Int64})
  tn = TransmissionNetwork(starting_states)
  for i in ids
    tx = generate(Transmission, tr, tnd, i)
    update!(tn, tx)
  end
  return tn
end

function generate(
  ::Type{TransmissionNetwork},
  tr::TransmissionRates,
  mcmc::MCMC,
  ids::Vector{Int64})
  return generate(TransmissionNetwork, mcmc.starting_states, tr, mcmc.transmission_network_prior, ids)
end

function generate(
  ::Type{TransmissionNetwork},
  lasttn::TransmissionNetwork,
  tr::TransmissionRates,
  tnd::Union{Nothing, TNDistribution},
  ids)
  tn = copy(lasttn)
  for i in ids
    tn.external[i] = 0
    tn.internal[:, i] .= 0
    tx = generate(Transmission, tr, tnd, i)
    update!(tn, tx)
  end
  return tn
end

function generate(
  ::Type{TransmissionNetwork},
  lasttn::TransmissionNetwork,
  tr::TransmissionRates,
  mcmc::MCMC,
  ids)
  return generate(TransmissionNetwork, lasttn, tr, mcmc.transmission_network_prior, ids)
end

function generate(
  ::Type{Events},
  obs::EventObservations{T, M},
  extents::EventExtents{T}) where {
  T <: DiseaseStateSequence,
  M <: ILM}
  events = Events{T}(individuals(obs))
  exposed_state = State_E ∈ T
  removed_state = State_R ∈ T
  for i = 1:individuals(obs)
    if obs.infection[i] == -Inf
      update!(events, Event{T}(-Inf, i, State_I))
      if exposed_state
        update!(events, Event{T}(-Inf, i, State_E))
      end
      if removed_state && obs.removal[i] == -Inf
        update!(events, Event{T}(-Inf, i, State_R))
      end
    elseif !isnan(obs.infection[i])
      i_lb = obs.infection[i] - extents.infection[2]
      i_ub = obs.infection[i] - extents.infection[1]
      i_time = rand(Uniform(i_lb, i_ub))
      update!(events, Event{T}(i_time, i, State_I))
      if exposed_state
        e_lb = i_time - extents.exposure[2]
        e_ub = i_time - extents.exposure[1]
        e_time = rand(Uniform(e_lb, e_ub))
        update!(events, Event{T}(e_time, i, State_E))
      end
      if removed_state && !isnan(obs.removal[i])
        r_lb = maximum([obs.infection[i]; obs.removal[i] - extents.removal[2]])
        r_ub = obs.removal[i] - extents.removal[1]
        r_time = rand(Uniform(r_lb, r_ub))
        update!(events, Event{T}(r_time, i, State_R))
      end
    end
  end
  return events
end

function generate(
  ::Type{Events},
  mcmc::MCMC{T}) where {
  T <: DiseaseStateSequence}
  return generate(Events, mcmc.event_observations, mcmc.event_extents)
end

function generate(
  ::Type{Event},
  last_event::Event{T},
  σ::Float64,
  extents::EventExtents{T},
  obs::EventObservations,
  events::Events{T}) where {
  T <: DiseaseStateSequence}
  lowerbound, upperbound = _bounds(last_event, extents, obs, events)
  time = rand(truncated(Normal(_time(last_event), σ),
                        lowerbound,
                        upperbound))
  return Event(time, last_event)
end

function generate(
  ::Type{Event},
  last_event::Event{T},
  σ::Float64,
  extents::EventExtents{T},
  obs::EventObservations,
  events::Events{T},
  network::TransmissionNetwork) where {
  T <: DiseaseStateSequence}
  lowerbound, upperbound = _bounds(last_event, extents, obs, events, network)
  time = rand(truncated(Normal(_time(last_event), σ),
                        lowerbound,
                        upperbound))
  return Event(time, last_event)
end

function generate(
  ::Type{RiskParameters{T}},
  rpriors::RiskPriors{T}) where {
  T <: DiseaseStateSequence}
  sparks = Float64[rand(x) for x in rpriors.sparks]
  susceptibility = Float64[rand(x) for x in rpriors.susceptibility]
  infectivity = Float64[rand(x) for x in rpriors.infectivity]
  transmissibility = Float64[rand(x) for x in rpriors.transmissibility]
  if State_E ∈ T
    latency = Float64[rand(x) for x in rpriors.latency]
  end
  if State_R ∈ T
    removal = Float64[rand(x) for x in rpriors.removal]
  end
  if T == SEIR
    return RiskParameters{T}(sparks,
                             susceptibility,
                             infectivity,
                             transmissibility,
                             latency,
                             removal)
  elseif T == SEI
    return RiskParameters{T}(sparks,
                             susceptibility,
                             infectivity,
                             transmissibility,
                             latency)
  elseif T == SIR
    return RiskParameters{T}(sparks,
                             susceptibility,
                             infectivity,
                             transmissibility,
                             removal)
  elseif T == SI
    return RiskParameters{T}(sparks,
                             susceptibility,
                             infectivity,
                             transmissibility)
  end
end

function generate(::Type{RiskParameters}, mcmc::MCMC{T}) where T <: DiseaseStateSequence
  return generate(RiskParameters{T}, mcmc.risk_priors)
end

function generate(
  ::Type{RiskParameters{T}},
  last_rparams::RiskParameters{T},
  Σ::Array{Float64, 2}) where {
  T <: DiseaseStateSequence}
  rparams_vector = rand(MvNormal(convert(Vector{Float64}, last_rparams), Σ))
  return _like(last_rparams, rparams_vector)
end

function generate(
  ::Type{RiskParameters{T}},
  mc::MarkovChain{T},
  Σ::Array{Float64, 2}) where {
  T <: DiseaseStateSequence}
  return generate(RiskParameters{T}, mc.risk_parameters[end], Σ)
end

function generate(
  ::Type{SM},
  priors::Vector{UnivariateDistribution}) where {
  SM <: NucleicAcidSubstitutionModel}
  return SM([rand(x) for x in priors])
end

function generate(
  ::Type{SM},
  last_sm::SM,
  Σ::Array{Float64, 2}) where {
  SM <: NucleicAcidSubstitutionModel}
  if size(Σ) == (0,0)
    return last_sm
  else
    return SM(rand(MvNormal([getproperty(last_sm, θ) for θ in [fieldnames(SM)...]], Σ)), false)
  end
end