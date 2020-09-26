function _initialization_attempt(mcmc::MCMC{S, M},
                                 max_lposterior::Float64) where {
                                 S <: DiseaseStateSequence,
                                 M <: ILM}
  events = generate(Events, mcmc)
  rparams = generate(RiskParameters, mcmc)
  lprior = logprior(rparams, mcmc.risk_priors)
  if M <: PhyloILM
    sm = generate(mcmc.substitution_model, mcmc.substitution_model_prior)
    lprior += logprior(sm, mcmc.substitution_model_prior)
  end
  llikelihood, tr, tx = loglikelihood(rparams,
                                      mcmc.risk_functions,
                                      events,
                                      mcmc.population,
                                      mcmc.starting_states,
                                      early_decision_value = max_lposterior - lprior)
  lposterior = llikelihood + lprior
  if lposterior > max_lposterior
    network = generate(TransmissionNetwork, tr, mcmc, tx)
    if M <: PhyloILM
      tree, obs_nodes = generate(Tree, events, mcmc.event_observations, network)
      leaf_data = Dict{Int64, RNASeq}()
      for k = eachindex(obs_nodes)
        if !isnothing(obs_nodes[k])
          leaf_data[obs_nodes[k]] = mcmc.event_observations.seq[k]
        end
      end
      lposterior += PhyloModels.loglikelihood(tree, sm, leaf_data)
    end
  end
  if lposterior > max_lposterior
    if M <: PhyloILM
      return MarkovChain{S, M}(events, network, rparams, sm, lposterior)
    elseif M <: TNILM
      return MarkovChain{S, M}(events, network, rparams, lposterior)
    end
  else
    return nothing
  end
end

function initialize(::Type{MarkovChain},
                    mcmc::MCMC{S, M},
                    progress_channel::RemoteChannel;
                    attempts::Int64=1000) where {
                    S <: DiseaseStateSequence,
                    M <: ILM}
  if attempts <= 0
    throw(ErrorException("Must have at least 1 initialization attempt"))
  end
  local markov_chain
  max_lposterior = -Inf
  for i in 1:attempts
    initialization = _initialization_attempt(mcmc, max_lposterior)
    if initialization !== nothing
      markov_chain = initialization
      max_lposterior = markov_chain.log_posterior[1]
    end
    put!(progress_channel, true)
  end
  if max_lposterior == -Inf
    @error "Failed to initialize Markov Chain"
    markov_chain = nothing
  end
  return markov_chain
end


function initialize(::Type{MarkovChain},
                    mcmc::MCMC{S, M};
                    attempts::Int64=1000) where {
                    S <: DiseaseStateSequence,
                    M <: ILM}
  if attempts <= 0
    throw(ErrorException("Must have at least 1 initialization attempt"))
  end
  max_lposterior = -Inf
  local markov_chain
  pmeter = Progress(attempts, "Initialization progress")
  for i in 1:attempts
    initialization = _initialization_attempt(mcmc, max_lposterior)
    if initialization !== nothing
      markov_chain = initialization
      max_lposterior = markov_chain.log_posterior[1]
    end
    next!(pmeter)
  end
  finish!(pmeter)
  if max_lposterior == -Inf
    @error "Failed to initialize Markov Chain"
    markov_chain = nothing
  end
  return markov_chain
end