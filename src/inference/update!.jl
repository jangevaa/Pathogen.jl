function update!(mc::MarkovChain{T, M},
                 mcmc::MCMC{T, M},
                 Σ::Array{Float64, 2},
                 σ::Float64;
                 condition_on_network::Bool=false,
                 event_batches::Int64=1,
                 transmission_rate_cache::Union{Nothing, TransmissionRateCache}=nothing) where {
                 T <: DiseaseStateSequence,
                 M <: TNILM}
  # Initialize
  new_tr_cache = transmission_rate_cache
  new_events = mc.events[end]
  new_events_array = new_events[convert(DiseaseStates, T)[2:end]]
  new_rp = mc.risk_parameters[end]
  new_network = mc.transmission_network[end]
  new_lposterior = mc.log_posterior[end]
  # Randomize event time augmentation
  event_indices = findall(new_events_array[:] .> -Inf)
  aug_order = sample(event_indices, length(event_indices), replace=false)
  @debug "$(length(aug_order)) events to augment"
  if event_batches < 0
    @error "Cannot have negative amount of event batches"
  end
  if event_batches > length(aug_order)
    @warn "More event batches than there are events to augment ($(event_batches) > $(length(aug_order))), setting to maximum ($(length(aug_order)))"
    event_batches = length(aug_order)
  end
  batch_size = fld(length(aug_order), event_batches)
  @debug "Performing data augmentation in batches of $batch_size events at a time"
  for i = 1:(1 + event_batches)
    if i == 1
      # Propose new risk parameters
      proposed_events = new_events
      proposed_events_array = new_events_array
      proposed_rp = generate(RiskParameters{T}, new_rp, Σ)
      @debug "Generating new TransmissionRateCache for TN-ILM parameter proposal"
      proposed_tr_cache = TransmissionRateCache(size(new_events_array, 1))
    else
      proposed_events_array = copy(new_events_array)
      proposed_events = Events{T}(proposed_events_array)
      for j = (batch_size*(i-2) + 1):min(batch_size*(i-1),length(aug_order))
        id, state_index = Tuple(CartesianIndices((individuals(new_events),
                                            length(convert(DiseaseStates, T)[2:end])))[aug_order[j]])
        new_state = convert(DiseaseStates, T)[state_index+1]
        time = new_events[new_state][id]

        # Conditioning on network means that only event times valid under current network will be proposed.
        # This is useful for models which may have additional contributions to the posterior based on network,
        # and require more modest proposals (e.g. phylodynamic models).
        if condition_on_network
          proposed_event = generate(Event,
                                    Event{T}(time, id, new_state),
                                    σ,
                                    mcmc.event_extents,
                                    mcmc.event_observations,
                                    proposed_events,
                                    new_network)
        else
          proposed_event = generate(Event,
                                    Event{T}(time, id, new_state),
                                    σ,
                                    mcmc.event_extents,
                                    mcmc.event_observations,
                                    proposed_events)
        end
        proposed_events_array = reshape([proposed_events_array[1:(aug_order[j]-1)]
                                         _time(proposed_event)
                                         proposed_events_array[(aug_order[j]+1):end]],
                                        size(new_events_array))
        proposed_events = Events{T}(proposed_events_array)
      end
      proposed_rp = new_rp
      proposed_tr_cache = new_tr_cache
    end
    proposed_lprior = logprior(proposed_rp, mcmc.risk_priors)
    # Based on the logprior and competiting MCMC iteration, this loglikelihood is required for acceptance
    # Calculating this in advance allows us to cut loglikelihood calculation short if it goes below threshold
    ll_acceptance_threshold = log(rand()) + new_lposterior - proposed_lprior
    if ll_acceptance_threshold < Inf
      proposed_llikelihood = loglikelihood(proposed_rp,
                                           mcmc.risk_functions,
                                           proposed_events,
                                           mcmc.population,
                                           mcmc.starting_states,
                                           transmission_rates_output = false,
                                           transmissions_output = false,
                                           early_decision_value = ll_acceptance_threshold,
                                           transmission_rate_cache = proposed_tr_cache)[1]
      proposed_lposterior = proposed_lprior + proposed_llikelihood
    else
      proposed_llikelihood = -Inf
      proposed_lposterior = -Inf
    end
    if proposed_llikelihood >= ll_acceptance_threshold
      @debug "MCMC proposal accepted (acceptance probability = $(round(min(1.0, exp(proposed_lposterior - new_lposterior)), digits=3)))"
      new_rp = proposed_rp
      new_tr_cache = proposed_tr_cache
      new_events = proposed_events
      new_events_array = proposed_events_array
      new_lposterior = proposed_lposterior
    else
      @debug "next!: MCMC proposal rejected (acceptance probability = $(round(min(1.0, exp(proposed_lposterior - new_lposterior)), digits=3)))"
    end
  end
  mc.iterations += 1
  tnr, tx = loglikelihood(new_rp,
                          mcmc.risk_functions,
                          new_events,
                          mcmc.population,
                          mcmc.starting_states,
                          loglikelihood_output = false,
                          transmission_rate_cache = new_tr_cache)[[2; 3]]
  new_network = generate(TransmissionNetwork, tnr, mcmc, tx)
  push!(mc.events, new_events)
  push!(mc.transmission_network, new_network)
  push!(mc.risk_parameters, new_rp)
  push!(mc.log_posterior, new_lposterior)
  return mc
end

function update!(mc::MarkovChain{T, M},
                 mcmc::MCMC{T, M},
                 Σrp::Array{Float64, 2},
                 Σsm::Array{Float64, 2},
                 σ::Float64;
                 event_batches::Int64=1,
                 transmission_rate_cache::Union{Nothing, TransmissionRateCache}=nothing,
                 DA_couple::Bool=false) where {
                 T <: DiseaseStateSequence,
                 M <: PhyloILM}
  # Initialize
  local tnr
  new_tr_cache = transmission_rate_cache
  new_events = mc.events[end]
  new_events_array = convert(Array{Float64, 2}, new_events)
  new_rp = mc.risk_parameters[end]
  new_sm = mc.substitution_model[end]
  new_network = mc.transmission_network[end]
  new_lposterior = mc.log_posterior[end]
  # Randomize event time augmentation
  event_indices = findall(new_events_array[:] .> -Inf)
  aug_order = sample(event_indices, length(event_indices), replace=false)
  @debug "$(length(aug_order)) events to augment"
  if event_batches < 0
    @error "Cannot have negative amount of event batches"
  end
  if event_batches > length(aug_order)
    @warn "More event batches than there are events to augment ($(event_batches) > $(length(aug_order))), setting to maximum ($(length(aug_order)))"
    event_batches = length(aug_order)
  end
  batch_size = fld(length(aug_order), event_batches)
  @debug "Performing data augmentation in batches of $batch_size events at a time"
  steps = DA_couple ? 1 + event_batches : 1 + event_batches + 1
  for i = 1:steps
    if i == 1
      # Propose new risk parameters
      proposed_events = new_events
      proposed_events_array = new_events_array
      proposed_rp = generate(RiskParameters{T}, new_rp, Σrp)
      proposed_sm = generate(mcmc.substitution_model, new_sm, Σsm)
      @debug "Generating new TransmissionRateCache for TN-ILM parameter proposal"
      proposed_tr_cache = TransmissionRateCache(size(new_events_array, 1))
    elseif i < (event_batches + 2)
      proposed_events_array = copy(new_events_array)
      proposed_events = Events{T}(proposed_events_array)
      t_ids = Int64[]
      for j = (batch_size*(i-2) + 1):min(batch_size*(i-1),length(aug_order))
        id, state_index = Tuple(CartesianIndices((individuals(new_events),
                                            length(convert(DiseaseStates, T)[2:end])))[aug_order[j]])
        # Generate new transmission source for transmission events
        state_index == 1 && push!(t_ids, id)
        new_state = convert(DiseaseStates, T)[state_index+1]
        time = new_events[new_state][id]
        proposed_event = generate(Event,
                                  Event{T}(time, id, new_state),
                                  σ,
                                  mcmc.event_extents,
                                  mcmc.event_observations,
                                  proposed_events,
                                  new_network)
        proposed_events_array = reshape([proposed_events_array[1:(aug_order[j]-1)]
                                         _time(proposed_event)
                                         proposed_events_array[(aug_order[j]+1):end]],
                                        size(new_events_array))
        proposed_events = Events{T}(proposed_events_array)
      end
      proposed_rp = new_rp
      proposed_sm = new_sm
      proposed_tr_cache = new_tr_cache
    elseif i == (event_batches + 2)
      proposed_events = new_events
      proposed_events_array = new_events_array
      proposed_rp = new_rp
      proposed_sm = new_sm
      proposed_tr_cache = new_tr_cache
    end
    proposed_lprior = logprior(proposed_rp, mcmc.risk_priors)
    proposed_lprior += logprior(proposed_sm, mcmc.substitution_model_prior)
    # Based on the logprior and competiting MCMC iteration, this loglikelihood is required for acceptance
    # Calculating this in advance allows us to cut loglikelihood calculation short if it goes below threshold
    ll_acceptance_threshold = log(rand()) + new_lposterior - proposed_lprior
    if ll_acceptance_threshold < Inf
      proposed_llikelihood, tnr, tx = loglikelihood(proposed_rp,
                                                    mcmc.risk_functions,
                                                    proposed_events,
                                                    mcmc.population,
                                                    mcmc.starting_states,
                                                    early_decision_value = ll_acceptance_threshold,
                                                    transmissions_output = true,
                                                    transmission_rate_cache = proposed_tr_cache)
      if proposed_llikelihood > ll_acceptance_threshold
        if DA_couple
          if i == 1 || length(t_ids) == 0
            proposed_network = new_network
          elseif 1 < i < event_batches + 2
            proposed_network = generate(TransmissionNetwork, new_network, tnr, mcmc, t_ids)
            @debug "Generating a Transmission Network proposal" ProposedNetwork = proposed_network
          end
        else
          if i < event_batches + 2
            proposed_network = new_network
          elseif i == event_batches + 2
            proposed_network = generate(TransmissionNetwork, new_network, tnr, mcmc, sample(tx))
            @debug "Generating a transmission network proposal" ProposedNetwork = proposed_network
          end
        end
        @debug "Generating a Tree for proposed network and event times"
        tree, obs_nodes = generate(Tree, proposed_events, mcmc.event_observations, proposed_network)
        leaf_data = Dict{Int64, RNASeq}()
        for k = eachindex(obs_nodes)
          if !isnothing(obs_nodes[k])
            leaf_data[obs_nodes[k]] = mcmc.event_observations.seq[k]
          end
        end
        proposed_llikelihood += PhyloModels.loglikelihood(tree, proposed_sm, leaf_data)
        proposed_lposterior = proposed_lprior + proposed_llikelihood
      end
    else
      proposed_llikelihood = -Inf
      proposed_lposterior = -Inf
    end
    if proposed_llikelihood >= ll_acceptance_threshold
      @debug "MCMC proposal accepted (acceptance probability = $(round(min(1.0, exp(proposed_lposterior - new_lposterior)), digits=3)))"
      new_rp = proposed_rp
      new_sm = proposed_sm
      new_tr_cache = proposed_tr_cache
      new_events = proposed_events
      new_events_array = proposed_events_array
      new_network = proposed_network
      new_lposterior = proposed_lposterior
    else
      @debug "MCMC proposal rejected (acceptance probability = $(round(min(1.0, exp(proposed_lposterior - new_lposterior)), digits=3)))"
    end
  end
  mc.iterations += 1
  push!(mc.events, new_events)
  push!(mc.transmission_network, new_network)
  push!(mc.risk_parameters, new_rp)
  push!(mc.substitution_model, new_sm)
  push!(mc.log_posterior, new_lposterior)
  return mc
end