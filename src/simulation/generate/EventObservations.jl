function generate(
  ::Type{EventObservations},
  sim::Simulation{S, M},
  delay_infection::UnivariateDistribution,
  delay_removal::UnivariateDistribution;
  seq_len::Union{Nothing, Int64}=nothing,
  force::Bool = false) where {
    S <: Union{SEIR, SIR},
    M <: ILM}
  infection = fill(NaN, individuals(sim.events))
  removal = fill(NaN, individuals(sim.events))
  if force
    for i in findall(sim.events.infection .!== NaN)
      if sim.events.infection[i] == -Inf
        infection[i] = -Inf
        if sim.events.removal[i] == -Inf
          removal[i] = -Inf
        else
          removal[i] = sim.events.removal[i] + rand(delay_removal)
        end
      else
        if isnan(sim.events.removal[i])
          infection[i] = sim.events.infection[i] + rand(delay_infection)
        elseif sim.events.removal[i] > -Inf
          infection_delay_ub = sim.events.removal[i] - sim.events.infection[i]
          infection[i] = sim.events.infection[i] +
            rand(Truncated(delay_infection, 0.0, infection_delay_ub))
          removal[i] = sim.events.removal[i] + rand(delay_removal)
        end
      end
      @debug "Infection observation of i = $i at t = $(round(infection[i], digits=3)) (actual infection onset at t = $(round(sim.events.infection[i], digits=3)))"
      @debug "Removal observation of i = $i at t = $(round(removal[i], digits=3)) (actual removal at t = $(round(sim.events.removal[i], digits=3)))"
    end
  else
    for i in findall(sim.events.infection .!== NaN)
      if sim.events.infection[i] == -Inf
        infection[i] = -Inf
        if sim.events.removal == -Inf
          removal[i] = -Inf
        else
          removal[i] = sim.events.removal[i] + rand(delay_removal)
        end
      else
        infection_delay = rand(delay_infection)
        if sim.events.removal[i] === NaN
          infection[i] = sim.events.infection[i] + infection_delay
        else
          if infection_delay + sim.events.infection[i] < sim.events.removal[i]
            infection[i] = sim.events.infection[i] + infection_delay
            removal[i] = sim.events.removal[i] + rand(delay_removal)
          end
        end
      end
      @debug "Infection observation of i = $i t = $(round(infection[i], digits=3)) (actual infection onset at t = $(round(sim.events.infection[i], digits=3)))"
      @debug "Removal observation of i = $i at t = $(round(removal[i], digits=3)) (actual removal at t = $(round(sim.events.removal[i], digits=3)))"
    end
  end
  if M == PhyloILM
    seq_len == nothing && error("Must specify a sequence length with $(M)s")
    tree, event_nodes = generate(Tree,
                                sim.events,
                                infection,
                                sim.transmission_network)
    seq_full = simulate(RNASeq, tree, sim.substitution_model, seq_len)
    seq_obs = Union{Nothing, RNASeq}[isnothing(event_nodes[i]) ? nothing : seq_full[event_nodes[i]] for i = eachindex(infection)]
    return EventObservations{S}(infection, removal, seq_obs)
  elseif M == TNILM
    return EventObservations{S}(infection, removal)
  end
end

function generate(
  ::Type{EventObservations},
  sim::Simulation{S, M},
  delay_infection::UnivariateDistribution;
  seq_len::Union{Nothing, Int64}=nothing) where {
    S <: Union{SEI, SI},
    M <: ILM}
  infection = fill(NaN, individuals(sim.events))
  @simd for i in findall(.!isnan.(sim.events.infection))
    if sim.events.infection[i] == -Inf
      infection[i] = -Inf
    else
      infection[i] = sim.events.infection[i] + rand(delay_infection)
    end
    @debug "Infection observation of i = $i at t = $(round(infection[i], digits=3))) (actual infection onset at t = $(round(sim.events.infection[i], digits=3)))"
  end
  if M == PhyloILM
    seq_len == nothing && error("Must specify a sequence length with $(M)s")
    tree, event_nodes = generate(Tree,
                                sim.events,
                                infection,
                                sim.transmission_network)
    seq_full = simulate(RNASeq, tree, sim.substitution_model, seq_len)
    seq_obs = Union{Nothing, RNASeq}[isnothing(event_nodes[i]) ? nothing : seq_full[event_nodes[i]] for i = eachindex(infection)]
    return EventObservations{S}(infection, seq_obs)
  elseif M == TNILM
    return EventObservations{S}(infection)
  end
end

observe(x,y,z; force::Bool=false, seq_len::Union{Nothing, Int64}=nothing) = generate(EventObservations, x, y, z, force=force, seq_len=seq_len)
observe(x,y; seq_len::Union{Nothing, Int64}=nothing) = generate(EventObservations, x, y, seq_len=seq_len)