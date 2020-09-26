function _bounds(
  id::Int64,
  new_state::DiseaseState,
  extents::EventExtents{T},
  obs::EventObservations{T, M},
  events::Events{T},
  network::TransmissionNetwork) where {
    T <: DiseaseStateSequence,
    M <: ILM}

  if new_state == State_E
    lowerbounds = [events.infection[id] - extents.exposure[2]]
    upperbounds = [events.infection[id] - extents.exposure[1]]
    path_to = _pathway_to(id, network, depth = 1)
    if path_to != nothing
      if path_to[2] != nothing
        parent_host = path_to[2]
        push!(lowerbounds, events.infection[parent_host])
        if (T == SEIR) && !isnan(events.removal[parent_host])
          push!(upperbounds, events.removal[parent_host])
        end
      end
    else
      error("Event inconsistent with TransmissionNetwork")
    end
  elseif new_state == State_I
    path_from = _pathway_from(id, network, depth = 1)
    if State_E ∈ T
      lowerbounds = [obs.infection[id] - extents.infection[2]
                     events.exposure[id] + extents.exposure[1]]
      upperbounds = [obs.infection[id] - extents.infection[1]
                     events.exposure[id] + extents.exposure[2]]
      if length(path_from) > 1
        child_hosts = path_from[2:end]
        append!(upperbounds, events.exposure[child_hosts])
      end
    elseif State_E ∉ T
      lowerbounds = [obs.infection[id] - extents.infection[2]]
      upperbounds = [obs.infection[id] - extents.infection[1]]
      path_to = _pathway_to(id, network, depth = 1)
      if length(path_from) > 1
        child_hosts = path_from[2:end]
        append!(upperbounds, events.infection[child_hosts])
      end
      if path_to != nothing
        if path_to[2] != nothing
          parent_host = path_to[2]
          push!(lowerbounds, events.infection[parent_host])
          if (T == SIR) && !isnan(events.removal[parent_host])
            push!(upperbounds, events.removal[parent_host])
          end
        end
      else
        error("Event inconsistent with TransmissionNetwork")
      end
    end
  elseif new_state == State_R
    path_from = _pathway_from(id, network, depth = 1)
    lowerbounds = [obs.removal[id] - extents.removal[2]
                   obs.infection[id]]
    upperbounds = [obs.removal[id] - extents.removal[1]]
    if length(path_from) > 1
      child_hosts = path_from[2:end]
      if T == SEIR
        append!(lowerbounds, events.exposure[child_hosts])
      elseif T == SIR
        append!(lowerbounds, events.infection[child_hosts])
      end
    end
  end
  @debug "Transition of $id into $new_state with bounds: \n  [max($(round.(lowerbounds, digits=3))),\n   min($(round.(upperbounds, digits=3)))]"
  lowerbound = maximum(lowerbounds)
  upperbound = minimum(upperbounds)
  if lowerbound >= upperbound
    @warn "Lower >= upper bound for transition of i=$id into $new_state" lower = lowerbounds upper = upperbounds
  end
  return lowerbound, upperbound
end


function _bounds(
  last_event::Event{T},
  extents::EventExtents{T},
  obs::EventObservations{T, M},
  events::Events{T},
  network::TransmissionNetwork) where {
    T <: DiseaseStateSequence,
    M <: ILM}

  return _bounds(last_event.individual,
                 last_event.new_state,
                 extents, obs, events, network)
end


function _bounds(
  id::Int64,
  new_state::DiseaseState,
  extents::EventExtents{T},
  obs::EventObservations{T, M},
  events::Events{T}) where {
    T <: DiseaseStateSequence,
    M<: ILM}

  if new_state == State_E
    lowerbounds = [events.infection[id] - extents.exposure[2]]
    upperbounds = [events.infection[id] - extents.exposure[1]]
  elseif new_state == State_I
    if State_E ∈ T
      lowerbounds = [obs.infection[id] - extents.infection[2]
                     events.exposure[id] + extents.exposure[1]]
      upperbounds = [obs.infection[id] - extents.infection[1]
                     events.exposure[id] + extents.exposure[2]]
    elseif State_E ∉ T
      lowerbounds = [obs.infection[id] - extents.infection[2]]
      upperbounds = [obs.infection[id] - extents.infection[1]]
    end
  elseif new_state == State_R
    lowerbounds = [obs.removal[id] - extents.removal[2]
                   obs.infection[id]]
    upperbounds = [obs.removal[id] - extents.removal[1]]
  end
  @debug "Transition of $id into $new_state with bounds: \n  [max($(round.(lowerbounds, digits=3))),\n   min($(round.(upperbounds, digits=3)))]"
  lowerbound = maximum(lowerbounds)
  upperbound = minimum(upperbounds)
  if lowerbound >= upperbound
    @warn "Lower >= upper bound for transition of i=$id into $new_state" lower = lowerbounds upper = upperbounds
  end
  return lowerbound, upperbound
end


function _bounds(
  last_event::Event{T},
  extents::EventExtents{T},
  obs::EventObservations{T, M},
  events::Events{T}) where {
    T <: DiseaseStateSequence,
    M <: ILM}

  return _bounds(last_event.individual,
                 last_event.new_state,
                 extents, obs, events)
end
