function _pathway_from(
  id::Int64,
  network::TransmissionNetwork;
  depth::Real=Inf,
  include_id::Bool=true)

  local path
  if network.external[id] || any(network.internal[:, id])
    path = Int64[]
    push!(path, id)
    path_lengths = [0]
    while (length(path) > path_lengths[end]) & (depth >= length(path_lengths))
      push!(path_lengths, length(path))
      for j in path[(path_lengths[end-1]+1):path_lengths[end]]
        append!(path, findall(network.internal[j, :]))
      end
    end
    if include_id == false
      popfirst!(path)
    end
  else
    path = nothing
  end
  @logmsg LogLevel(-5000) "Transmission pathway from Individual $id: $path"
  return path
end


function _pathway_from(
  id::Nothing,
  network::TransmissionNetwork;
  depth::Real=Inf,
  include_id::Bool=true)

  path = Union{Nothing, Int64}[nothing]
  append!(path, findall(network.external))
  path_lengths = [0; 1]
  while (length(path) > path_lengths[end]) & (depth >= length(path_lengths)-1)
    push!(path_lengths, length(path))
    for j in path[(path_lengths[end-1]+1):path_lengths[end]]
      append!(path, findall(network.internal[j, :]))
    end
  end
  if include_id == false
    popfirst!(path)
  end
  @logmsg LogLevel(-5000) "Transmission pathway from external source: $path"
  return path
end

function _pathway_to(
  id::Int64,
  network::TransmissionNetwork;
  depth::Real=Inf,
  include_id::Bool=true)

  local path
  if network.external[id] || any(network.internal[:, id])
    path = Union{Nothing, Int64}[id]
    next_id = findfirst(network.internal[:, id])
    while (length(path) < depth) && (next_id !== nothing)
      push!(path, next_id)
      next_id = findfirst(network.internal[:, path[end]])
    end
    push!(path, next_id)
    if include_id == false
      popfirst!(path)
    end
  else
    path = nothing
  end
  @logmsg LogLevel(-5000) "Transmission pathway to Individual $id: $path"
  return path
end

function _source(id::Int64, network::TransmissionNetwork)
  source = findfirst(network.internal[:, id])
  @logmsg LogLevel(-5000) "Transmission source lookup" Individual = id Source = source
  return source
end

function generate(
  ::Type{Tree},
  events::Events{S},
  obs::Vector{Float64},
  network::TransmissionNetwork) where {
  S <: DiseaseStateSequence}

  if length(obs) != individuals(events)
    throw(ErrorException("Infection observation vector does not match number of individuals in population"))
  end

  # Initialization
  expos = 0
  tree = Tree()
  addnode!(tree) # Root/MRCA node
  local event_times
  local parent_node
  local branch_length

  if State_E in S
    event_times = [events.exposure obs]
  else
    event_times = [events.infection obs]
  end
  event_nodes = Array{Union{Nothing, Int64}, 2}(fill(nothing, size(event_times)))
  event_order = sortperm(event_times[:])
  event_lookup = CartesianIndices(size(event_times))
  @logmsg LogLevel(-5000) "event_nodes, event_order, and event_lookup created" ArraySize = size(event_times)

  pathways = [_pathway_from(i, network, include_id = true) for i = 1:individuals(events)]

  # Determine significant transmissions (results in a observation in the future)
  significant_transmissions = (event_times[:, 1] .== -Inf) .| [pathways[i] !== nothing && any(event_times[pathways[i], 2] .>= event_times[i, 1]) for i = 1:individuals(events)]
  @logmsg LogLevel(-5000) "Significant transmitters" significant_transmissions

  start_time = minimum(events)

  # Iterate through all events to build tree
  for i = eachindex(event_order)

    event_time = event_times[event_order[i]]

    # Stop when no remaining valid event times
    isnan(event_time) && break

    # Determine the individual and event type
    event_individual, event_type = Tuple(event_lookup[event_order[i]])

    @debug "Building phylogenetic tree..." Event = i Individual = event_individual Type = event_type == 1 ? "Transmission" : "Observation" Time = event_time Index = event_order[i]

    # For purposes of tree generation, replace -Inf times (initial condition coding) with start_time
    if event_time == -Inf
      @logmsg LogLevel(-5000) "Setting -Inf event time to $(start_time)"
      event_time = start_time
    end

    # For transmission events...
    if event_type == 1

      # Add transmission event to tree only if it is significant
      if significant_transmissions[event_individual]

        # For significant external transmission events...
        if network.external[event_individual]
          @logmsg LogLevel(-5000) "Event $i is a significant external transmission"
          parent_node = 1
          expos += 1
          if expos > 1
            @warn "Current support is intended for a single intial event with no external exposures" ExternalExposures = expos
          end
          branch_length = event_time - start_time
          @logmsg LogLevel(-5000) "Set parent node to the root" ParentNode = parent_node BranchLength = branch_length
          branch!(tree, parent_node, branch_length)

          # Record tree and node ID
          event_nodes[event_order[i]] = length(tree.nodes)

        # For other significant transmission events...
        else

          # Determine source of transmission
          source = _source(event_individual, network)

          # If initial condition, no source
          if source === nothing
            @warn "Shouldn't be!"

          # Previous significant transmissions from this source of transmission
          else
            @logmsg LogLevel(-5000) "Event $i is a significant internal transmission" Source = source
            source_prior_transmissions = findall(network.internal[source, :][:] .& (event_times[:, 1] .< event_time) .& significant_transmissions)

            # For when the source of transmission has yet to be detected...
            if isnan(event_times[source, 2]) || event_times[source, 2] > event_times[event_individual, 1]
              @logmsg LogLevel(-5000) "Individual $source, the transmission source for individual $event_individual had not yet been detected"
              # When there have been prior transmissions from this undetected source
              if length(source_prior_transmissions) > 0
                @logmsg LogLevel(-5000) "Individual $source has had earlier transmissions" PriorTransmissionIds = source_prior_transmissions PriorTransmissionTimes = event_times[source_prior_transmissions, 1] PriorTransmissioNodes = event_nodes[source_prior_transmissions, 1]
                parent_node = event_nodes[source_prior_transmissions[argmax(event_times[source_prior_transmissions, 1])], 1]
                branch_length = event_time - maximum(event_times[source_prior_transmissions, 1])

              # When there have not been any prior transmissions from this undetected source
              else
                parent_node = event_nodes[source, 1]
                branch_length = event_time - max(event_times[source, 1], start_time)
                @logmsg LogLevel(-5000) "Individual $source has had no earlier transmissions" Source = source ParentNode1 = event_nodes[source, 1] ParentNode2 = parent_node BranchLength = branch_length EventTime = event_time ParentNodeTime = event_times[source, 1]
              end

            # For when the source of transmission has been detected...
            else
              @logmsg LogLevel(-5000) "Individual $source, the transmission source for individual $event_individual had been already detected"
              # And detection of source of transmission is most recent, relevant event...
              if length(source_prior_transmissions) == 0 || all(event_times[source, 2] .> event_times[source_prior_transmissions, 1])
                @logmsg LogLevel(-5000) "The detection of individual $source is the recentmost relevant event to transmission of Individual $event_individual"
                # Determine parent node of the observation leaf node
                parent_node = parentnode(tree, event_nodes[source, 2])
                branch_length = event_time - max(event_times[source, 2], start_time)

              # Or some prior transmission is most recent, relevant event...
              else
                @logmsg LogLevel(-5000) "Individual $source has had earlier transmissions, after detection" PriorTransmissionIds = source_prior_transmissions PriorTransmissionTimes = event_times[source_prior_transmissions, 1] PriorTransmissioNodes = event_nodes[source_prior_transmissions, 1]
                parent_node = event_nodes[source_prior_transmissions[argmax(event_times[source_prior_transmissions, 1])], 1]
                branch_length = event_time - maximum(event_times[source_prior_transmissions, 1])
              end
            end
            # Add branch and node
            @logmsg LogLevel(-5000) "branch! information" Individual = event_individual Source = source ParentNode = parent_node Node = length(tree.nodes) + 1 BranchLength = branch_length
            branch!(tree, parent_node, branch_length)

            # Record node id
            event_nodes[event_order[i]] = length(tree.nodes)
          end
        end
      else
        @logmsg LogLevel(-5000) "Event $i is an insignificant transmission"
      end

    # Infection observation event
    elseif event_type == 2
      # Determine if individual has had significant transmissions prior to being observed as infected
      prior_transmissions = findall(network.internal[event_individual, :][:] .&
                            (event_times[:, 1] .< event_time) .&
                            significant_transmissions)

      if length(prior_transmissions) > 0
        # Parent node will be the internal node associated with the final pre-observation transmission
        @logmsg LogLevel(-5000) "Individual $event_individual transmitted to others prior to their detection" PriorTransmissionIds = prior_transmissions PriorTransmissionTimes = event_times[prior_transmissions, 1] PriorTransmissioNodes = event_nodes[prior_transmissions, 1]
        parent_id = prior_transmissions[argmax(event_times[prior_transmissions, 1])]
        parent_node = event_nodes[parent_id, 1]
        branch_length = event_time - event_times[parent_id, 1]
      else
        # Individual has not transmitted to others prior to being observed as infected
        @logmsg LogLevel(-5000) "Individual $event_individual has not transmitted to others prior to their detection"
        parent_node = event_nodes[event_individual, 1]
        branch_length = event_time - max(event_times[event_individual, 1], start_time)
      end

      if any(network.internal[event_individual, :][:] .& (event_times[:, 1] .>= event_time) .& significant_transmissions)
        # Individual goes on to have significant transmissions
        # Internal node for observation
        @logmsg LogLevel(-5000) "Individual $event_individual has significant future transmissions, add internal node for observation" ParentNode = parent_node BranchLength = branch_length
        branch!(tree, parent_node, branch_length)
        # Add zero-length branch so observation event will be a leaf node
        @logmsg LogLevel(-5000) "Add a zero length branch to ensure observation of individual $event_individual is a leaf node" ParentNode = length(tree.nodes) BranchLength = 0.0
        branch!(tree, length(tree.nodes), 0.0)

      else
        # Individual does not go one to have any significant transmissions
        # so will surely be a leaf node
        @logmsg LogLevel(-5000) "Individual $event_individual does not have significant future transmissions - their infection observation represents a leaf node" ParentNode = parent_node BranchLength = branch_length
        branch!(tree, parent_node, branch_length)
      end
      # Record leaf node id
      event_nodes[event_order[i]] = length(tree.nodes)
    end
  end
  return tree, event_nodes[:,2]
end

function generate(
  ::Type{Tree},
  events::Events{S},
  obs::EventObservations{S, M},
  network::TransmissionNetwork) where {
  S <: DiseaseStateSequence,
  M <:ILM}
  return generate(Tree, events, obs.infection, network)
end