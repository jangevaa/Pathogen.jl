"""
infer.jl
"""

function surveil(population::Population, ν::Float64)
  """
  Gather surveillance data on specific individuals in a population, with an exponentially distributed detection lag with rate ν
  """
  exposed_actual = fill(NaN, length(population.events)-1)
  infectious_actual = fill(NaN, length(population.events)-1)
  infectious_observed = fill(NaN, length(population.events)-1)
  removed_actual = fill(NaN, length(population.events)-1)
  removed_observed = fill(NaN, length(population.events)-1)
  covariates_actual = fill(fill(NaN, length(population.history[2][1][1])),  length(population.events)-1)
  covariates_observed = fill(fill(NaN, length(population.history[2][1][1])),  length(population.events)-1)
  seq_actual = convert(Vector{Any}, fill(NaN, length(population.events)-1))
  seq_observed = convert(Vector{Any}, fill(NaN, length(population.events)-1))

  for i = 2:length(population.events)
    # Initial conditions (assumed to be constant and observed without error)
    covariates_actual[i-1] = population.history[i][1][1]
    covariates_observed[i-1] = population.history[i][1][1]

    # Exposure time (unobservable)
    if length(population.events[i][1]) > 0
      exposed_actual[i-1] = population.events[i][1][1]
    end

    # Infectious time (observed with latency)
    if length(population.events[i][3]) > 0
      infectious_actual[i-1] = population.events[i][3][1]
      seq_actual[i-1] = convert(Vector{Int64}, population.history[i][2][find(infectious_actual[i-1] .>= population.events[i][6])[end]])
      if ν < Inf
        infectious_observed[i-1] = infectious_actual[i-1] + rand(Exponential(1/ν))
      elseif ν == Inf
        infectious_observed[i-1] = infectious_actual[i-1]
      end

      if length(population.events[i][4]) > 0 && infectious_observed[i-1] >= population.events[i][4][1]
        infectious_observed[i-1] = NaN
      else
        seq_observed[i-1] = convert(Vector{Int64}, population.history[i][2][find(infectious_observed[i-1] .>= population.events[i][6])[end]])
      end
    end

    # Removal time (observed with latency)
    if length(population.events[i][4]) > 0
      removed_actual[i-1] = population.events[i][4][1]
      if !isnan(infectious_observed[i-1])
        if ν < Inf
          removed_observed[i-1] = removed_actual[i-1] + rand(Exponential(1/ν))
        elseif ν == Inf
          removed_observed[i-1] = removed_actual[i-1]
        end
      end
    end
  end
  return SEIR_actual(exposed_actual, infectious_actual, removed_actual, covariates_actual, seq_actual), SEIR_observed(infectious_observed, removed_observed, covariates_observed, seq_observed)
end


surveil(population, Inf) = surveil(population::Population)


function augment(ρ::Float64, ν::Float64, network::Array{Bool, 2}, obs::SEIR_observed, debug=false::Bool)
  """
  Augments surveilance data, organizes observations, based on ρ, ν, and a transmission network
  """
  exposed_augmented = fill(NaN, length(obs.infectious))
  infectious_augmented = fill(NaN, length(obs.infectious))
  removed_augmented = fill(NaN, length(obs.removed))

  # Determine which event times have additional restrictions...
  unrestricted = find(network[1,:])
  restricted = find(network[unrestricted[1], :])
  restricted_lengths = [0]
  for i = unrestricted[2:end]
    append!(restricted, find(network[i, :]))
  end
  while length(restricted)-restricted_lengths[end] > 0
    push!(restricted_lengths, length(restricted))
    for i = restricted[(restricted_lengths[end-1]+1):end]
      append!(restricted, find(network[i, :]))
    end
  end

  for i = unrestricted
    if ν < Inf
      infectious_augmented[i] = obs.infectious[i] - rand(Exponential(1/ν))
    elseif ν == Inf
      infectious_augmented[i] = obs.infectious[i]
    end
    exposed_augmented[i] = infectious_augmented[i] - rand(Exponential(1/ρ))
    if ν < Inf
      removed_augmented[i] = obs.removed[i] - rand(Truncated(Exponential(1/ν), -Inf, obs.removed[i] - obs.infectious[i]))
    elseif ν == Inf
      removed_augmented[i] = obs.removed[i]
    end
  end

  for i = restricted
    source = findfirst(network[:,i])
    if ν < Inf
      infectious_augmented[i] = obs.infectious[i] - rand(Truncated(Exponential(1/ν), -Inf, obs.infectious[i] - infectious_augmented[source]))
    elseif ν == Inf
      infectious_augmented[i] = obs.infectious[i]
    end
    if isnan(removed_augmented[source])
      exposed_augmented[i] = infectious_augmented[i] - rand(Truncated(Exponential(1/ρ), -Inf, infectious_augmented[i]-infectious_augmented[source]))
    else
      exposed_augmented[i] = infectious_augmented[i] - rand(Truncated(Exponential(1/ρ), infectious_augmented[i]-removed_augmented[source], infectious_augmented[i]-infectious_augmented[source]))
    end
    if ν < Inf
      removed_augmented[i] = obs.removed[i] - rand(Truncated(Exponential(1/ν), -Inf, obs.removed[i] - obs.infectious[i]))
    elseif ν == Inf
      removed_augmented[i] = obs.removed[i]
    end
  end
  if debug
    println("DATA AUGMENTATION")
    println("Infected individuals: $(find(sum(network,1)))")
    println("Unrestricted individuals: $unrestricted")
    println("Restricted individuals: $restricted")
    println("Augmented exposure times ($(sum(!isnan(exposed_augmented))) total): $(round(exposed_augmented,3))")
    println("Augmented infection times ($(sum(!isnan(infectious_augmented))) total): $(round(infectious_augmented,3))")
    println("Augmented removal times ($(sum(!isnan(removed_augmented))) total): $(round(removed_augmented,3))")
    println("")
  end
  if debug
    for i = 1:length(obs.infectious)
      @assert(!(isnan(exposed_augmented[i]) && !isnan(obs.infectious[i])), "Data augmentation error: Could not generate exposure event $i")
      @assert(!(isnan(infectious_augmented[i]) && !isnan(obs.infectious[i])), "Data augmentation error: Could not generate infectious event $i")
      @assert(!(isnan(removed_augmented[i]) && !isnan(obs.removed[i])), "Data augmentation error: Could not generate exposure event $i")
    end
  end
  return SEIR_augmented(exposed_augmented, infectious_augmented, removed_augmented)
end


augment(ρ, Inf, network, obs, debug) = augment(ρ::Float64, network::array{Bool, 2}, obs::SEIR_observed, debug=false::Bool)


function augment(ρ::Float64, ν::Float64, obs::SEIR_observed, debug=false::Bool)
  """
  Augments surveilance data, organizes observations
  """
  exposed_augmented = fill(NaN, length(obs.infectious))
  infectious_augmented = fill(NaN, length(obs.infectious))
  removed_augmented = fill(NaN, length(obs.removed))
  for i = 1:length(obs.infectious)
    if !isnan(obs.infectious[i])
      if ν < Inf
        infectious_augmented[i] = obs.infectious[i] - rand(Exponential(1/ν))
      elseif ν == Inf
        infectious_augmented[i] = obs.infectious[i]
      end
      exposed_augmented[i] = infectious_augmented[i] - rand(Exponential(1/ρ))
      if !isnan(obs.removed[i])
        if ν < Inf
          removed_augmented[i] = obs.removed[i] - rand(Truncated(Exponential(1/ν), -Inf, obs.removed[i] - obs.infectious[i]))
        elseif ν == Inf
          removed_augmented[i] = obs.removed[i]
        end
      end
    end
  end
  if debug
    println("DATA AUGMENTATION")
    println("Augmented exposure times ($(sum(!isnan(exposed_augmented))) total): $(round(exposed_augmented,3))")
    println("Augmented infection times ($(sum(!isnan(removed_augmented))) total): $(round(infectious_augmented,3))")
    println("Augmented removal times ($(sum(!isnan(removed_augmented))) total): $(round(removed_augmented,3))")
    println("")
  end
  return SEIR_augmented(exposed_augmented, infectious_augmented, removed_augmented)
end


augment(ρ, Inf, obs, debug) = augment(ρ::Float64, obs::SEIR_observed, debug=false::Bool)


function logprior(priors::Priors, params::Vector{Float64})
  """
  Calculate the logprior from prior distributions
  """
  lprior = 0.
  for i = 1:length(params)
    lprior += logpdf(priors.(names(priors)[i]), params[i])
  end
  return lprior
end


function randprior(priors::Priors)
  """
  Randomly generate a parameter vector from specified priors
  """
  params = [rand(priors.(names(priors)[1]))]
  for i = 2:length(names(priors))
    push!(params, rand(priors.(names(priors)[i])))
  end
  return params
end


function propose_network(network_rates::Array{Float64, 2}, uniform=true::Bool, debug=false::Bool)
  """
  Propose a network based on network_rates
  """
  network = fill(false, size(network_rates))
  if uniform
    for i = 1:size(network_rates,2)
      rate_sum = sum(network_rates[:,i])
      if rate_sum > 0.
        network[sample(find(network_rates[:,i] .> 0.)), i] = true
      end
    end
  else
    for i = 1:size(network_rates,2)
      rate_sum = sum(network_rates[:,i])
      if rate_sum > 0.
        network[findfirst(rand(Multinomial(1, network_rates[:,i]/rate_sum))), i] = true
      end
    end
  end
  if debug
    println("NETWORK PROPOSAL")
    println("Individual exposure rate sums:")
    println("$(round(sum(network_rates, 1), 3))")
    println("Network proposal ($(sum(network)) infections total):")
    println("$(0 + network)")
    println("")
  end
  return network
end


function seq_distances(obs::SEIR_observed, aug::SEIR_augmented, infected::Vector{Int64}, network::Array{Bool, 2}, debug=false::Bool)
  """
  For a given transmission network, find the time between the pathogen sequences between every individuals i and j
  """
  if debug
    println("SEQUENCE DISTANCES")
    println("Pathogen sequences collected from $infected individuals")
    println("")
  end

  pathway = [infected[1]]

  if debug
    while pathway[end] > 0
      println("SEQUENCE DISTANCES")
      @assert(findfirst(network[:,pathway[end]]) > 0, "Network does not account for all exposure events")
      println("Adding individual $(findfirst(network[:,pathway[end]])-1) to individual $(infected[1])'s transmission pathway")
      push!(pathway, findfirst(network[:,pathway[end]])-1)
    end
  else
    while pathway[end] > 0
      push!(pathway, findfirst(network[:,pathway[end]])-1)
    end
  end
  pathways = Vector[pathway]

  for i = 2:length(infected)
    pathway = [infected[i]]
    if debug
      while pathway[end] > 0
        @assert(findfirst(network[:,pathway[end]]) > 0, "Network does not account for all exposure events")
        println("Adding individual $(findfirst(network[:,pathway[end]])-1) to individual $(infected[i])'s transmission pathway")
        push!(pathway, findfirst(network[:,pathway[end]])-1)
      end
    else
      while pathway[end] > 0
        push!(pathway, findfirst(network[:,pathway[end]])-1)
      end
    end
    push!(pathways, pathway)
  end

  seq_dist = fill(0., (size(network, 2), size(network, 2)))

  for i = 1:length(infected)
    if debug
      println("")
      println("SEQUENCE DISTANCES")
      println("Infection of individual $(infected[i]) observed at $(obs.infectious[infected[i]])")
      println("Infection pathway of individual $(infected[i]) is $(pathways[i])")
    end
    for j = 1:(i-1)
      k = 1
      while length(pathways[i]) > k && length(pathways[j]) > k && pathways[i][end - k] == pathways[j][end - k]
        k += 1
      end
      if debug
        println("Infection of individual $(infected[j]) observed at $(obs.infectious[infected[j]])")
        println("Infection pathway of individual $(infected[j]) is $(pathways[j])")
        if k == length(pathways[i]) || k == length(pathways[j])
          println("Linear infection pathway between individual $(infected[i]) and individual $(infected[j])")
        else
          println("Most recent common infection source of individuals $(infected[i]) and $(infected[j]) is $(pathways[i][end - k + 1])")
          println("The infection pathway of $(infected[i]) and $(infected[j]) diverged with $(pathways[i][end - k]) and $(pathways[j][end - k])")
        end
      end

      if k == length(pathways[i])
        seq_dist[infected[i],infected[j]] += obs.infectious[infected[j]] - aug.exposed[pathways[j][end - k]]
        seq_dist[infected[i],infected[j]] += abs(aug.exposed[pathways[j][end - k]] - obs.infectious[infected[i]])
      elseif k == length(pathways[j])
        seq_dist[infected[i],infected[j]] += obs.infectious[infected[i]] - aug.exposed[pathways[i][end - k]]
        seq_dist[infected[i],infected[j]] += abs(aug.exposed[pathways[i][end - k]] - obs.infectious[infected[j]])
      else
        seq_dist[infected[i],infected[j]] += obs.infectious[infected[i]] - aug.exposed[pathways[i][end - k]]
        seq_dist[infected[i],infected[j]] += obs.infectious[infected[j]] - aug.exposed[pathways[j][end - k]]
        seq_dist[infected[i],infected[j]] += abs(aug.exposed[pathways[j][end - k]] - aug.exposed[pathways[i][end - k]])
      end
    end
    if debug
      println("")
    end
  end

  return seq_dist += transpose(seq_dist)
end


function network_loglikelihood(obs::SEIR_observed, aug::SEIR_augmented, network::Array{Bool, 2}, substitution_matrix::Array{Float64, 2}, debug=false::Bool)
  """
  Loglikelihood for an entire transmission network
  """
  if debug
    @assert(size(substitution_matrix) == (4,4), "Network loglikelihood: Invalid substitution_matrix")
  end

  ll = 0.
  infected = find(isseq(obs.seq))
  seq_dist = seq_distances(obs, aug, infected, network, debug)
  for i = 1:length(infected)
    for j = 1:(i-1)
      ll += sum(log(expm(substitution_matrix*seq_dist[i,j]))[sub2ind((4,4), obs.seq[infected[i]], obs.seq[infected[j]])])
    end
  end
  return ll
end


function SEIR_loglikelihood(α::Float64, β::Float64, η::Float64, ρ::Float64, γ::Float64, aug::SEIR_augmented, obs::SEIR_observed, dist::Metric, debug=false::Bool)
  """
  Calculate the loglikelihood and return an exposure network array under specified parameters values and observations

  α, β: powerlaw exposure kernel parameters
  η: external pressure rate
  ρ: infectivity rate (1/mean latent period)
  γ: removal rate (1/mean infectious period)
  """
  # Initiate an exposure network
  network_rates = fill(0., (1 + length(obs.covariates), length(obs.covariates)))

  # Start with loglikelihood of 0
  ll = 0.

  # Create event timing array
  event_times = [aug.exposed aug.infectious aug.removed]

  # Find event order
  event_order = sortperm(event_times[:])

  # Create empty rate array
  rate_array = fill(0., (1 + length(obs.covariates) + 2, length(obs.covariates)))

  # First row is external pressure rate
  rate_array[1,:] = η

  # Loglikelihood of events that have been observed
  for i = 1:length(event_order)

    # Stop loglikelihood calculation after last event considered
    isnan(event_times[event_order[i]]) && break

    # Stop loglikelihood calculation anytime the loglikelihood goes to -Inf
    ll == -Inf && break

    # Convert linear index to an event tuple (individual, event type)
    id = ind2sub(size(event_times), event_order[i])

    # Don't consider likelilihood contribution of first event
    if i > 1
      # Find the "master rate" through the sum of rate_array
      master_rate = sum(rate_array)

      # loglikelihood of event time with master rate
      ll += logpdf(Exponential(1/master_rate), event_times[event_order[i]] - event_times[event_order[i-1]])

      # loglikelihood of the one event that did occur
      ll += log(sum(rate_array[:, id[1]]) / master_rate)
    end

    # Exposure event
    if id[2] == 1

      # Record exposure rates at time of exposure
      network_rates[:, id[1]] = copy(rate_array[1:(length(obs.covariates)+1), id[1]])

      # Update exposure rates
      rate_array[1:(1 + size(rate_array, 2)), id[1]] = 0.

      # Update infectivity rate
      rate_array[1 + size(rate_array, 2) + 1, id[1]] = ρ

    # Infectiousness event
    elseif id[2] == 2
      # Update infectivity rate
      rate_array[1 + size(rate_array, 2) + 1, id[1]] = 0.

      # Update removal rate
      rate_array[1 + size(rate_array, 2) + 2, id[1]] = γ

      # Update exposure rates for rest of susceptible population
      for j = 1:size(rate_array, 2)
        if j != id[1] && rate_array[1, j] > 0.
          rate_array[id[1] + 1, j] = α*evaluate(dist, obs.covariates[id[1]], obs.covariates[j])^-β
        end
      end

    # Removal event
    elseif id[2] == 3
      # Update removal rate
      rate_array[1 + size(rate_array, 2) + 2, id[1]] = 0.

      # Update exposure rates
      rate_array[id[1] + 1, :] = 0.
    end

    # Provide loop position when loglikelihood goes to -Inf when debugging
    if debug && ll == -Inf
      println("SEIR LOGLIKELIHOOD")
      if id[2] == 1
        println("Event $i (exposure of individual $(id[1])) caused loglikelihood to go to -Inf")
      elseif id[2] == 2
        println("Event $i (infection of individual $(id[1])) caused loglikelihood to go to -Inf")
      elseif id[2] == 3
        println("Event $i (removal of individual $(id[1])) caused loglikelihood to go to -Inf")
      end
      println("")
    end
  end
  return ll, network_rates
end


function initialize(ilm_priors::SEIR_priors, mutation_priors::JC69_priors, detection_priors::Lag_priors, obs::SEIR_observed, limit=500::Int, debug=false::Bool, dist=Euclidean())
  """
  Initiate an Trace object by sampling from specified prior distributions
  """
  ilm_params = randprior(ilm_priors)
  mutation_params = randprior(mutation_priors)
  detection_params = randprior(detection_priors)
  aug = augment(ilm_params[4], detection_params[1], obs, debug)
  lp1, network_rates = SEIR_loglikelihood(ilm_params[1], ilm_params[2], ilm_params[3], ilm_params[4], ilm_params[5], aug, obs, dist, debug)
  count = 1

  # Retry initialization until non-negative infinity loglikelihood
  while lp1 == -Inf && count < limit
    count += 1
    ilm_params = randprior(ilm_priors)
    mutation_params = randprior(mutation_priors)
    detection_params = randprior(detection_priors)
    aug = augment(ilm_params[4], detection_params[1], obs, debug)
    lp1, network_rates = SEIR_loglikelihood(ilm_params[1], ilm_params[2], ilm_params[3], ilm_params[4], ilm_params[5], aug, obs, dist, debug)
  end

  if count < limit
    lp1 += logprior(ilm_priors, ilm_params) + logprior(detection_priors, detection_params)
    network = propose_network(network_rates, false, debug)
    lp2 = network_loglikelihood(obs, aug, network, jc69([mutation_params[1]]), debug)
    lp2 += logprior(mutation_priors, mutation_params)

    while lp1 + lp2 == -Inf && count < limit
      count += 1
      lp1, network_rates = SEIR_loglikelihood(ilm_params[1], ilm_params[2], ilm_params[3], ilm_params[4], ilm_params[5], aug, obs, dist, debug)
      lp1 += logprior(ilm_priors, ilm_params) + logprior(detection_priors, detection_params)
      network = propose_network(network_rates, false, debug)
      lp2 = network_loglikelihood(obs, aug, network, jc69([mutation_params[1]]), debug)
      lp2 += logprior(mutation_priors, mutation_params)
    end

    if count < limit && lp1 + lp2 > -Inf
      println("INITIALIZATON")
      println("Successful on attempt $count (lp1 = $lp1, lp2 = $lp2)")
      return SEIR_trace([ilm_params[1]], [ilm_params[2]], [ilm_params[3]], [ilm_params[4]], [ilm_params[5]], [aug], Array[network_rates], Array[network], [lp1], [lp2]), Lag_trace([detection_params[1]]), JC69_trace([mutation_params[1]])
    else
      println("INITIALIZATON")
      println("Failed to initialize after $count attempts (lp1 = $lp1, lp2 = $lp2)")
    end
    println("")
  end
end


# function initialize(ilm_priors::SEIR_priors, mutation_priors::JC69_priors, obs::SEIR_observed, limit=500::Int, debug=false::Bool, dist=Euclidean())
#   """
#   Initiate an Trace object by sampling from specified prior distributions
#   """
#   ilm_params = randprior(ilm_priors)
#   mutation_params = randprior(mutation_priors)
#   aug = augment(ilm_params[4], obs)
#   ll, network = SEIR_loglikelihood(ilm_params[1], ilm_params[2], ilm_params[3], ilm_params[4], ilm_params[5], aug, obs, dist, debug)
#   ll += network_loglikelihood(obs, aug, network, jc69([mutation_params[1]]), debug)
#   count = 1

#   # Retry initialization until non-negative infinity loglikelihood
#   while ll == -Inf && count < limit
#     count += 1
#     ilm_params = randprior(ilm_priors)
#     mutation_params = randprior(mutation_priors)
#     aug = augment(ilm_params[4], obs)
#     ll, network = SEIR_loglikelihood(ilm_params[1], ilm_params[2], ilm_params[3], ilm_params[4], ilm_params[5], aug, obs, dist, debug)
#     ll += network_loglikelihood(obs, aug, network, jc69([mutation_params[1]]), debug)
#   end

#   if count < limit
#     print("Successfully initialized on attempt $count")
#     lp = ll + logprior(ilm_priors, ilm_params) + logprior(mutation_priors, mutation_params)
#       return SEIR_trace([ilm_params[1]], [ilm_params[2]], [ilm_params[3]], [ilm_params[4]], [ilm_params[5]], [aug], Array[network], [lp]), JC69_trace([mutation_params[1]])
#   else
#     print("Failed to initialize after $count attempts")
#   end
# end

# function initialize(ilm_priors::SEIR_priors, detection_priors::Lag_priors, obs::SEIR_observed, limit=500::Int, debug=false::Bool, dist=Euclidean())
#   """
#   Initiate an Trace object by sampling from specified prior distributions
#   """
#   ilm_params = randprior(ilm_priors)
#   detection_params = randprior(detection_priors)
#   aug = augment(ilm_params[4], detection_params[1], obs)
#   ll, network = SEIR_loglikelihood(ilm_params[1], ilm_params[2], ilm_params[3], ilm_params[4], ilm_params[5], aug, obs, dist, debug)
#   count = 1

#   # Retry initialization until non-negative infinity loglikelihood
#   while ll == -Inf && count < limit
#     count += 1
#     ilm_params = randprior(ilm_priors)
#     detection_params = randprior(detection_priors)
#     aug = augment(ilm_params[4], detection_params[1], obs)
#     ll, network = SEIR_loglikelihood(ilm_params[1], ilm_params[2], ilm_params[3], ilm_params[4], ilm_params[5], aug, obs, dist, debug)
#   end

#   if count < limit
#     print("Successfully initialized on attempt $count")
#     lp = ll + logprior(ilm_priors, ilm_params) + logprior(detection_priors, detection_params)
#       return SEIR_trace([ilm_params[1]], [ilm_params[2]], [ilm_params[3]], [ilm_params[4]], [ilm_params[5]], [aug], Array[network], [lp]), Lag_trace([detection_params[1]])
#   else
#     print("Failed to initialize after $count attempts")
#   end
# end


# function initialize(ilm_priors::SEIR_priors, obs::SEIR_observed, limit=500::Int, debug=false::Bool, dist=Euclidean())
#   """
#   Initiate an Trace object by sampling from specified prior distributions
#   """
#   ilm_params = randprior(ilm_priors)
#   aug = augment(ilm_params[4], obs)
#   ll, network_rates = SEIR_loglikelihood(ilm_params[1], ilm_params[2], ilm_params[3], ilm_params[4], ilm_params[5], aug, obs, dist, debug)
#   count = 1

#   # Retry initialization until non-negative infinity loglikelihood
#   while ll == -Inf && count < limit
#     count += 1
#     ilm_params = randprior(ilm_priors)
#     aug = augment(ilm_params[4], obs)
#     ll, network_rates = SEIR_loglikelihood(ilm_params[1], ilm_params[2], ilm_params[3], ilm_params[4], ilm_params[5], aug, obs, dist, debug)
#   end

#   if count < limit
#     print("Successfully initialized on attempt $count")
#     network = propose_network(network_rates, false)
#     lp = ll + logprior(ilm_priors, ilm_params)
#       return SEIR_trace([ilm_params[1]], [ilm_params[2]], [ilm_params[3]], [ilm_params[4]], [ilm_params[5]], [aug], Array[network_rates], Array[network], [lp], [0.])
#   else
#     print("Failed to initialize after $count attempts")
#   end
# end


function MCMC(n::Int64,
              transition_cov::Array{Float64},
              ilm_trace::SEIR_trace,
              detection_trace::Lag_trace,
              mutation_trace::JC69_trace,
              ilm_priors::SEIR_priors,
              detection_priors::Lag_priors,
              mutation_priors::JC69_priors,
              obs::SEIR_observed,
              debug=false::Bool,
              progress=true::Bool,
              dist=Euclidean())
  """
  Performs `n` data-augmented metropolis hastings within Gibbs MCMC iterations. Initiates a single chain by sampling from prior distribution
  """

  @assert(size(transition_cov) == (7,7), "transition_cov must be a 7x7 matrix")

  for i = 1:n

    # Create and incremenet progress bar
    if progress
      if i == 1
        progressbar = Progress(n, 5, "Performing $n MCMC iterations...", 30)
      else
        next!(progressbar)
      end
    end

    # Step 1a: Metropolis-Hastings proposal
    # Only generate valid proposals
    step = rand(MvNormal(transition_cov))
    ilm_proposal = [ilm_trace.α[end], ilm_trace.β[end], ilm_trace.η[end], ilm_trace.ρ[end], ilm_trace.γ[end]] .+ step[1:5]
    detection_proposal = [detection_trace.ν[end]] .+ step[6]
    mutation_proposal = [mutation_trace.λ[end]] .+ step[7]

    lp1_proposal = logprior(ilm_priors, ilm_proposal)
    lp1_proposal += logprior(detection_priors, detection_proposal)
    lp2_proposal = logprior(mutation_priors, mutation_proposal)

    while lp1_proposal + lp2_proposal == -Inf
      step = rand(MvNormal(transition_cov))
      ilm_proposal = [ilm_trace.α[end], ilm_trace.β[end], ilm_trace.η[end], ilm_trace.ρ[end], ilm_trace.γ[end]] .+ step[1:5]
      detection_proposal = [detection_trace.ν[end]] .+ step[6]
      mutation_proposal = [mutation_trace.λ[end]] .+ step[7]
      lp1_proposal = logprior(ilm_priors, ilm_proposal)
      lp1_proposal += logprior(detection_priors, detection_proposal)
      lp2_proposal = logprior(mutation_priors, mutation_proposal)
    end

    # Step 1b: Gibbs within Metropolis data augmentation
    aug_proposal = augment(ilm_proposal[4], detection_proposal[1], ilm_trace.network[end], obs, debug)

    # Step 1c: loglikelihood calculation for Metropolis-Hastings step
    ll_proposal, network_rates_proposal = SEIR_loglikelihood(ilm_proposal[1], ilm_proposal[2], ilm_proposal[3], ilm_proposal[4], ilm_proposal[5], aug_proposal, obs, dist)
    lp1_proposal += ll_proposal

    # Step 1d: accept/reject based on logposterior comparison
    reject = true
    if lp1_proposal > ilm_trace.logposterior_1[end]
      reject = false
    elseif exp(lp1_proposal - ilm_trace.logposterior_1[end]) > rand()
      reject = false
    end

    if debug
      if lp1_proposal > ilm_trace.logposterior_1[end]
        println("MCMC")
        println("Accepted ILM proposal ($lp1_proposal > $(ilm_trace.logposterior_1[end])) on $(i)th iteration")
      elseif reject == false
        println("MCMC")
        println("Accepted ILM proposal (with probability $(exp(lp1_proposal - ilm_trace.logposterior_1[end])) on $(i)th iteration)")
      else
        println("MCMC")
        println("Rejected ILM proposal (with probability $(1-exp(lp1_proposal - ilm_trace.logposterior_1[end])) on $(i)th iteration")
      end
      println("")
    end

    if reject
      ilm_proposal[1] = ilm_trace.α[end]
      ilm_proposal[2] = ilm_trace.β[end]
      ilm_proposal[3] = ilm_trace.η[end]
      ilm_proposal[4] = ilm_trace.ρ[end]
      ilm_proposal[5] = ilm_trace.γ[end]
      aug_proposal = ilm_trace.aug[end]
      network_rates_proposal = ilm_trace.network_rates[end]
      lp1_proposal = ilm_trace.logposterior_1[end]
      detection_proposal = detection_trace.ν[end]
    end

    # Step 1e: Update chain
    push!(ilm_trace.α, ilm_proposal[1])
    push!(ilm_trace.β, ilm_proposal[2])
    push!(ilm_trace.η, ilm_proposal[3])
    push!(ilm_trace.ρ, ilm_proposal[4])
    push!(ilm_trace.γ, ilm_proposal[5])
    push!(ilm_trace.aug, aug_proposal)
    push!(ilm_trace.network_rates, network_rates_proposal)
    push!(ilm_trace.logposterior_1, lp1_proposal)
    push!(detection_trace.ν, detection_proposal[1])

    # Step 3a: Independence sampling of network
    network_proposal = propose_network(ilm_trace.network_rates[end], false, debug)
    lp2_proposal += network_loglikelihood(obs, ilm_trace.aug[end], network_proposal, jc69([mutation_proposal[1]]), debug)

    lp2 = logprior(mutation_priors, [mutation_trace.λ[end]])
    lp2 += network_loglikelihood(obs, ilm_trace.aug[end], ilm_trace.network[end], jc69([mutation_trace.λ[end]]), debug)

    # Step 3b: Accept/reject proposal
    reject = true
    if lp2_proposal > lp2
      reject = false
    elseif exp(lp2_proposal - lp2) > rand()
      reject = false
    end

    if reject
      lp2_proposal = lp2
      network_proposal = ilm_trace.network[end]
      mutation_proposal = mutation_trace.λ[end]
    end

    if debug
      if lp2_proposal > lp2
        println("MCMC")
        println("Accepted ILM proposal ($lp2_proposal > $lp2) on $(i)th iteration")
      elseif reject == false
        println("MCMC")
        println("Accepted ILM proposal (with probability $(exp(lp2_proposal - lp2)) on $(i)th iteration)")
      else
        println("MCMC")
        println("Rejected ILM proposal (with probability $(1-exp(lp2_proposal - lp2)) on $(i)th iteration")
      end
      println("")
    end

    push!(ilm_trace.logposterior_2, lp2_proposal)
    push!(ilm_trace.network, network_proposal)
    push!(mutation_trace.λ, mutation_proposal[1])
  end

  return ilm_trace, detection_trace, mutation_trace
end


function MCMC(n::Int64,
              transition_cov::Array{Float64},
              ilm_priors::SEIR_priors,
              detection_priors::Lag_priors,
              mutation_priors::JC69_priors,
              obs::SEIR_observed,
              debug=false::Bool,
              progress=true::Bool,
              dist=Euclidean(),
              init_limit=500::Int64)
  """
  Performs `n` data-augmented metropolis hastings within Gibbs MCMC iterations. Initiates a single chain by sampling from prior distribution
  """

  @assert(size(transition_cov) == (7,7), "transition_cov must be a 7x7 matrix")

  ilm_trace, detection_trace, mutation_trace = initialize(ilm_priors, mutation_priors, detection_priors, obs, init_limit, debug, dist)

  return MCMC(n,
              transition_cov,
              ilm_trace,
              detection_trace,
              mutation_trace,
              ilm_priors,
              detection_priors,
              mutation_priors,
              obs,
              debug,
              progress,
              dist)

  end


function MCMC(n::Int64,
              ilm_priors::SEIR_priors,
              detection_priors::Lag_priors,
              mutation_priors::JC69_priors,
              obs::SEIR_observed,
              debug=false::Bool,
              progress=true::Bool,
              dist=Euclidean(),
              init_limit=500::Int64)
  """
  Performs `n` data-augmented metropolis hastings within Gibbs MCMC iterations. Initiates a single chain by sampling from prior distribution
  """

  transition_cov = diagm([var(ilm_priors.α),
                          var(ilm_priors.β),
                          var(ilm_priors.η),
                          var(ilm_priors.ρ),
                          var(ilm_priors.γ),
                          var(detection_priors.ν),
                          var(mutation_priors.λ)])

  return MCMC(n,
              transition_cov,
              ilm_priors,
              detection_priors,
              mutation_priors,
              obs,
              debug,
              progress,
              dist,
              init_limit)
end
