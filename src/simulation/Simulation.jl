mutable struct Simulation{T <: DiseaseStateSequence}
  time::Float64
  iterations::Int64
  population::Population
  risk_functions::RiskFunctions{T}
  risk_parameters::RiskParameters{T}
  disease_states::Vector{DiseaseState}
  transmission_rates::TransmissionRates
  event_rates::EventRates{T}
  events::Events{T}
  transmission_network::TransmissionNetwork

  function Simulation(pop::Population,
                      rf::RiskFunctions{T},
                      rp::RiskParameters{T}) where T <: DiseaseStateSequence
    states = fill(State_S, individuals(pop))
    tr = initialize(TransmissionRates, states, pop, rf, rp)
    rates = initialize(EventRates, tr, states, pop, rf, rp)
    events = Events{T}(individuals(pop))
    net = TransmissionNetwork(states)
    return new{T}(0.0, 0, pop, rf, rp, states, tr, rates, events, net)
  end

  function Simulation(pop::Population,
                      states::Vector{DiseaseState},
                      rf::RiskFunctions{T},
                      rp::RiskParameters{T};
                      time::Float64 = 0.0,
                      skip_checks::Bool=false) where T <: DiseaseStateSequence
    @debug "Initializing $T Simulation with the following starting states:" states
    if !skip_checks
      if length(states) != individuals(pop)
        @error "Length of initial disease state vector must match number of individuals"
      elseif !all(in.(states, Ref(convert(DiseaseStates, T))))
        @error "All states in initial disease state vector must be valid within specified epidemic model"
      end
    end
    tr = initialize(TransmissionRates, states, pop, rf, rp)
    rates = initialize(EventRates, tr, states, pop, rf, rp)
    events = Events{T}(states)
    net = TransmissionNetwork(states)
    return new{T}(time, 0, pop, rf, rp, copy(states), tr, rates, events, net)
  end
end

function Base.show(io::IO, x::Simulation{T}) where T <: DiseaseStateSequence
  y = "$T epidemic simulation @ time = $(round(x.time, digits=2))\n"
  for s in convert(DiseaseStates, T)
    y *= "\n$s = $(sum(x.disease_states .== Ref(s)))"
  end
  return print(io, y)
end
