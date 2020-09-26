function generate(::Type{Event},
                  rates::EventRates{T},
                  time::Float64) where T <: DiseaseStateSequence
  totals = Weights([sum(rates[state]) for state in convert(DiseaseStates, T)[2:end]])
  if sum(totals) == Inf
    new_state = sample(convert(DiseaseStates, T)[2:end], totals)
    id = sample(1:individuals(rates), Weights(rates[new_state]))
    return Event{T}(time, id, new_state)
  elseif sum(totals) == 0.0
    return NoEvent()
  else
    # Generate new state
    new_state = sample(convert(DiseaseStates, T)[2:end], totals)
    # Generate event individual
    id = sample(1:individuals(rates), Weights(rates[new_state]))
    return Event{T}(time + rand(Exponential(1.0 / sum(totals))), id, new_state)
  end
end