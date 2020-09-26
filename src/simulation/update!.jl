function update!(
  s::Simulation{T},
  event::AbstractEvent) where {
    T <: DiseaseStateSequence}
  if _time(event) < s.time
    error("Time of a new event must be >= that of the previous event")
  elseif _time(event) < Inf
    update!(s.events, event)
    update!(s.disease_states, event)
    update!(s.transmission_network,
            generate(Transmission,
                     s.transmission_rates,
                     event))
    update!(s.transmission_rates,
            s.event_rates,
            event,
            s.disease_states,
            s.population,
            s.risk_functions,
            s.risk_parameters)
    # Count the iteration as having occurred.
    s.iterations += 1
  end
  # Update simulation time
  s.time = _time(event)
  # Return updated Simulation object
  return s
end