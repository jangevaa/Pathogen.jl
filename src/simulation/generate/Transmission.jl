function generate(::Type{Transmission},
                  tr::TransmissionRates,
                  event::Event{T}) where T <: DiseaseStateSequence
  id = event.individual
  if _new_transmission(event)
    external_or_internal = Weights([tr.external[id]; sum(tr.internal[:,id])])
    if sum(external_or_internal) == 0.0
      @error "All transmission rates = 0.0, No transmission can be generated"
      return NoTransmission()
    elseif sample([true; false], external_or_internal)
      @debug "Exogenous tranmission generated"
      return ExogenousTransmission(id)
    else
      source = sample(1:individuals(tr), Weights(tr.internal[:, id]))
      @debug "Endogenous transmission generated (source id = $source)"
      return EndogenousTransmission(id, source)
    end
  else
    @debug "No transmission generated"
    return NoTransmission()
  end
end