function _ids_by_state(events::Events{T}, state::DiseaseState, time::Float64) where T <: EpidemicModel
  if time < 0.0
    error("Time must be ≥ 0.0")
  end
  if state == _state_progressions[T][1]
    nextstate = _state_progressions[T][2]
    return findall((events[nextstate] .> time) .| isnan.(events[nextstate]))
  elseif state ∈ _state_progressions[T][2:end-1]
    nextstate = _state_progressions[T][findfirst(_state_progressions[T] .== state) + 1]
    return findall((events[state] .≤ time) .& ((events[nextstate] .> time) .| isnan.(events[nextstate])))
  elseif state == _state_progressions[T][end]
    return findall(events[state] .≤ time)
  else
    error("Invalid state specified")
  end
end
