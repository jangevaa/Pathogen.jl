function _count_by_state(
  events::EventObservations{T, M},
  state::DiseaseState,
  time::Float64) where {
    T <: DiseaseStateSequence,
    M <: ILM}
  n_ids = 0
  if state == State_I && State_R ∈ T
  # E/I at or before time and I/R after time or never
    for i = 1:individuals(events)
      n_ids += events[state][i] <= time && ((events[State_R][i] > time) || (events[State_R][i] === NaN))
    end
  else
    for i = 1:individuals(events)
      n_ids += events[state][i] <= time
    end
  end
  @debug "$n_ids individual(s) observered in state $state at t = $time"
  return n_ids
end

function _obs_curve(
  events::EventObservations{T, M},
  state::DiseaseState,
  min::Float64,
  max::Float64) where {
    T <: DiseaseStateSequence,
    M <: ILM}
  if min >= max
    @error "Minimum time must be less than maximum time"
  end
  local times
  if state == State_I && State_R ∈ T
    times = events[[State_I; State_R]][:]
  else
    times = events[state]
  end
  times = times[Ref(min) .< times .< Ref(max)]
  sort!(times)
  insert!(times, 1, min)
  push!(times, max)
  counts = _count_by_state.(Ref(events), Ref(state), times)
  return times, counts
end

@recipe function f(
  events::EventObservations{T, M},
  state::DiseaseState,
  min::Float64,
  max::Float64) where {
    T <: DiseaseStateSequence,
    M <: ILM}
  xguide --> "Time"
  yguide --> "N"
  linewidth --> 2.0
  linecolor --> _state_color(state)
  label --> ""
  seriestype --> :steppost
  _obs_curve(events, state, min, max)
end

@recipe function f(
  events::EventObservations{T, M},
  state::DiseaseState) where {
    T <: DiseaseStateSequence,
    M <: ILM}
  events, state, 0.0, maximum(events)
end

@recipe function f(
  events::EventObservations{T, M},
  min::Float64,
  max::Float64) where {
    T <: DiseaseStateSequence,
    M <: ILM}
  @series begin
    linecolor --> _state_color(State_I)
    label --> "I"
    events, State_I, min, max
  end
  if State_R ∈ T
    @series begin
      linecolor --> _state_color(State_R)
      label --> "R"
      events, State_R, min, max
    end
  end
end

@recipe function f(
  events::EventObservations{T, M}) where {
    T <: DiseaseStateSequence,
    M <: ILM}
  events, minimum(events), maximum(events)
end
