abstract type AbstractEvent end

struct NoEvent <: AbstractEvent
end

struct Event{T <: DiseaseStateSequence} <: AbstractEvent
  time::Float64
  individual::Int64
  new_state::DiseaseState

  function Event{T}(time::Float64, id::Int64, new_state::DiseaseState) where T <: DiseaseStateSequence
    if new_state in convert(DiseaseStates, T)
      return new{T}(time, id, new_state)
    else
      @error "Invalid disease state provided for $(T) Epidemic Models."
    end
  end

  function Event(new_time::Float64, event::Event{T}) where T <: DiseaseStateSequence
    return new{T}(new_time, event.individual, event.new_state)
  end
end

function _time(x::Event{T}) where T <: DiseaseStateSequence
  return x.time
end

function _time(x::NoEvent)
  return Inf
end

# function Base.copy(x::Event{T}) where T <: DiseaseStateSequence
#   return Event{T}(copy(x.time), copy(x.id), copy(x.new_state))
# end

function _new_transmission(e::Event{T}) where T <: Union{SEIR, SEI}
  return e.new_state == State_E
end

function _new_transmission(e::Event{T}) where T <: Union{SIR, SI}
  return e.new_state == State_I
end

function _individual(x::Event{T}) where T <: DiseaseStateSequence
  return x.individual
end

function _new_state(x::Event{T}) where T <: DiseaseStateSequence
  return x.new_state
end

function Base.show(io::IO, x::Event{T}) where T <: DiseaseStateSequence
  return print(io, "Transition of individual $(x.individual) into $(x.new_state) state at t = $(x.time)")
end
