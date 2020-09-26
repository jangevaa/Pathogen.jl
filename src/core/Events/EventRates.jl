struct EventRates{T <: DiseaseStateSequence}
  exposure::Union{Nothing, Vector{Float64}}
  infection::Vector{Float64}
  removal::Union{Nothing, Vector{Float64}}

  function EventRates{T}(n::Int64) where T <: SEIR
    return new{T}(fill(0.0, n), fill(0.0, n), fill(0.0, n))
  end

  function EventRates{T}(n::Int64) where T <: SEI
    return new{T}(fill(0.0, n), fill(0.0, n), nothing)
  end
  
  function EventRates{T}(n::Int64) where T <: SIR
    return new{T}(nothing, fill(0.0, n), fill(0.0, n))
  end
  
  function EventRates{T}(n::Int64) where T <: SI
    return new{T}(nothing, fill(0.0, n), nothing)
  end
end

function individuals(x::EventRates{T}) where {
                     T <: DiseaseStateSequence}
  return length(x.infection)
end

function Base.getindex(x::EventRates{T}, new_state::DiseaseState) where T <: DiseaseStateSequence
  if new_state == State_E
    return x.exposure
  elseif new_state == State_I
    return x.infection
  elseif new_state == State_R
    return x.removal
  else
    @error "Unrecognized indexing disease state"
  end
end

function Base.getindex(x::EventRates{T}, new_states::DiseaseStates) where T <: DiseaseStateSequence
  y = x[new_states[1]]
  for i=2:length(new_states)
    y = hcat(y, x[new_states[i]])
  end
  return y
end

function Base.sum(x::EventRates{T}) where T <: DiseaseStateSequence
  return sum([sum(x[i]) for i in convert(DiseaseStates, T)[2:end]])
end
