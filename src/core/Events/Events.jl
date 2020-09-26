struct Events{T <: DiseaseStateSequence}
  exposure::Union{Nothing, Vector{Float64}}
  infection::Vector{Float64}
  removal::Union{Nothing, Vector{Float64}}

  function Events{T}(e::V, i::V, r::V) where {T <: SEIR, V <: Vector{Float64}}
    if length(unique((length.([e; i; r])))) != 1
      @error "Length of event time vectors must be equal"
    end
    return new{T}(e, i, r)
  end

  function Events{T}(n::Integer) where T <: SEIR
    return new{T}(fill(NaN, n), fill(NaN, n), fill(NaN, n))
  end

  function Events{T}(e::V, i::V) where {T <: SEI, V <: Vector{Float64}}
    if length(unique((length.([e; i])))) != 1
      @error "Length of event time vectors must be equal"
    end
    return new{T}(e, i, nothing)
  end

  function Events{T}(n::Integer) where T <: SEI
    return new{T}(fill(NaN, n), fill(NaN, n), nothing)
  end

  function Events{T}(i::V, r::V) where {T <: SIR, V <: Vector{Float64}}
    if length(unique((length.([i; r])))) != 1
      @error "Length of event time vectors must be equal"
    end
    return new{T}(nothing, i, r)
  end

  function Events{T}(n::Integer) where T <: SIR
    return new{T}(nothing, fill(NaN, n), fill(NaN, n))
  end

  function Events{T}(i::V) where {T <: SI, V <: Vector{Float64}}
    return new{T}(nothing, i, nothing)
  end

  function Events{T}(n::Integer) where T <: SI
    return new{T}(nothing, fill(NaN, n), nothing)
  end
end

function Events{T}(a::Array{Float64,2}) where T <: DiseaseStateSequence
  if size(a, 2) == 3
    return Events{T}(a[:,1], a[:,2], a[:,3])
  elseif size(a, 2) == 2
    return Events{T}(a[:,1], a[:,2])
  elseif size(a, 2) == 1
    return Events{T}(a[:,1])
  else
    @error "Invalid array size for construction of an $(Events{T}) object"
  end
end

function Events{T}(x::DiseaseStates) where T <: DiseaseStateSequence
  events = Events{T}(length(x))
  for i = 1:length(x)
    if x[i] in convert(DiseaseStates, T)
      for j = convert(DiseaseStates, T)[2:findfirst(Ref(x[i]) .== convert(DiseaseStates, T))]
        events[j][i] = -Inf
      end
    else
      @error "Invalid initial state for individual $i"
    end
  end
  return events
end

function individuals(x::Events{M}) where{
                     M <: DiseaseStateSequence}
  return length(x.infection)
end

function Base.show(io::IO, x::Events{T}) where T <: DiseaseStateSequence
  return print(io, "$T model event times (n=$(individuals(x)))")
end

function Base.getindex(x::Events{T}, new_state::DiseaseState) where T <: DiseaseStateSequence
  if new_state == State_E
    return x.exposure
  elseif new_state == State_I
    return x.infection
  elseif new_state == State_R
    return x.removal
  else
    error("Unsupported indexing disease state")
  end
end

function Base.getindex(x::Events{T}, states::DiseaseStates) where T <: DiseaseStateSequence
  y = x[states[1]]
  for i = 2:length(states)
    y = hcat(y, x[states[i]])
  end
  return y
end

function Base.convert(::Type{Array{Float64, 2}}, x::Events{T}) where T <: DiseaseStateSequence
  return x[convert(DiseaseStates, T)[2:end]]
end

function Base.convert(::Type{Vector{Float64}}, x::Events{T}) where T <: DiseaseStateSequence
  return x[convert(DiseaseStates, T)[2:end]][:]
end

function Base.convert(::Type{Array{Float64, 2}}, x::Array{Events{T}, 1}) where T <: DiseaseStateSequence
  y = convert(Vector{Float64}, x[1])'
  for i = 2:length(x)
    y = vcat(y, convert(Vector{Float64}, x[i])')
  end
  return y
end

function Base.minimum(x::Events{T}) where T <: DiseaseStateSequence
  y = convert(Array{Float64, 2}, x)
  return minimum(y[(y .!== NaN) .& (y .> -Inf)])
end

function Base.maximum(x::Events{T}) where T <: DiseaseStateSequence
  y = convert(Array{Float64, 2}, x)
  return maximum(y[(y .!== NaN) .& (y .> -Inf)])
end

function mean(x::Vector{Events{T}}) where T <: DiseaseStateSequence
  return Events{T}(mean([convert(Array{Float64, 2}, i) for i in x]))
end