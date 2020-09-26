abstract type AbstractRisk{T <: DiseaseStateSequence} end


function Base.copy(x::AbstractRisk{SEIR})
  return AbstractRisk{SEIR}(copy(x.sparks),
                            copy(x.susceptibility),
                            copy(x.infectivity),
                            copy(x.transmissibility),
                            copy(x.latency),
                            copy(x.removal))
end

function Base.copy(x::AbstractRisk{SEI})
  return AbstractRisk{SEI}(copy(x.sparks),
                           copy(x.susceptibility),
                           copy(x.infectivity),
                           copy(x.transmissibility),
                           copy(x.latency))
end

function Base.copy(x::AbstractRisk{SIR})
  return AbstractRisk{SIR}(copy(x.sparks),
                           copy(x.susceptibility),
                           copy(x.infectivity),
                           copy(x.transmissibility),
                           copy(x.removal))
end

function Base.copy(x::AbstractRisk{SI})
  return AbstractRisk{SI}(copy(x.sparks),
                          copy(x.susceptibility),
                          copy(x.infectivity),
                          copy(x.transmissibility))
end


function _indices(x::AbstractRisk{T}; zeros::Bool=true, cumulative::Bool=true) where T <: DiseaseStateSequence
  indices = [length(x.sparks)
             length(x.susceptibility)
             length(x.infectivity)
             length(x.transmissibility)]
  if State_E ∈ T
    push!(indices, length(x.latency))
  elseif zeros
    push!(indices, 0)
  end
  if State_R ∈ T
    push!(indices, length(x.removal))
  elseif zeros
    push!(indices, 0)
  end
  if cumulative
    return cumsum(indices)
  else
    return indices
  end
end


function Base.length(x::AbstractRisk{T}) where T <: DiseaseStateSequence
  return _indices(x)[end]
end


function Base.getindex(x::AbstractRisk{T},
                       i::Int64) where T <: DiseaseStateSequence
  indices = _indices(x, zeros = true)
  riskfunc = findfirst(i .<= indices)
  return getfield(x, riskfunc)[end - (indices[riskfunc] - i)]
end
