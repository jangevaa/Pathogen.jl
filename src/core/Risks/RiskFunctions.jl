struct RiskFunctions{T <: DiseaseStateSequence}
  sparks::Union{Nothing, Function}
  susceptibility::Union{Nothing, Function}
  infectivity::Union{Nothing, Function}
  transmissibility::Union{Nothing, Function}
  latency::Union{Nothing, Function}
  removal::Union{Nothing, Function}

  RiskFunctions{SEIR}(ϵ, Ωs, κ, Ωt, Ωl, Ωr) = new{SEIR}(ϵ, Ωs, κ, Ωt, Ωl, Ωr)
  RiskFunctions{SEI}(ϵ, Ωs, κ, Ωt, Ωl) = new{SEI}(ϵ, Ωs, κ, Ωt, Ωl, nothing)
  RiskFunctions{SIR}(ϵ, Ωs, κ, Ωt, Ωr) = new{SIR}(ϵ, Ωs, κ, Ωt, nothing, Ωr)
  RiskFunctions{SI}(ϵ, Ωs, κ, Ωt) = new{SI}(ϵ, Ωs, κ, Ωt, nothing, nothing)
end


function Base.show(io::IO, x::RiskFunctions{T}) where T <: DiseaseStateSequence
  return print(io, "$T model risk functions")
end
