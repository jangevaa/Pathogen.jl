struct RiskPriors{T} <: AbstractRisk{T}
  sparks::Union{Nothing, Vector{UnivariateDistribution}}
  susceptibility::Union{Nothing, Vector{UnivariateDistribution}}
  infectivity::Union{Nothing, Vector{UnivariateDistribution}}
  transmissibility::Union{Nothing, Vector{UnivariateDistribution}}
  latency::Union{Nothing, Vector{UnivariateDistribution}}
  removal::Union{Nothing, Vector{UnivariateDistribution}}

  RiskPriors{SEIR}(pϵ, ps, pκ, pt, pl, pr) = new{SEIR}(pϵ, ps, pκ, pt, pl, pr)
  RiskPriors{SEI}(pϵ, ps, pκ, pt, pl) = new{SEI}(pϵ, ps, pκ, pt, pl, nothing)
  RiskPriors{SIR}(pϵ, ps, pκ, pt, pr) = new{SIR}(pϵ, ps, pκ, pt, nothing, pr)
  RiskPriors{SI}(pϵ, ps, pκ, pt) = new{SI}(pϵ, ps, pκ, pt, nothing, nothing)
end


function Base.show(io::IO, x::RiskPriors{T}) where T <: DiseaseStateSequence
  return print(io, "$T model risk function priors")
end
