_tupleme(x::Float64) = x < 0.0 ? error("Max delay must be positive") : (0.0, x)
_tupleme(x::Tuple{Float64, Float64}) = x[1] >= x[2] ? error("Specify delay as (min, max)") : x

struct EventExtents{T <: DiseaseStateSequence}
  exposure::Union{Nothing, Tuple{Float64, Float64}}
  infection::Union{Nothing, Tuple{Float64, Float64}}
  removal::Union{Nothing, Tuple{Float64, Float64}}
  EventExtents{SEIR}(e, i, r) = new{SEIR}(_tupleme(e), _tupleme(i), _tupleme(r))
  EventExtents{SEI}(e, i)     = new{SEI}(_tupleme(e), _tupleme(i), nothing)
  EventExtents{SIR}(i, r)     = new{SIR}(nothing, _tupleme(i), _tupleme(r))
  EventExtents{SI}(i)         = new{SI}(nothing, _tupleme(i), nothing)
end

function Base.show(io::IO, x::EventExtents{T}) where T <: DiseaseStateSequence
  print(io, "$T model event extents")
end
