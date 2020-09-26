struct Population
  risks::AbstractDataFrame
  distances::Union{Nothing, AbstractArray}
end

Population(risks) = Population(risks, nothing)

function individuals(x::Population)
  return size(x.risks, 1)
end

function Base.show(io::IO, x::Population)
  return print(io, "Population object (n=$(individuals(x)))")
end