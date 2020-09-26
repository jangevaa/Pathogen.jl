function _outdegree(x::TN) where {TN <: Union{TNDistribution, TransmissionNetwork}}
  infected = [x.external[i] > 0 || any(x.internal[:, i] .> 0) for i = 1:individuals(x)]
  return sum(x.internal[infected, infected], dims=2)
end

@recipe function plot(x::TN) where {TN <: Union{TNDistribution, TransmissionNetwork}}
  seriestype := :histogram
  legend --> :none
  xguide --> "Out degree"
  yguide --> "Frequency"
  od = _outdegree(x)
  bins --> 0:ceil(max(od))
  od
end