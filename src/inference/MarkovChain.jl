mutable struct MarkovChain{T <: DiseaseStateSequence}
  iterations::Int64
  events::Vector{Events{T}}
  transmission_network::Vector{TransmissionNetwork}
  risk_parameters::Vector{RiskParameters{T}}
  log_posterior::Vector{Float64}
  cov::OnlineStats.CovMatrix

  function MarkovChain(events::Events{T},
                       tn::TransmissionNetwork,
                       rp::RiskParameters{T},
                       lp::Float64) where  T <: DiseaseStateSequence
    return new{T}(0, [events], [tn], [rp], [lp], OnlineStats.CovMatrix())
  end

  function MarkovChain{T}() where  T <: DiseaseStateSequence
    return new{T}(0, Events{T}[], TransmissionNetwork[], RiskParameters{T}[], Float64[], OnlineStats.CovMatrix())
  end
end

function Base.length(x::MarkovChain{T}) where T <: DiseaseStateSequence
  return x.iterations
end

function Base.show(io::IO, x::MarkovChain{T}) where T <: DiseaseStateSequence
  return print(io, "$T model Markov chain (iterations = $(length(x)))")
end

function TransmissionNetworkDistribution(x::MarkovChain; burnin::Int64=0, thin::Int64=1)
  return TNDistribution(x.transmission_network[1+burnin:thin:end])
end