mutable struct MCMC{T <: EpidemicModel}
  observations::EventObservations{T}
  event_extents::EventExtents{T}
  population::DataFrame
  risk_functions::RiskFunctions{T}
  risk_priors::RiskPriors{T}
  traces::Vector{Trace{T}}

  function MCMC(obs::EventObservations{T},
                ee::EventExtents{T},
                pop::DataFrame,
                rf::RiskFunctions{T},
                rp::RiskPriors{T}) where T <: EpidemicModel
    return new{T}(obs, ee, pop, rf, rp, Trace{T}[])
  end
end