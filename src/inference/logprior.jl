function logprior(
  rparams::RiskParameters{S},
  rpriors::RiskPriors{S}) where {
  S <: DiseaseStateSequence}
  lprior = 0.0
  for i in 1:length(rpriors.sparks)
    lprior == -Inf && break
    lprior += Distributions.logpdf(rpriors.sparks[i], rparams.sparks[i])
  end
  for i in 1:length(rpriors.susceptibility)
    lprior == -Inf && break
    lprior += Distributions.logpdf(rpriors.susceptibility[i], rparams.susceptibility[i])
  end
  for i in 1:length(rpriors.infectivity)
    lprior == -Inf && break
    lprior += Distributions.logpdf(rpriors.infectivity[i], rparams.infectivity[i])
  end
  for i in 1:length(rpriors.transmissibility)
    lprior == -Inf && break
    lprior += Distributions.logpdf(rpriors.transmissibility[i], rparams.transmissibility[i])
  end
  if State_E ∈ S
    for i in 1:length(rpriors.latency)
      lprior == -Inf && break
      lprior += Distributions.logpdf(rpriors.latency[i], rparams.latency[i])
    end
  end
  if State_R ∈ S
    for i in 1:length(rpriors.removal)
      lprior == -Inf && break
      lprior += Distributions.logpdf(rpriors.removal[i], rparams.removal[i])
    end
  end
  return lprior
end

function logprior(
  sm::NucleicAcidSubstitutionModel,
  sm_priors::Vector{UnivariateDistribution})
  lprior = 0.0
  for i in zip([getproperty(sm, θ) for θ in [propertynames(sm)...]], sm_priors)
    lprior += Distributions.logpdf(i[2], i[1])
    lprior == -Inf && break
  end
  return lprior
end