module Pathogen

  # Dependencies
  using Distributed,
        DataFrames,
        Distributions,
        RecipesBase,
        Logging,
        StatsBase,
        Statistics,
        ProgressMeter,
        LinearAlgebra,
        OnlineStats,
        PhyloModels

  # Methods for functions not in Base
  import StatsBase.mean,
         StatsBase.mode,
         ProgressMeter.next!

  # Core
  include("core/IndividualLevelModel.jl")
  include("core/DiseaseStateSequence.jl")
  include("core/DiseaseState.jl")
  include("core/Population.jl")
  include("core/Events/Event.jl")
  include("core/Events/Events.jl")
  include("core/Events/EventRates.jl")
  include("core/Events/EventObservations.jl")
  include("core/Risks/AbstractRisk.jl")
  include("core/Risks/RiskFunctions.jl")
  include("core/Risks/RiskParameters.jl")
  include("core/Transmissions/Transmission.jl")
  include("core/Transmissions/AbstractTransmissionNetwork.jl")
  include("core/Transmissions/TransmissionNetwork.jl")
  include("core/Transmissions/TransmissionRates.jl")
  include("core/Transmissions/TransmissionRateCache.jl")
  include("core/initialize.jl")
  include("core/update!.jl")

  # Simulation
  include("simulation/Simulation.jl")
  include("simulation/generate/Event.jl")
  include("simulation/generate/Tree.jl")
  include("simulation/generate/EventObservations.jl")
  include("simulation/generate/Transmission.jl")
  include("simulation/simulate!.jl")
  include("simulation/update!.jl")

  # Inference
  include("inference/TransmissionNetworkDistribution.jl")
  include("inference/RiskPriors.jl")
  include("inference/EventExtents.jl")
  include("inference/MarkovChain.jl")
  include("inference/MCMC.jl")
  include("inference/_accept.jl")
  include("inference/_bounds.jl")
  include("inference/generate.jl")
  include("inference/initialize.jl")
  include("inference/update!.jl")
  include("inference/loglikelihood.jl")
  include("inference/logprior.jl")
  include("inference/start!.jl")
  include("inference/iterate!.jl")
  include("inference/summary.jl")

  # Visualization
  include("visualization/epidemic_curve.jl")
  include("visualization/epidemic_curve_distribution.jl")
  include("visualization/observations.jl")
  include("visualization/population.jl")
  include("visualization/transmission_network.jl")
  include("visualization/transmission_network_distribution.jl")
  include("visualization/degree_distribution.jl")
  include("visualization/trace.jl")

  # re-export PhyloModels, DataFrames, Distributions
  for pkg in [PhyloModels, DataFrames, Distributions]
    for name in names(pkg)
      @eval export $(name)
    end
  end

  export
    ILM, TransmissionNetworkILM, TNILM, PhylodynamicILM, PhyloILM,
    SEIR, SEI, SIR, SI, DiseaseStateSequence,
    DiseaseState, DiseaseStates,
    State_S, State_E, State_I, State_R,
    RiskFunctions, RiskParameters, RiskPriors,
    Population,
    AbstractTransmissionNetwork,
    TransmissionNetwork,
    TransmissionNetworkDistribution, TNDistribution,
    TransmissionNetworkPrior, TNPrior,
    TransmissionNetworkPosterior, TNPosterior,
    EndogenousTransmission, ExogenousTransmission, NoTransmission,
    Simulation,
    next!, simulate!, generate,
    Events, EventObservations, EventExtents,
    observe,
    MCMC, start!, iterate!,
    summary, mean, mode, innerjoin, leftantijoin
end
