__precompile__()

module Pathogen

  # Dependencies
  using DataFrames, Distributions, RecipesBase, ProgressMeter

  import Base.length, Base.convert, Base.show, Base.copy
  import RecipesBase.plot
  import ProgressMeter.next!

  # Types
  include("types/EpidemicModel.jl")
  include("types/DiseaseState.jl")
  include("types/Events/Event.jl")
  include("types/Events/Events.jl")
  include("types/Events/EventExtents.jl")
  include("types/Events/EventObservations.jl")
  include("types/Events/EventRates.jl")
  include("types/Risks/RiskFunctions.jl")
  include("types/Risks/RiskParameters.jl")
  include("types/Risks/RiskPriors.jl")
  include("types/Transmissions/Transmission.jl")
  include("types/Transmissions/TransmissionNetwork.jl")
  include("types/Transmissions/TransmissionRates.jl")
  include("types/Simulation.jl")
  include("types/MarkovChain.jl")
  include("types/MCMC.jl")

  # Functions
  include("functions/_pathway_from.jl")
  include("functions/_pathway_to.jl")
  include("functions/_accept.jl")
  include("functions/_bounds.jl")
  include("functions/generate.jl")
  include("functions/initialize.jl")
  include("functions/observe.jl")
  include("functions/update!.jl")
  include("functions/next!.jl")
  include("functions/simulate!.jl")
  include("functions/loglikelihood.jl")
  include("functions/logpriors.jl")
  include("functions/start!.jl")
  include("functions/iterate!.jl")

  # Visualization
  include("visualization/_count_by_state.jl")
  include("visualization/_ids_by_state.jl")
  include("visualization/_epidemic_curve.jl")
  include("visualization/_population_plot.jl")
  include("visualization/_pathway_to.jl")
  include("visualization/plot.jl")

  export
    SEIR, SEI, SIR, SI,
    DiseaseState, DiseaseStates,
    State_S, State_E, State_I, State_R,
    RiskFunctions, RiskParameters, RiskPriors,
    Simulation,
    next!, simulate!,
    EventObservations, EventExtents,
    observe,
    MCMC, start!, iterate!,
    plot
end
