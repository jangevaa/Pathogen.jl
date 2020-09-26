function _transmission_pathway(pop::Population,
                               i::Int64, j::Int64)
  return pop.risks[[i; j], :x], pop.risks[[i; j], :y]
end

@recipe function f(tn::TransmissionNetwork,
                   pop::Population,
                   events::Events{T},
                   time::Float64;
                   show_individuals=true,
                   show_states=true) where T <: DiseaseStateSequence
  if show_states & !show_individuals
    @warn "Will not show disease states without also showing individuals"
  end
  xguide --> ""
  yguide --> ""
  legend --> :none
  xlims --> extrema(pop.risks[!, :x]) .+ (sum(extrema(pop.risks[!, :x]).*(-1,1)) .* (-0.05, 0.05))
  ylims --> extrema(pop.risks[!, :y]) .+ (sum(extrema(pop.risks[!, :y]).*(-1,1)) .* (-0.05, 0.05))
  aspect_ratio := :equal
  linecolor --> :black
  axis --> nothing
  titlefontcolor --> :black
  framestyle --> :box
  n = individuals(pop)
  susceptible = _ids_by_state(events, State_S, time)
  for i = 1:n
    if i ∉ susceptible
      source = findfirst(tn.internal[:, i])
      if source != nothing
        @series begin
          x, y = _transmission_pathway(pop, i, source)
          seriestype := :path
          label := ""
          x, y
        end
      end
    end
  end
  if show_individuals
    if show_states
      @series begin
        pop, events, time
      end
    else
      @series begin
        pop
      end
    end
  end
end

@recipe function f(
  sim::Simulation{T, M},
  time::Float64;
  show_individuals=true,
  show_states=true) where {
  T <: DiseaseStateSequence,
  M <: ILM}
  show_individuals := show_individuals
  show_states := show_states
  sim.transmission_network, sim.population, sim.events, time
end

@recipe function f(tn::TransmissionNetwork,
                   pop::Population;
                   show_individuals=true)
  xguide --> ""
  yguide --> ""
  legend --> :none
  xlims --> extrema(pop.risks[!, :x]) .+ (sum(extrema(pop.risks[!, :x]).*(-1,1)) .* (-0.05, 0.05))
  ylims --> extrema(pop.risks[!, :y]) .+ (sum(extrema(pop.risks[!, :y]).*(-1,1)) .* (-0.05, 0.05))
  aspect_ratio := :equal
  axis --> nothing
  titlefontcolor --> :black
  linecolor --> :black
  framestyle --> :box
  n = individuals(pop)
  for i = 1:n
    if !tn.external[i] & any(tn.internal[:, i])
      @series begin
        source = findfirst(tn.internal[:, i])
        x, y = _transmission_pathway(pop, i, source)
        seriestype := :path
        label := ""
        x, y
      end
    end
  end
  if show_individuals
    @series begin
      pop
    end
  end
end
