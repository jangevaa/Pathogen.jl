function _transmission_pathway(pop::Population,
                               i::Int64, j::Int64)
  return pop.risks[[i; j], :x], pop.risks[[i; j], :y]
end

@recipe function f(tn::TransmissionNetwork,
                   pop::Population,
                   events::Events{T},
                   time::Float64;
                   show_states=true::Bool) where T <: EpidemicModel
  xguide --> ""
  yguide --> ""
  legend --> :none
  xlims --> extrema(pop.risks[!, :x]) .+ (sum(extrema(pop.risks[!, :x]).*(-1,1)) .* (-0.05, 0.05))
  ylims --> extrema(pop.risks[!, :y]) .+ (sum(extrema(pop.risks[!, :y]).*(-1,1)) .* (-0.05, 0.05))
  aspect_ratio := :equal
  markerstrokecolour --> :black
  axis --> nothing
  titlefontcolor --> :black
  n = pop.individuals
  susceptible = _ids_by_state(events, State_S, time)
  for i = 1:n
    if i ∉ susceptible
      source = findfirst(tn.internal[:, i])
      if source != nothing
        @series begin
          x, y = _transmission_pathway(pop, i, source)
          seriestype := :path
          seriescolor --> :black
          label := ""
          x, y
        end
      end
    end
  end
  if show_states
    @series begin
      pop, events, time
    end
  end
end

@recipe function f(tn::TransmissionNetwork,
                   pop::Population) where T <: EpidemicModel
  xguide --> ""
  yguide --> ""
  legend --> :none
  xlims --> extrema(pop.risks[!, :x]) .+ (sum(extrema(pop.risks[!, :x]).*(-1,1)) .* (-0.05, 0.05))
  ylims --> extrema(pop.risks[!, :y]) .+ (sum(extrema(pop.risks[!, :y]).*(-1,1)) .* (-0.05, 0.05))
  aspect_ratio := :equal
  axis --> nothing
  titlefontcolor --> :black
  n = pop.individuals
  for i = 1:n
    if !tn.external[i] & any(tn.internal[:, i])
      @series begin
        source = findfirst(tn.internal[:, i])
        x, y = _transmission_pathway(pop, i, source)
        seriestype := :path
        seriescolor --> :black
        label := ""
        x, y
      end
    end
  end
end
