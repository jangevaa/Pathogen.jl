function _transmission_pathway_density(pop::Population,
                                       tn::TNDistribution,
                                       i::I, j::I) where {I <: Integer}
  x = pop.risks[[i; j], :x]
  y = pop.risks[[i; j], :y]
  dens = tn.internal[j, i]
  @debug "Transmission pathway density for individuals $i and $j" xy1 = (x[1], y[1]) xy2 = (x[2], y[2]) density = dens
  return x, y, dens
end

@recipe function f(tn::TNDistribution,
                   pop::Population;
                   show_individuals=true::Bool) where T <: DiseaseStateSequence
  xguide --> ""
  yguide --> ""
  legend --> :none
  xlims --> extrema(pop.risks[!, :x]) .+ (sum(extrema(pop.risks[!, :x]).*(-1,1)) .* (-0.05, 0.05))
  ylims --> extrema(pop.risks[!, :y]) .+ (sum(extrema(pop.risks[!, :y]).*(-1,1)) .* (-0.05, 0.05))
  aspect_ratio := :equal
  axis --> nothing
  framestyle --> :box
  linecolor --> :black
  n = individuals(pop)
  for i = 1:n, j = 1:n
    if tn.internal[j, i] > 0.0
      @series begin
        x, y, density = _transmission_pathway_density(pop, tn, i, j)
        seriestype --> :path
        linealpha --> density
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
