function _epidemic_curve_distribution(
  x::Vector{Events{S}},
  state::DiseaseState,
  times;
  credibleinterval::Float64=0.95) where {
    S <: DiseaseStateSequence}
  !(0.0 <= credibleinterval <= 1.0) && error("`credibleinterval` must be between 0 and 1")
  ci = (1-credibleinterval)/2
  
  lwr  = Array{Tuple{Float64,Float64}, 1}(undef, length(times))
  upr  = Array{Tuple{Float64,Float64}, 1}(undef, length(times))
  for t in eachindex(times)
    y = Array{Int64, 1}(undef, length(x))
    @simd for i in eachindex(x)
      y[i] = _count_by_state(x[i], state, times[t])
    end
    q = quantile(y, [ci, 1-ci])
    lwr[t] = (times[t], q[1])
    upr[t] = (times[t], q[2])
  end
  return [lwr; reverse(upr)]
end

@recipe function f(
  x::Vector{Events{S}},
  state::DiseaseState,
  times;
  credibleinterval=0.95,
  thin=1,
  burnin=0)  where {
    S <: DiseaseStateSequence}
  xguide --> "Time"
  yguide --> "N"
  label --> convert(Char, state)
  legend --> :topright
  fillcolor --> _state_color(state)
  linecolor --> _state_color(state)
  fillalpha --> 0.5
  linealpha --> 0.7
  seriestype := :shape
  _epidemic_curve_distribution(x[1+burnin:thin:end], state, times, credibleinterval=credibleinterval)
end

@recipe function f(
  x::Vector{Events{S}},
  times;
  credibleinterval=0.95,
  thin=1,
  burnin=0) where {
    S <: DiseaseStateSequence}
  for s in convert(Tuple, S)
    @series begin
      credibleinterval := credibleinterval
      burnin := burnin
      thin := thin
      x, s, times
    end
  end
end

@recipe function f(
  x::MCMC,
  state::DiseaseState,
  times;
  credibleinterval=0.95,
  thin=1,
  burnin=0,
  bychain=false)
  if !bychain
    y = x.markov_chains[1].events[1+burnin:thin:end]
    for i in 2:length(x.markov_chains)
      append!(y, x.markov_chains[i].events[1+burnin:thin:end])
    end
    @series begin
      xguide --> "Time"
      yguide --> "N"
      legend --> :topright
      label --> convert(Char, state)
      fillcolor --> _state_color(state)
      linecolor --> _state_color(state)
      fillalpha --> 0.5
      linealpha --> 0.7
      seriestype := :shape
      _epidemic_curve_distribution(y, state, times, credibleinterval=credibleinterval)
    end
  else
    label --> convert(Char, state)
    for i in eachindex(x.markov_chains)
      @series begin
        burnin := burnin
        thin := thin
        credibleinterval := credibleinterval
        label := :none
        x.markov_chains[i].events, state, times
      end
    end
  end
end


@recipe function f(
  x::MCMC{S, M},
  times;
  credibleinterval=0.95,
  thin=1,
  burnin=0,
  bychain=false) where {
    S <: DiseaseStateSequence,
    M <: ILM}
  for s in convert(Tuple, S)
    @series begin
      credibleinterval := credibleinterval
      burnin := burnin
      thin := thin
      bychain := bychain
      x, s, times
    end
  end
end