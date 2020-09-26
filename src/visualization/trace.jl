function Base.length(x::NASM)
  return length(fieldnames(typeof(x)))
end

function Base.getindex(x::NASM, i::Int64)
  return getproperty(x, fieldnames(typeof(x))[i])
end

@recipe function f(iter,
                   x::Vector{RiskParameters{T}}) where {T <: DiseaseStateSequence}
  xguide --> "Iteration"
  yguide --> "Value"
  legend --> :none
  seriestype := :line
  # label --> permutedims(_parameters(x[1]))
  iter, convert(Array{Float64, 2}, x[iter])
end

@recipe function f(iter, x::Vector{T}) where {T <: NASM}
  xguide --> "Iteration"
  yguide --> "Value"
  legend --> :none
  seriestype := :line
  # label --> permutedims(_parameters(x[1]))
  iter, convert(Array{Float64, 2}, x[iter])
end

@recipe function f(iter, p::Int64, x::MCMC{S, M}) where {S <: DiseaseStateSequence, M <: ILM}
  np = length(x.markov_chains[1].risk_parameters[1])
  sp = 0
  if M == PhyloILM
    sp += length(x.markov_chains[1].substitution_model[1])
  end
  xguide --> "Iteration"
  yguide --> "Value"
  guidefontsize --> 8
  seriestype --> :line
  if 0 <= p <= np
    for i in eachindex(x.markov_chains)
      @series begin
        label --> "Markov chain $i"
        iter, getindex.(x.markov_chains[i].risk_parameters[iter], p)
      end
    end
  elseif np < p <= (np + sp)
    for i in eachindex(x.markov_chains)
      @series begin
        label --> "Markov chain $i"
        iter, getindex.(x.markov_chains[i].substitution_model[iter], p-np)
      end
    end
  else
    error("Invalid parameter index")
  end
end

@recipe function f(iter, x::MCMC{S, M}) where {S <: DiseaseStateSequence, M <: ILM}
  np = length(x.markov_chains[1].risk_parameters[1])
  if M == PhyloILM
    np += length(x.markov_chains[1].substitution_model[1])
  end
  layout --> (np, 1)
  legend --> :none
  link --> :y
  for i = 1:np
    @series begin
      subplot := i
      iter, i, x
    end
  end
end