abstract type AbstractTransmissionNetwork end

const AbstractTN = AbstractTransmissionNetwork

function Base.sum(x::TN) where {TN <: AbstractTN}
  return sum(x.external) + sum(x.internal)
end

function Base.show(io::IO, x::TN) where {TN <: AbstractTN}
  return print(io, "$TN Σexternal = $(sum(x.external)) Σinternal = $(sum(x.internal))" )
end

function Base.copy(x::TN) where {TN <: AbstractTN}
  return TN(copy(x.external),
            copy(x.internal))
end

function individuals(x::TN) where {TN <: AbstractTN}
  return length(x.external)
end

function Base.sum(x::TN, i::Integer) where {TN <: AbstractTN}
  if i < 1 || i > individuals(x)
    @error "Invalid individual identifier" i
  end
  return x.external[i] + sum(x.internal[:, i])
end