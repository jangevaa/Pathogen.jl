primitive type DiseaseState 8 end

DiseaseState(x::UInt8) = reinterpret(DiseaseState, x)
UInt8(x::DiseaseState) = reinterpret(UInt8, x)

# `DieaseState` const representations
const State_S = DiseaseState(0b0001)
const State_E = DiseaseState(0b0010)
const State_I = DiseaseState(0b0100)
const State_R = DiseaseState(0b1000)

const DiseaseStates = Vector{DiseaseState}

const _DiseaseStateVector = [State_S, State_E, State_I, State_R]
const _DiseaseStateCharVector = ['S', 'E', 'I', 'R']

# `DiseaseState` conversion to/from `Char`
function Base.convert(::Type{Char}, x::DiseaseState)
  return _DiseaseStateCharVector[_DiseaseStateVector .== Ref(x)][1]
end

function Base.show(io::IO, x::DiseaseState)
  return print(io, convert(Char, x))
end

Base.convert(::Type{DiseaseStates}, x::Type{SEIR}) = [State_S, State_E, State_I, State_R]
Base.convert(::Type{DiseaseStates}, x::Type{SEI}) = [State_S, State_E, State_I]
Base.convert(::Type{DiseaseStates}, x::Type{SIR}) = [State_S, State_I, State_R]
Base.convert(::Type{DiseaseStates}, x::Type{SI}) = [State_S, State_I]

Base.convert(::Type{Tuple}, x::Type{SEIR}) = (State_S, State_E, State_I, State_R)
Base.convert(::Type{Tuple}, x::Type{SEI}) = (State_S, State_E, State_I)
Base.convert(::Type{Tuple}, x::Type{SIR}) = (State_S, State_I, State_R)
Base.convert(::Type{Tuple}, x::Type{SI}) = (State_S, State_I)

Base.iterate(::Type{S}, v...) where {S <: DiseaseStateSequence} = iterate(convert(Tuple, S), v...)