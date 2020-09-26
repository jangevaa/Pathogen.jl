abstract type IndividualLevelModel end
const ILM = IndividualLevelModel

abstract type TransmissionNetworkILM <: ILM end
const TNILM = TransmissionNetworkILM

abstract type PhylodynamicILM <: ILM end
const PhyloILM = PhylodynamicILM