module RandomTools

export ARS, ARMS, IA2RMS, ARSProposal, ARMSProposal, draw, drawn, findAbscissa

include("Utilities.jl")

#Adaptive proposals
include("AdaptiveProposal.jl")
include("PiecewiseLogLinearDensityProposal.jl")
include("ARSProposal.jl")
include("ARMSProposal.jl")

#Adaptive sampling regimes
include("AdaptiveSamplingRegime.jl")
include("ARS.jl")
include("ARMS.jl")
include("IA2RMS.jl")

end
