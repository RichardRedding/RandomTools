struct ARMS <: AdaptiveSamplingRegime end

function draw(r :: ARMS, d :: AdaptiveProposal, targetlogmpdf :: Function, init :: Number) :: Tuple{Number, AdaptiveProposal}
    candidate, d2, ldmtargetCandidate, ldmproposalCandidate = arsStep(d, targetlogmpdf)
    ldmtargetInit = targetlogmpdf(init)
    ldmproposalInit = logmpdf(d2, init)
    proposalLogRatio = min(ldmtargetCandidate, ldmproposalCandidate) - min(ldmtargetInit, ldmproposalInit)
    stationaryDensityLogRatio = ldmtargetCandidate - ldmtargetInit
    metropolisAccept = Base.rand() <= exp(stationaryDensityLogRatio - proposalLogRatio)
    result = metropolisAccept ? candidate : init
    (result, d2)
end