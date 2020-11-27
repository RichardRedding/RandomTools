struct IA2RMS <: AdaptiveSamplingRegime end

function draw(r :: IA2RMS, d :: AdaptiveProposal, targetlogmpdf :: Function, init :: Float64) :: Tuple{Float64, AdaptiveProposal}
    candidate, d2, ldmtargetCandidate, ldmproposalCandidate = arsStep(d, targetlogmpdf)
    ldmtargetInit = targetlogmpdf(init)
    ldmproposalInit = logmpdf(d2, init)
    proposalLogRatio = min(ldmtargetCandidate, ldmproposalCandidate) - min(ldmtargetInit, ldmproposalInit)
    stationaryDensityLogRatio = ldmtargetCandidate - ldmtargetInit
    metropolisAccept = Base.rand() <= exp(stationaryDensityLogRatio - proposalLogRatio)
    if metropolisAccept
        if Base.rand() > exp(ldmproposalInit - ldmtargetInit)
            (candidate, update(d2, targetlogmpdf, init))
        else
            (candidate, d2)
        end
    else
        if Base.rand() > exp(ldmproposalCandidate - ldmtargetCandidate)
            (init, update(d2, targetlogmpdf, candidate))
        else
            (init, d2)
        end
    end
end