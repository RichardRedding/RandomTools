using Base: rand

struct ARS <: AdaptiveSamplingRegime end

function arsStep(d :: AdaptiveProposal, targetlogmpdf :: Function) :: Tuple{Number, AdaptiveProposal, Float64, Float64}
    accept = false
    result = 0
    ldmtargetNew = 0
    ldmproposalNew = 0
    while(!accept)
        result = draw(d)
        ldmtargetNew = targetlogmpdf(result)
        ldmproposalNew = logmpdf(d, result)
        acceptanceRatio = min(1, exp(ldmtargetNew - ldmproposalNew))
        accept = Base.rand() <= acceptanceRatio
        d = accept ? d : update(d, targetlogmpdf, result)
    end
    (result, d, ldmtargetNew, ldmproposalNew)
end

draw(r :: ARS, d :: AdaptiveProposal, targetlogmpdf :: Function) :: Tuple{Number, AdaptiveProposal} = arsStep(d, targetlogmpdf)[1:2]