abstract type AdaptiveSamplingRegime end

function draw(r :: AdaptiveSamplingRegime, d :: AdaptiveProposal, targetlogmpdf :: Function) :: Tuple{Float64, AdaptiveProposal}
    #draw where proposal sampling is not based on markov chain
    error("no method for draw($(typeof(r)), $(typeof(d)), Function) :: Tuple{Float64, AdaptiveProposal}")
end

function draw(r :: AdaptiveSamplingRegime, d :: AdaptiveProposal, targetlogmpdf :: Function, init :: Float64) :: Tuple{Float64, AdaptiveProposal}
    #draw where proposal sampling is based on markov chain
    draw(r, d, targetlogmpdf)
end

function drawn(r :: AdaptiveSamplingRegime, d :: AdaptiveProposal, targetlogmpdf :: Function, n :: Int) :: Tuple{Array{Float64, 1}, AdaptiveProposal}
    #drawn where proposal sampling is not based on markov chain
    result = Array{Float64, 1}(undef, n)
    for i in 1:n
        y, d1 = draw(r, d, targetlogmpdf)
        result[i] = y
        d = d1
    end
    (result, d)
end

function drawn(r :: AdaptiveSamplingRegime, d :: AdaptiveProposal, targetlogmpdf :: Function, init :: Float64, n :: Int) :: Tuple{Array{Float64, 1}, AdaptiveProposal}
    #drawn where proposal sampling is based on markov chain
    result = Array{Float64, 1}(undef, n)
    for i in 1:n
        y, d1 = draw(r, d, targetlogmpdf, init)
        init = y
        result[i] = y
        d = d1
    end
    (result, d)
end
