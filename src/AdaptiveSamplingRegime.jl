abstract type AdaptiveSamplingRegime end

function draw(r :: AdaptiveSamplingRegime, d :: AdaptiveProposal, targetlogmpdf :: Function) :: Tuple{Number, AdaptiveProposal}
    #draw where proposal sampling is not based on markov chain
    error("no method for draw($(typeof(r)), $(typeof(d)), Function) :: Tuple{Number, AdaptiveProposal}")
end

function draw(r :: AdaptiveSamplingRegime, d :: AdaptiveProposal, targetlogmpdf :: Function, init :: Number) :: Tuple{Number, AdaptiveProposal}
    #draw where proposal sampling is based on markov chain
    draw(r, d, targetlogmpdf)
end

function drawn(r :: AdaptiveSamplingRegime, d :: AdaptiveProposal, targetlogmpdf :: Function, n :: Int) :: Tuple{Array{Number, 1}, AdaptiveProposal}
    #drawn where proposal sampling is not based on markov chain
    result = Array{Number, 1}(undef, n)
    for i in 1:n
        y, d1 = draw(r, d, targetlogmpdf)
        result[i] = y
        d = d1
    end
    (result, d)
end

function drawn(r :: AdaptiveSamplingRegime, d :: AdaptiveProposal, targetlogmpdf :: Function, init :: Number, n :: Int) :: Tuple{Array{Number, 1}, AdaptiveProposal}
    #drawn where proposal sampling is based on markov chain
    result = Array{Number, 1}(undef, n)
    for i in 1:n
        y, d1 = draw(r, d, targetlogmpdf, init)
        init = y
        result[i] = y
        d = d1
    end
    (result, d)
end
