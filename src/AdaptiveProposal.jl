abstract type AdaptiveProposal end

function draw(d :: AdaptiveProposal) :: Float64
    error("no method for draw($(typeof(d)))")
end

function update(d :: AdaptiveProposal, targetlogmpdf :: Function, x :: Float64) :: AdaptiveProposal
    error("no method for update($(typeof(d)), Function, Float64) :: AdaptiveProposal")
end

function logmpdf(d :: AdaptiveProposal, x :: Float64) :: Float64
    error("no method for logmpdf($(typeof(d)), Float64) :: Float64")
end