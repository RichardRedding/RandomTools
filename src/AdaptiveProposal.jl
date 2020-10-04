abstract type AdaptiveProposal end

function draw(d :: AdaptiveProposal) :: Number
    error("no method for draw($(typeof(d))) :: Number")
end

function update(d :: AdaptiveProposal, targetlogmpdf :: Function, x :: Number) :: AdaptiveProposal
    error("no method for update($(typeof(d)), Function, Number) :: AdaptiveProposal")
end

function logmpdf(d :: AdaptiveProposal, x :: Number) :: Real
    error("no method for logmpdf($(typeof(d)), Number) :: Real")
end