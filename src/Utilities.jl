function findAbscissa(targetlogmpdf :: Function, lower :: Float64, upper :: Float64) :: Array{Float64, 1}
    shrinkRate = 0.05
    flower, fupper = targetlogmpdf(lower), targetlogmpdf(upper)
    while !(isfinite(flower) && isfinite(fupper))
        d = (upper - lower) * shrinkRate
        if flower > fupper
            upper -= d * 0.5
            fupper = targetlogmpdf(upper)
        else
            lower += d * 0.5
            flower = targetlogmpdf(lower)
        end
    end
    
    middle = (lower + upper) * 0.5
    while !((flower < targetlogmpdf(middle) > fupper) || (abs((middle - lower) / (upper - lower) - 0.5) > 0.45))
        d = flower > fupper ? (lower - middle) : (upper - middle)
        middle += d * 0.5
    end

    [lower, middle, upper]
end