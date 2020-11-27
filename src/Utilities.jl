function findAbscissa(targetlogmpdf :: Function, lower :: Float64, upper :: Float64) :: Array{Float64, 1}
    shrinkRate = 0.05
    flower, fupper = targetlogmpdf(lower), targetlogmpdf(upper)
    i = 0
    while !(isfinite(flower) && isfinite(fupper)) && i < 10
        d = (upper - lower) * shrinkRate
        lower += d * 0.5
        upper -= d * 0.5
        flower, fupper = targetlogmpdf(lower), targetlogmpdf(upper)
        i+=1
    end
    
    i = 0
    middle = (lower + upper) * 0.5
    while !((flower < targetlogmpdf(middle) > fupper) || (0.95 < (middle - lower) / (upper - lower) < 0.05)) && i < 10
        println((middle - lower) / (upper - lower), " ", flower < targetlogmpdf(middle) > fupper)
        d = flower > fupper ? (lower - middle) : (upper - middle)
        middle += d * 0.5
        i+=1
    end

    [lower, middle, upper]
end