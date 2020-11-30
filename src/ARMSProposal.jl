struct ARMSProposal <: AdaptiveProposal
    dist :: PiecewiseLogLinearDensityProposal
    abscissa :: Array{Float64, 1}
end

draw(d :: ARMSProposal) :: Float64 = draw(d.dist)

logmpdf(d :: ARMSProposal, x :: Float64) :: Float64 = logmpdf(d.dist, x)

function update(d :: ARMSProposal, targetlogmpdf :: Function, x :: Float64) :: ARMSProposal
    i = findfirst(a -> x <= a, d.abscissa)
    newAbscissa = x > last(d.abscissa) ? push!(copy(d.abscissa), x) : insert!(copy(d.abscissa), i, x)
    newOrdinate = map(targetlogmpdf, newAbscissa)
    ARMSProposal(newAbscissa, newOrdinate, lower(d.dist), upper(d.dist))
end

isConvex(slopePrevious :: Float64, slope :: Float64, slopeNext :: Float64) :: Bool = slopeNext <= slope && slope <= slopePrevious

function calcGradient(x1 :: Float64, x2 :: Float64, y1 :: Float64, y2 :: Float64) :: Float64
    numerator = y2 - y1
    denominator = x2 - x1
    denominator == 0 ? 0.0 : numerator / denominator
end

calcIntercept(x :: Float64, y :: Float64, g :: Float64) :: Float64 = y - x * g

calcIntersection(intercept0 :: Float64, intercept1 :: Float64, slope0 :: Float64, slope1 :: Float64) :: Float64 = (intercept1 - intercept0) / (slope0 - slope1)

function ProduceARMSIntervalData(interval :: Int, abscissa :: Array{Float64, 1}, ordinate :: Array{Float64, 1}, lower :: Float64, upper :: Float64) :: NamedTuple{(:lower, :upper, :intercept, :slope), Tuple{Float64, Float64, Float64, Float64}}
    numAbscissa = length(abscissa)
    numIntervals = numAbscissa < 3 ? 1 + numAbscissa : 2 * (numAbscissa - 1)
    lineIdxIfConcave = ceil(Int64, interval / 2)
    previousSlope = lineIdxIfConcave == 1 ? Inf : calcGradient(abscissa[lineIdxIfConcave - 1], abscissa[lineIdxIfConcave], ordinate[lineIdxIfConcave - 1], ordinate[lineIdxIfConcave])
    thisSlope = calcGradient(abscissa[lineIdxIfConcave], abscissa[lineIdxIfConcave + 1], ordinate[lineIdxIfConcave], ordinate[lineIdxIfConcave + 1])
    nextSlope = lineIdxIfConcave == (numAbscissa - 1) ? -Inf : calcGradient(abscissa[lineIdxIfConcave + 1], abscissa[lineIdxIfConcave + 2], ordinate[lineIdxIfConcave + 1], ordinate[lineIdxIfConcave + 2])
    intervalIsConvex = isConvex(previousSlope, thisSlope, nextSlope)
    lineIdxIfConvex = max(1, min(numAbscissa - 1, 1 + ceil(Int64, interval / 2) - (interval % 2) * 2))
    lineIdx = intervalIsConvex ? lineIdxIfConvex : lineIdxIfConcave
    slope = intervalIsConvex ? calcGradient(abscissa[lineIdx], abscissa[lineIdx + 1], ordinate[lineIdx], ordinate[lineIdx + 1]) : thisSlope
    intercept = calcIntercept(abscissa[lineIdx + 1], ordinate[lineIdx + 1], slope)
    if interval == 1
        (lower = lower, upper = abscissa[1], intercept = intercept, slope = slope)
    elseif interval == 2
        (lower = abscissa[1], upper = abscissa[2], intercept = intercept, slope = slope)
    elseif interval == (numIntervals - 1)
        (lower = abscissa[numAbscissa - 1], upper = last(abscissa), intercept = intercept, slope = slope)
    elseif interval == numIntervals
        (lower = last(abscissa), upper = upper, intercept = intercept, slope = slope)
    else
        lowerIsAbscissa = isodd(interval)
        offsetDirection = (lowerIsAbscissa * 2 - 1)
        bound1 = if intervalIsConvex
            intersectingLineIdx = lineIdx + offsetDirection * 2
            slope0 = calcGradient(abscissa[intersectingLineIdx], abscissa[intersectingLineIdx + 1], ordinate[intersectingLineIdx], ordinate[intersectingLineIdx + 1])
            intercept0 = calcIntercept(abscissa[intersectingLineIdx + 1], ordinate[intersectingLineIdx + 1], slope0)
            slope0 != slope ? calcIntersection(intercept0, intercept, slope0, slope) : (abscissa[lineIdx + offsetDirection] + abscissa[1 + lineIdx + offsetDirection]) / 2
        else
            (abscissa[lineIdx + 1] + abscissa[lineIdx]) / 2
        end
        bound2 = abscissa[ceil(Int64, (interval + 1) / 2)]
        lowerIsAbscissa ? (lower = bound2, upper = bound1, intercept = intercept, slope = slope) : (lower = bound1, upper = bound2, intercept = intercept, slope = slope)
    end
end

function ProduceARMSProposalData(abscissa :: Array{Float64, 1}, ordinate :: Array{Float64, 1}, lower :: Float64, upper :: Float64) :: Array{NamedTuple{(:m, :d), Tuple{Float64, BoundedLogLinear}}}
    numAbscissa = length(abscissa)
    numIntervals = numAbscissa < 3 ? 1 + numAbscissa : 2 * (numAbscissa - 1)
    proposalData = Array{NamedTuple{(:m, :d), Tuple{Float64, BoundedLogLinear}}}(undef, numIntervals)
    w = 0.0
    for i in 1:numIntervals
        (l, u, intercept, slope) = ProduceARMSIntervalData(i, abscissa, ordinate, lower, upper)
        d = BoundedLogLinear(l, u, slope)
        w += if slope == 0
            exp(intercept) * (u - l)
        else
            frac = slope > 0 ? 1 - exp(-slope * (u - l)) : 1 - exp(slope * (u - l))
            k = slope > 0 ? exp(intercept + slope * u - log(abs(slope))) : exp(intercept + slope * l - log(abs(slope)))
            frac * k
        end
        proposalData[i] = (m = w, d = d)
    end
    proposalData
end

function ARMSProposalUnsafe(abscissa :: Array{Float64, 1}, ordinate :: Array{Float64, 1}, lower :: Float64, upper :: Float64) :: ARMSProposal
    ARMSProposal(PiecewiseLogLinearDensityProposal(ProduceARMSProposalData(abscissa, ordinate, lower, upper)), abscissa)
end

function ARMSProposal(abscissa :: Array{Float64, 1}, targetlogmpdf :: Function, lower :: Float64, upper :: Float64) :: ARMSProposal
    ordinate = map(targetlogmpdf, abscissa)
    if proposalInputsAreInvalid(abscissa, ordinate, lower, upper)
        error("More than 3 abscissa with finite ordinate must be provided")
    else
        ARMSProposalUnsafe(abscissa, ordinate, lower, upper)
    end
end

function ARMSProposal(abscissa :: Array{Float64, 1}, ordinate :: Array{Float64, 1}, lower :: Float64, upper :: Float64) :: ARMSProposal
    if proposalInputsAreInvalid(abscissa, ordinate, lower, upper)
        error("More than 3 abscissa with finite ordinate must be provided")
    else
        ARMSProposalUnsafe(abscissa, ordinate, lower, upper)
    end
end