struct ARMSProposal <: AdaptiveProposal
    dist :: PiecewiseLogLinearDensityProposal
    abscissa :: Array{Float64, 1}
end

draw(d :: ARMSProposal) :: Float64 = draw(d.dist)

logmpdf(d :: ARMSProposal, x :: Number) :: Float64 = logmpdf(d.dist, x)

function update(d :: ARMSProposal, targetlogmpdf :: Function, x :: Number) :: ARMSProposal
    i = findfirst(a -> x <= a, d.abscissa)
    newAbscissa = x > last(d.abscissa) ? push!(copy(d.abscissa), x) : insert!(copy(d.abscissa), i, x)
    newOrdinate = map(targetlogmpdf, newAbscissa)
    ARMSProposal(newAbscissa, newOrdinate, lower(d.dist), upper(d.dist))
end

isConvex(slopePrevious :: Number, slope :: Number, slopeNext :: Number) :: Bool = slopeNext <= slope && slope <= slopePrevious

function ProduceARMSIntervalData(interval :: Int, abscissa :: Array{T, 1}, ordinate :: Array{S, 1}, lower :: T, upper :: T) :: NamedTuple{(:lower, :upper, :intercept, :slope), Tuple{Number, Number, Number, Number}} where {T <: Number, S <: Number}
    numAbscissa = length(abscissa)
    numIntervals = numAbscissa < 3 ? 1 + numAbscissa : 2 * (numAbscissa - 1)
    lineIdxIfConcave = ceil(Int64, interval / 2)
    previousSlope = lineIdxIfConcave == 1 ? Inf : (ordinate[lineIdxIfConcave] - ordinate[lineIdxIfConcave - 1]) / (abscissa[lineIdxIfConcave] - abscissa[lineIdxIfConcave - 1])
    thisSlope = (ordinate[lineIdxIfConcave + 1] - ordinate[lineIdxIfConcave]) / (abscissa[lineIdxIfConcave + 1] - abscissa[lineIdxIfConcave])
    nextSlope = lineIdxIfConcave == (numAbscissa - 1) ? -Inf : (ordinate[lineIdxIfConcave + 2] - ordinate[lineIdxIfConcave + 1]) / (abscissa[lineIdxIfConcave + 2] - abscissa[lineIdxIfConcave + 1])
    intervalIsConvex = isConvex(previousSlope, thisSlope, nextSlope)
    lineIdxIfConvex = max(1, min(numAbscissa - 1, 1 + ceil(Int64, interval / 2) - (interval % 2) * 2))
    lineIdx = intervalIsConvex ? lineIdxIfConvex : lineIdxIfConcave
    slope = intervalIsConvex ? (ordinate[lineIdx + 1] - ordinate[lineIdx]) / (abscissa[lineIdx + 1] - abscissa[lineIdx]) : thisSlope
    intercept = ordinate[lineIdx + 1] - abscissa[lineIdx + 1] * slope
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
            slope0 = (ordinate[intersectingLineIdx + 1] - ordinate[intersectingLineIdx]) / (abscissa[intersectingLineIdx + 1] - abscissa[intersectingLineIdx])
            intercept0 = ordinate[intersectingLineIdx + 1] - abscissa[intersectingLineIdx + 1] * slope0
            slope0 != slope ? (intercept - intercept0) / (slope0 - slope) : (abscissa[lineIdx + offsetDirection] + abscissa[1 + lineIdx + offsetDirection]) / 2
        else
            (abscissa[lineIdx + 1] + abscissa[lineIdx]) / 2
        end
        bound2 = abscissa[ceil(Int64, (interval + 1) / 2)]
        lowerIsAbscissa ? (lower = bound2, upper = bound1, intercept = intercept, slope = slope) : (lower = bound1, upper = bound2, intercept = intercept, slope = slope)
    end
end

function ProduceARMSProposalData(abscissa :: Array{T, 1}, ordinate :: Array{S, 1}, lower :: T, upper :: T) :: Array{NamedTuple{(:m, :d), Tuple{Number, BoundedLogLinear}}} where {T <: Number, S <: Number}
    numAbscissa = length(abscissa)
    numIntervals = numAbscissa < 3 ? 1 + numAbscissa : 2 * (numAbscissa - 1)
    proposalData = Array{NamedTuple{(:m, :d), Tuple{Number, BoundedLogLinear}}}(undef, numIntervals)
    w = 0.0
    for i in 1:numIntervals
        (l, u, intercept, slope) = ProduceARMSIntervalData(i, abscissa, ordinate, lower, upper)
        w = w + (abs(slope) > 0 ? sign(slope) * (exp(intercept + slope * u - log(abs(slope))) - exp(intercept + slope * l - log(abs(slope)))) : exp(intercept) * (u - l))
        proposalData[i] = (m = w, d = BoundedLogLinear(l, u, slope))
    end
    proposalData
end

function ARMSProposalUnsafe(abscissa :: Array{T, 1}, ordinate :: Array{S, 1}, lower :: T, upper :: T) :: ARMSProposal where {T <: Number, S <: Number}
    ARMSProposal(PiecewiseLogLinearDensityProposal(ProduceARMSProposalData(abscissa, ordinate, lower, upper)), abscissa)
end

function ARMSProposal(abscissa :: Array{T, 1}, targetlogmpdf :: Function, lower :: T, upper :: T) :: ARMSProposal where {T <: Number, S <: Number}
    ordinate = map(targetlogmpdf, abscissa)
    if proposalInputsAreInvalid(abscissa, ordinate, lower, upper)
        error("More than 3 abscissa with finite ordinate must be provided")
    else
        ARMSProposalUnsafe(abscissa, ordinate, lower, upper)
    end
end

function ARMSProposal(abscissa :: Array{T, 1}, ordinate :: Array{S, 1}, lower :: T, upper :: T) :: ARMSProposal where {T <: Number, S <: Number}
    if proposalInputsAreInvalid(abscissa, ordinate, lower, upper)
        error("More than 3 abscissa with finite ordinate must be provided")
    else
        ARMSProposalUnsafe(abscissa, ordinate, lower, upper)
    end
end