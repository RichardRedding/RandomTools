struct ARSProposal <: AdaptiveProposal
    dist :: PiecewiseLogLinearDensityProposal
    abscissa :: Array{Float64, 1}
end

draw(d :: ARSProposal) :: Float64 = draw(d.dist)

logmpdf(d :: ARSProposal, x :: Number) :: Float64 = logmpdf(d.dist, x)

function update(d :: ARSProposal, targetlogmpdf :: Function, x :: Number) :: ARSProposal
    if !(x in d.abscissa)
        currentSize = length(d.abscissa)
        idx = (x > last(d.abscissa)) | (currentSize == 0) ? currentSize + 1 : findfirst(a -> x <= a, d.abscissa)
        newAbscissa = idx > currentSize ? push!(copy(d.abscissa), x) : insert!(copy(d.abscissa), idx, x)
        newOrdinate = map(targetlogmpdf, newAbscissa)
        ARSProposalUnsafe(newAbscissa, newOrdinate, lower(d.dist), upper(d.dist))
    else
        d
    end
end

function proposalInputsAreInvalid(abscissa :: Array{T, 1}, ordinate :: Array{S, 1}, lower :: T, upper :: T) :: Bool where {T <: Number, S <: Number}
    #Check that there are at least 3 finite abscissa within the lower and upper bounds that have a finite ordinate
    i = 0
    invalid = true
    iter = zip(Iterators.Stateful(abscissa), Iterators.Stateful(ordinate))
    while invalid && !isempty(iter)
        (abs, ord) = first(iter)
        i = i + (abs <= upper && abs >= lower && isfinite(abs) && isfinite(ord))
        invalid = i < 3
    end
    invalid
end

function ProduceARSIntervalData(interval :: Int, abscissa :: Array{T, 1}, ordinate :: Array{S, 1}, lower :: T, upper :: T) :: NamedTuple{(:lower, :upper, :intercept, :slope), Tuple{Number, Number, Number, Number}} where {T <: Number, S <: Number}
    numAbscissa = length(abscissa)
    numIntervals = numAbscissa < 3 ? 1 + numAbscissa : 2 * (numAbscissa - 1)
    lineIdx = max(1, min(numAbscissa - 1, 1 + ceil(Int64, interval / 2) - (interval % 2) * 2))
    slope = (ordinate[lineIdx + 1] - ordinate[lineIdx]) / (abscissa[lineIdx + 1] - abscissa[lineIdx])
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
        intersectingLineIdx = lineIdx + offsetDirection * 2
        slope0 = (ordinate[intersectingLineIdx + 1] - ordinate[intersectingLineIdx]) / (abscissa[intersectingLineIdx + 1] - abscissa[intersectingLineIdx])
        intercept0 = ordinate[intersectingLineIdx + 1] - abscissa[intersectingLineIdx + 1] * slope0
        bound1 = slope0 != slope ? (intercept - intercept0) / (slope0 - slope) : (abscissa[lineIdx + offsetDirection] + abscissa[lineIdx + 1 + offsetDirection]) / 2
        bound2 = abscissa[ceil(Int64, (interval + 1) / 2)]
        lowerIsAbscissa ? (lower = bound2, upper = bound1, intercept = intercept, slope = slope) : (lower = bound1, upper = bound2, intercept = intercept, slope = slope)
    end
end

function ProduceARSProposalData(abscissa :: Array{T, 1}, ordinate :: Array{S, 1}, lower :: T, upper :: T) :: Array{NamedTuple{(:m, :d), Tuple{Number, BoundedLogLinear}}} where {T <: Number, S <: Number}
    numAbscissa = length(abscissa)
    numIntervals = numAbscissa < 3 ? 1 + numAbscissa : 2 * (numAbscissa - 1)
    proposalData = Array{NamedTuple{(:m, :d), Tuple{Number, BoundedLogLinear}}}(undef, numIntervals)
    w = 0.0
    for i in 1:numIntervals
        (l, u, intercept, slope) = ProduceARSIntervalData(i, abscissa, ordinate, lower, upper)
        w = w + (abs(slope) > 0 ? sign(slope) * (exp(intercept + slope * u - log(abs(slope))) - exp(intercept + slope * l - log(abs(slope)))) : exp(intercept) * (u - l))
        proposalData[i] = (m = w, d = BoundedLogLinear(l, u, slope))
    end
    proposalData
end

function ARSProposalUnsafe(abscissa :: Array{T, 1}, ordinate :: Array{S, 1}, lower :: T, upper :: T) :: ARSProposal where {T <: Number, S <: Number}
    ARSProposal(PiecewiseLogLinearDensityProposal(ProduceARSProposalData(abscissa, ordinate, lower, upper)), abscissa)
end

function ARSProposal(abscissa :: Array{T, 1}, targetlogmpdf :: Function, lower :: T, upper :: T) :: ARSProposal where {T <: Number, S <: Number}
    ordinate = map(targetlogmpdf, abscissa)
    if proposalInputsAreInvalid(abscissa, ordinate, lower, upper)
        error("More than 3 abscissa with finite ordinate must be provided")
    else
        ARSProposalUnsafe(abscissa, ordinate, lower, upper)
    end
end

function ARSProposal(abscissa :: Array{T, 1}, ordinate :: Array{S, 1}, lower :: T, upper :: T) :: ARSProposal where {T <: Number, S <: Number}
    if proposalInputsAreInvalid(abscissa, ordinate, lower, upper)
        error("More than 3 abscissa with finite ordinate must be provided")
    else
        ARSProposalUnsafe(abscissa, ordinate, lower, upper)
    end
end