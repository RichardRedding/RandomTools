struct BoundedLogLinear
    lower :: Float64
    upper :: Float64
    slope :: Float64
end

struct PiecewiseLogLinearDensityProposal
    dist :: Array{NamedTuple{(:m, :d), Tuple{Float64, BoundedLogLinear}}}
end

function logpdf(d :: BoundedLogLinear, x :: Number) :: Float64
    abs(d.slope) > 0 ? x * d.slope + log(abs(d.slope)) - log(abs(exp(d.slope * d.upper) - exp(d.slope * d.lower))) : -log(d.upper - d.lower)
end

function quantile(d :: BoundedLogLinear, q :: Number) :: Float64
    abs(d.slope) > 0 ? log(q * (exp(d.slope * d.upper) - exp(d.slope * d.lower)) + exp(d.slope * d.lower)) / d.slope : d.lower + q * (d.upper - d.lower)
end

function cdf(d :: BoundedLogLinear, x :: Number) :: Float64
    abs(d.slope) > 0 ? (exp(d.slope * x) - exp(d.slope * d.lower)) / (exp(d.slope * d.upper) - exp(d.slope * d.lower)) : (x - d.lower) / (d.upper - d.lower)
end

lower(d :: PiecewiseLogLinearDensityProposal) = first(d.dist).d.lower

upper(d :: PiecewiseLogLinearDensityProposal) = last(d.dist).d.upper

draw(d :: PiecewiseLogLinearDensityProposal) :: Number = quantile(d, Base.rand())

findIntervalInDomain(d :: PiecewiseLogLinearDensityProposal, x :: Number) :: Int = findfirst(y -> x <= y.d.upper, d.dist)

function findIntervalInCodomain(d :: PiecewiseLogLinearDensityProposal, p :: Number) :: Int
    findfirst(x -> x.m == Inf || p * last(d.dist).m <= x.m, d.dist)
end

function logmpdf(d :: PiecewiseLogLinearDensityProposal, x :: Number) :: Float64
    if lower(d) <= x && x <= upper(d)
        i = findIntervalInDomain(d, x)
        mi, di = d.dist[i]
        wl = i == 1 ? 0 : d.dist[i-1].m
        logpdf(di, x) + log(mi - wl)
    else
        -Inf
    end
end

logpdf(d :: PiecewiseLogLinearDensityProposal, x :: Number) :: Float64 = logmpdf(d, x) - log(last(d.dist).m)

function quantile(d :: PiecewiseLogLinearDensityProposal, q :: Number) :: Number
    if q <= 0
        lower(d)
    elseif q >= 1
        upper(d)
    else
        i = findIntervalInCodomain(d, q)
        mi, di = d.dist[i]
        wl = i == 1 ? 0 : d.dist[i-1].m
        q2 = (q * last(d.dist).m - wl) / (mi - wl)
        quantile(di, q2)
    end
end

function cdf(d :: PiecewiseLogLinearDensityProposal, x :: Number) :: Float64
    if x < lower(d)
        0.0
    elseif x > upper(d)
        1.0
    else
        i = findIntervalInDomain(d, x)
        mi, di = d.dist[i]
        wl = i == 1 ? 0 : d.dist[i-1].m
        p = cdf(di, x)
        (wl + p * (mi - wl)) / last(d.dist).m
    end
end
