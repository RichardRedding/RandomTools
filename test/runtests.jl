println("Running tests")

include("Utilities.jl")

using Test, Random, Distributions, RandomTools

convexTargets = [
[(1.0,Normal(0,1))]
]

nonConvexTargets = [
[(0.3, Normal(-5,1)), (0.3, Normal(1,1)), (0.4, Normal(7,1))]
]

testpdf(x :: Array{Tuple{Float64, Normal{Float64}}, 1}, y :: Number) = foldl((l, (w, d)) -> l + exp(log(w) + logpdf(d, y)), x, init = 0)
testlogpdf(x :: Array{Tuple{Float64, Normal{Float64}}, 1}, y :: Number) = log(testpdf(x, y))
testcdf(x :: Array{Tuple{Float64, Normal{Float64}}, 1}, y :: Number) = foldl((l, (w, d)) -> l + w * cdf(d, y), x, init = 0)

convexSamplers = [
(target, n) -> drawn(ARS(), ARSProposal(findAbscissa(target, -20.0, 20.0), target, -20.0, 20.0), target, 0, n),
(target, n) -> drawn(ARMS(), ARMSProposal(findAbscissa(target, -20.0, 20.0), target, -20.0, 20.0), target, 0, n),
(target, n) -> drawn(IA2RMS(), ARMSProposal(findAbscissa(target, -20.0, 20.0), target, -20.0, 20.0), target, 0, n),
]

nonConvexSamplers = [
(target, n) -> drawn(ARMS(), ARMSProposal(findAbscissa(target, -20.0, 20.0), target, -20.0, 20.0), target, 0, n),
(target, n) -> drawn(IA2RMS(), ARMSProposal(findAbscissa(target, -20.0, 20.0), target, -20.0, 20.0), target, 0, n)
]

convexTest = map(convexTargets) do target
    map(convexSamplers) do sampler
        ps = [0.005, 0.1, 0.25, 0.5, 0.75, 0.9, 0.995]
        n = Int64(1e4)
        x, prop = sampler(x -> testlogpdf(target, x), n)
        map(ps) do p
            estimate = p
            exact = testcdf(target, quantile(x, p))
            (abs(estimate - exact) / sqrt(exact * (1 - exact) / n)) > 2
            #"estimate = $(p), exact = $(exact), target = $(repr(target)), sampler = $(repr(typeof(prop)))"
        end
    end
end |> unfold

nonConvexTest = map(nonConvexTargets) do target
    map(nonConvexSamplers) do sampler
        ps = [0.005, 0.1, 0.25, 0.5, 0.75, 0.9, 0.995]
        n = Int64(1e4)
        x, prop = sampler(x -> testlogpdf(target, x), n)
        map(ps) do p
            estimate = p
            exact = testcdf(target, quantile(x, p))
            (abs(estimate - exact) / sqrt(exact * (1 - exact) / n)) > 2
            #"estimate = $(p), exact = $(exact), target = $(repr(target)), sampler = $(repr(typeof(prop)))"
        end
    end
end |> unfold

@test sum(convexTest) == 0
@test sum(nonConvexTest) == 0