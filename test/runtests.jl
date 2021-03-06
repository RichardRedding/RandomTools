println("Running tests")

include("Utilities.jl")

using Test, Random, Distributions, RandomTools

convexTargets = [
[(1.0,Normal(0,1))]
]

nonConvexTargets = [
[(0.3, Normal(-5,1)), (0.3, Normal(1,1)), (0.4, Normal(7,1))]
]

testpdf(x :: Array{Tuple{Float64, Normal{Float64}}, 1}, y :: Number) = foldl((l, (w, d)) -> l + exp(log(w) + Distributions.logpdf(d, y)), x, init = 0)
testlogpdf(x :: Array{Tuple{Float64, Normal{Float64}}, 1}, y :: Number) = log(testpdf(x, y))
testcdf(x :: Array{Tuple{Float64, Normal{Float64}}, 1}, y :: Number) = foldl((l, (w, d)) -> l + w * Distributions.cdf(d, y), x, init = 0)

### Test bounded domain
boundedConvexSamplers = [
(target, n) -> drawn(ARS(), ARSProposal(findAbscissa(target, -20.0, 20.0), target, -20.0, 20.0), target, 0.0, n),
(target, n) -> drawn(ARMS(), ARMSProposal(findAbscissa(target, -20.0, 20.0), target, -20.0, 20.0), target, 0.0, n),
(target, n) -> drawn(IA2RMS(), ARMSProposal(findAbscissa(target, -20.0, 20.0), target, -20.0, 20.0), target, 0.0, n),
]

boundedNonConvexSamplers = [
(target, n) -> drawn(ARMS(), ARMSProposal(findAbscissa(target, -20.0, 20.0), target, -20.0, 20.0), target, 0.0, n),
(target, n) -> drawn(IA2RMS(), ARMSProposal(findAbscissa(target, -20.0, 20.0), target, -20.0, 20.0), target, 0.0, n)
]

boundedConvexTest = map(convexTargets) do target
    map(boundedConvexSamplers) do sampler
        ps = [0.005, 0.1, 0.25, 0.5, 0.75, 0.9, 0.995]
        n = Int64(1e5)
        x, prop = sampler(x -> testlogpdf(target, x), n)
        map(ps) do p
            estimate = p
            exact = testcdf(target, Distributions.quantile(x, p))
            (abs(estimate - exact) / sqrt(exact * (1 - exact) / n)) > 2
        end
    end
end |> unfold

boundedNonConvexTest = map(nonConvexTargets) do target
    map(boundedNonConvexSamplers) do sampler
        ps = [0.005, 0.1, 0.25, 0.5, 0.75, 0.9, 0.995]
        n = Int64(1e5)
        x, prop = sampler(x -> testlogpdf(target, x), n)
        map(ps) do p
            estimate = p
            exact = testcdf(target, Distributions.quantile(x, p))
            (abs(estimate - exact) / sqrt(exact * (1 - exact) / n)) > 2
        end
    end
end |> unfold

@test sum(boundedConvexTest) == 0
@test sum(boundedNonConvexTest) == 0

### Test unbounded domain
unboundedConvexSamplers = [
(target, n) -> drawn(ARS(), ARSProposal(findAbscissa(target, -20.0, 20.0), target, -Inf, Inf), target, 0.0, n),
(target, n) -> drawn(ARMS(), ARMSProposal(findAbscissa(target, -20.0, 20.0), target, -Inf, Inf), target, 0.0, n),
(target, n) -> drawn(IA2RMS(), ARMSProposal(findAbscissa(target, -20.0, 20.0), target, -Inf, Inf), target, 0.0, n),
]

unboundedNonConvexSamplers = [
(target, n) -> drawn(ARMS(), ARMSProposal(findAbscissa(target, -20.0, 20.0), target, -Inf, Inf), target, 0.0, n),
(target, n) -> drawn(IA2RMS(), ARMSProposal(findAbscissa(target, -20.0, 20.0), target, -Inf, Inf), target, 0.0, n)
]

unboundedConvexTest = map(convexTargets) do target
    map(unboundedConvexSamplers) do sampler
        ps = [0.005, 0.1, 0.25, 0.5, 0.75, 0.9, 0.995]
        n = Int64(1e5)
        x, prop = sampler(x -> testlogpdf(target, x), n)
        map(ps) do p
            estimate = p
            exact = testcdf(target, Distributions.quantile(x, p))
            (abs(estimate - exact) / sqrt(exact * (1 - exact) / n)) > 2
        end
    end
end |> unfold

unboundedNonConvexTest = map(nonConvexTargets) do target
    map(unboundedNonConvexSamplers) do sampler
        ps = [0.005, 0.1, 0.25, 0.5, 0.75, 0.9, 0.995]
        n = Int64(1e5)
        x, prop = sampler(x -> testlogpdf(target, x), n)
        map(ps) do p
            estimate = p
            exact = testcdf(target, Distributions.quantile(x, p))
            (abs(estimate - exact) / sqrt(exact * (1 - exact) / n)) > 2
        end
    end
end |> unfold

@test sum(unboundedConvexTest) == 0
@test sum(unboundedNonConvexTest) == 0