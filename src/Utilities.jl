function randomSearch(success :: Function, lower :: Number, upper :: Number)
    stop = false
    result = 0.0
    while !stop
        result = lower + (upper - lower) * Base.rand()
        stop = success(result)
    end
    result
end

function findAbscissa(targetlogmpdf :: Function, lower :: T, upper :: T) :: Array{T, 1} where {T <: Number}
    sort([randomSearch(x -> isfinite(targetlogmpdf(x)), lower, upper), randomSearch(x -> isfinite(targetlogmpdf(x)), lower, upper), randomSearch(x -> isfinite(targetlogmpdf(x)), lower, upper)])
end