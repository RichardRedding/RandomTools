function unfold(A)
    V = []
    for x in A
        if !(typeof(x) <: AbstractArray)
            push!(V, x)
        else
            append!(V, unfold(x))
        end
    end
    V
end