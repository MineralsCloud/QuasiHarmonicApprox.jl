"""
# module Tools



# Examples

```jldoctest
julia>
```
"""
module Tools

export differentiate,
    deepflatten

function differentiate(x::T, f::T)::T where {T <: AbstractVector}
    length(x) == length(f) ? n = length(x) : throw(DimensionMismatch("The two arguments must have the same length!"))

    derivative = zeros(n)
    derivative[1] = (f[2] - f[1]) / (x[2] - x[1])
    for i in 2:(n - 1)
        derivative[i] = (f[i + 1] - f[i - 1]) / (x[i + 1] - x[i - 1])
    end
    derivative[n] = (f[n] - f[n - 1]) / (x[n] - x[n - 1])
    derivative
end
function differentiate(x::T, f::T, dim::Int)::T where {T <: AbstractMatrix}
    m, n = size(f)
    derivative = zeros(m, n)

    if dim == 1
        for j = 1:n
            derivative[:, j] = differentiate(x[:, j], f[:, j])  # for each column of `x`
        end
    elseif dim == 2
        for i = 1:m
            derivative[i, :] = differentiate(x[i, :], f[i, :])  # for each row of `x`
        end
    else
        throw(DomainError(dim, "The `dim` variable must be `1` or `2`!"))
    end
    derivative
end

function deepflatten(arr::Vector{<: Vector})
    dim = [1]

    function recursiveflatten(arr, dim)
        if arr isa Vector{<: Vector}
            recursiveflatten(collect(Iterators.flatten(arr)),
                pushfirst!(dim, length(arr) / prod(dim)))
        else
            arr, pushfirst!(dim, length(arr) / prod(dim))
        end
    end

    flattened, dim = recursiveflatten(arr, dim)
    reshape(flattened, reverse(dim)[2:end]...)
end

end
