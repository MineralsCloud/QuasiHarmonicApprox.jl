using ConstructionBase: constructorof

export Dimension, BidimensionalData, hasdim, dimnum, dims, isdimequal

abstract type Dimension{T,A<:AbstractVector{T}} <: AbstractVector{T} end
function (::Type{T})(data) where {S,A,T<:Dimension{<:S,<:A}}
    # Can't use `S` since if `T` has no type parameter, `S` would not be defined!
    data = map(Base.Fix1(convert, eltype(T)), data)
    return constructorof(T){eltype(T),typeof(data)}(data)
end
abstract type BidimensionalData{X<:Dimension,Y<:Dimension,T,Z<:AbstractMatrix{T}} <:
              AbstractMatrix{T} end

function hasdim(A::BidimensionalData, dim::Type{<:Dimension{<:T,<:S}}) where {T,S}
    return A.x isa dim || A.y isa dim
end
hasdim(A::BidimensionalData, dim::Dimension) = A.x == dim || A.y == dim

function dimnum(A::BidimensionalData, dim::Integer)
    @assert 0 < dim <= ndims(A) "A doesn't have $dim dimensions!"
    return dim
end
function dimnum(A::BidimensionalData, dim::Type{<:Dimension{<:T,<:S}}) where {T,S}
    if A.x isa dim
        return 1
    elseif A.y isa dim
        return 2
    else
        throw(DimensionMismatch("A doesn't have dimension `$dim`!"))
    end
end
function dimnum(A::BidimensionalData, dim::Dimension)
    if A.x == dim
        return 1
    elseif A.y == dim
        return 2
    else
        throw(DimensionMismatch("A doesn't have dimension `$dim`!"))
    end
end

dims(A::BidimensionalData) = (A.x, A.y)
dims(A::BidimensionalData, dim) = dims(A)[dimnum(A, dim)]

function isdimequal(A::BidimensionalData, Bs::BidimensionalData...)
    return foldl(&, all(hasdim(B, dim) for dim in dims(A)) for B in Bs; init=true)
end

Base.size(A::Dimension) = size(parent(A))
# See https://github.com/rafaqz/DimensionalData.jl/blob/bd28d08/src/array/array.jl#L67
Base.size(A::BidimensionalData) = size(parent(A))
# See https://github.com/rafaqz/DimensionalData.jl/blob/bd28d08/src/array/array.jl#L74
@inline Base.size(A::BidimensionalData, dim) = size(parent(A), dimnum(A, dim))  # Here, `parent(A)` is necessary to avoid `StackOverflowError`.

# See https://github.com/rafaqz/DimensionalData.jl/blob/bd28d08/src/array/array.jl#L68
Base.axes(A::BidimensionalData) = axes(parent(A))
# See https://github.com/rafaqz/DimensionalData.jl/blob/bd28d08/src/array/array.jl#L73
@inline Base.axes(A::BidimensionalData, dim) = axes(parent(A), dimnum(A, dim))

Base.getindex(A::Dimension, i...) = getindex(parent(A), i...)
Base.getindex(A::BidimensionalData, i...) = getindex(parent(A), i...)

function Base.transpose(A::BidimensionalData)
    return constructor(A)(A.y, A.x, transpose(parent(A)))  # `constructorof` is needed!
end

Base.parent(A::Dimension) = A.data
Base.parent(A::BidimensionalData) = A.z

# See https://github.com/rafaqz/DimensionalData.jl/blob/bd28d08/src/array/methods.jl#L95-L99
function Base.map(f, As::Dimension...)
    A = first(As)
    @assert foldl(&, map(==(A), As))
    newdata = map(f, map(parent, As)...)
    return constructor(A)(newdata)  # `constructorof` is needed!
end
function Base.map(f, As::BidimensionalData...)
    A = first(As)
    @assert isdimequal(As...)
    newdata = map(f, map(parent, As)...)
    return constructor(A)(dims(A)..., newdata)  # `constructorof` is needed!
end

function Base.:(==)(A::Dimension, B::Dimension)
    return constructor(A) == constructor(B) && parent(A) == parent(B)
end
# See https://github.com/rafaqz/DimensionalData.jl/blob/bd28d08/src/array/array.jl#L76-L77
function Base.:(==)(A::BidimensionalData, B::BidimensionalData)
    return dims(A) == dims(B) && parent(A) == parent(B)
end

# See https://github.com/rafaqz/DimensionalData.jl/blob/bd28d08/src/array/methods.jl#L101-L105
function Base.mapslices(f, A::BidimensionalData; dims)
    dimnums = collect(dimnum(A, dim) for dim in dims)
    newdata = mapslices(f, parent(A); dims=dimnums)
    return constructor(A)(dims(A)..., newdata)
end

function Base.eachslice(A::BidimensionalData; dims)
    if dims isa Tuple && length(dims) != 1
        throw(ArgumentError("only single dimensions are supported"))
    end
    dim = dimnum(A, dims)
    idx1, idx2 = ntuple(d -> (:), dim - 1), ntuple(d -> (:), ndims(A) - dim)
    return (view(A, idx1..., i, idx2...) for i in axes(A, dim))
end

constructor(x) = constructorof(typeof(x))
