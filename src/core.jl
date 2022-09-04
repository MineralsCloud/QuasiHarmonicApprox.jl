using ConstructionBase: constructorof

export Dimension, BidimensionalData

abstract type Dimension{T,A<:AbstractVector{T}} <: AbstractVector{T} end
abstract type BidimensionalData{X<:Dimension,Y<:Dimension,T,Z<:AbstractMatrix{T}} <:
              AbstractMatrix{T} end

hasdim(A::BidimensionalData, dim::Type{<:Dimension}) = A.x isa dim || A.y isa dim
hasdim(A::BidimensionalData, dim::Dimension) = A.x == dim || A.y == dim

function dimnum(A::BidimensionalData, dim::Type{<:Dimension})
    @assert hasdim(A, dim)
    return A.x isa dim ? 1 : 2
end
function dimnum(A::BidimensionalData, dim::Dimension)
    @assert hasdim(A, dim)
    return A.x == dim ? 1 : 2
end

Base.size(A::Dimension) = size(A.data)
Base.size(A::BidimensionalData) = size(A.z)
Base.size(A::BidimensionalData, dim::Type{<:Dimension}) = size(A.z, dimnum(A, dim))

Base.IndexStyle(::Type{<:Dimension}) = IndexLinear()
Base.IndexStyle(::Type{<:BidimensionalData}) = IndexLinear()

Base.getindex(A::Dimension, i...) = getindex(A.data, i...)
Base.getindex(A::BidimensionalData, i...) = getindex(A.z, i...)

function Base.transpose(A::BidimensionalData)
    return constructorof(typeof(A))(A.y, A.x, transpose(A.z))
end
