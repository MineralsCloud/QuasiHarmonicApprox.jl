using ConstructionBase: constructorof

abstract type Dimension{T} <: AbstractVector{T} end
abstract type BidimensionalData{X<:Dimension,Y<:Dimension,Z} <: AbstractMatrix{Z} end
(func::Type{<:BidimensionalData})(x::X, y::Y, z::Z) where {X,Y,Z} = func{X,Y,Z}(x, y, z)

hasdim(A::BidimensionalData, dim::Type{<:Dimension}) = A.x isa dim || A.y isa dim

function dimnum(A::BidimensionalData, dim::Type{<:Dimension})
    @assert hasdim(A, dim)
    return A.x isa dim ? 1 : 2
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
