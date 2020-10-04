module SingleConfig

using DimensionalData:
    AbstractDimArray,
    AbstractDimMatrix,
    AbstractDimVector,
    DimArray,
    Dim,
    dims,
    swapdims,
    hasdim
using OptionalArgChecks: @argcheck
using Unitful: Temperature, Energy

import DimensionalData
import ..StatMech: ho_free_energy, ho_internal_energy, ho_entropy, ho_vol_sp_ht

export Wavevector,
    Branch,
    Temp,
    Vol,
    Press,
    FreeEnergy,
    TempIndepNormalModes,
    TempDepNormalModes,
    sample_bz

const Wavevector = Dim{:Wavevector}  # TODO: Should I add more constraints?
const Branch = Dim{:Branch}
const Temp = Dim{:Temp}
const Vol = Dim{:Vol}
const Press = Dim{:Press}
const NormalModes = AbstractDimMatrix{
    T,
    <:Union{Tuple{Wavevector,Branch},Tuple{Branch,Wavevector}},
} where {T}
const TempVolOrVolTemp = Union{Tuple{Temp,Vol},Tuple{Vol,Temp}}
const TempIndepNormalModes = AbstractDimVector{<:NormalModes,<:Tuple{Vol}}
const TempDepNormalModes = AbstractDimMatrix{<:NormalModes,<:TempVolOrVolTemp}
const TempVolOrVolTempField = AbstractDimMatrix{T,<:TempVolOrVolTemp} where {T}
const FreeEnergy = DimArray{<:Energy,2,<:TempVolOrVolTemp}
const InternalEnergy = DimArray{<:Energy,2,<:TempVolOrVolTemp}
const Entropy = DimArray{T,2,<:TempVolOrVolTemp} where {T}
const VolSpHt = DimArray{T,2,<:TempVolOrVolTemp} where {T}

function testconverge(t, ωs, wₖs, N = 3)
    perm = sortperm(wₖs; by = length)
    fe = map(ωs[perm[(end-N+1):end]], wₖs[perm[(end-N+1):end]]) do ω, wₖ
        ho_free_energy(t, ω, wₖ)
    end
    return all(y / x < 1 for (x, y) in zip(fe, fe[2:end]))
end

function TempIndepNormalModes(ω::AbstractDimArray{T,3})::TempIndepNormalModes where {T}
    @argcheck all(hasdim(ω, (Branch, Wavevector, Vol)))
    return DimArray([ωᵥ for ωᵥ in eachslice(ω; dims = Vol)], dims(ω, Vol))
end

function TempDepNormalModes(ω::AbstractDimArray{T,4})::TempDepNormalModes where {T}
    @argcheck all(hasdim(ω, (Branch, Wavevector, Vol, Temp)))
    M, N = map(Base.Fix1(size, ω), (Temp, Vol))
    return DimArray([ω[Temp(i), Vol(j)] for i in 1:M, j in 1:N], dims(ω, (Temp, Vol)))
end

for (T, f) in zip(
    (:FreeEnergy, :InternalEnergy, :Entropy, :VolSpHt),
    (:ho_free_energy, :ho_internal_energy, :ho_entropy, :ho_vol_sp_ht),
)
    expr = quote
        function $T(ω::TempIndepNormalModes, wₖ, axes)::TempVolOrVolTempField
            t, v = dims(axes, Temp), Base.Fix2(dims, Vol)(hasdim(axes, Vol) ? axes : ω)
            M, N = map(length, (t, v))
            arr = [sample_bz(x -> $f(t[i], x), ω[j], wₖ) for i in 1:M, j in 1:N]  # Slower than `eachslice(ω; dims = Vol)`
            return swapdims(DimArray(reshape(arr, (M, N)), (t, v)), axes)
        end
        function $T(
            ω::TempDepNormalModes,
            wₖ,
            axes = dims(ω, (Temp, Vol)),
        )::TempVolOrVolTempField
            t, v = dims(axes, (Temp, Vol))
            M, N = map(length, (t, v))
            arr = [sample_bz(x -> $f(t[i], x), ω[i, j], wₖ) for i in 1:M, j in 1:N]  # `eachslice` is not easy to use here
            return swapdims(DimArray(arr, (t, v)), axes)
        end
    end
    eval(expr)
end

# Relax the constraint on wₖ, it can even be a 2×1 matrix!
function sample_bz(f, ω::NormalModes, wₖ)  # Scalar
    wₖ = wₖ ./ sum(wₖ)  # Normalize weights
    fₙₖ = f.(ω)  # Physical property on each harmonic oscillator
    return sample_bz(fₙₖ, wₖ)  # Scalar
end
function sample_bz(fₙₖ::AbstractDimMatrix{T,<:Tuple{Branch,Wavevector}}, wₖ) where {T}
    if any(wₖ .<= zero(eltype(wₖ)))  # Must hold, or else wₖ is already wrong
        throw(DomainError("All the values of the weights should be greater than 0!"))
    end
    return sum(fₙₖ * collect(wₖ))  # `collect` allows wₖ to be a tuple
end
sample_bz(fₖₙ::AbstractDimMatrix{T,<:Tuple{Wavevector,Branch}}, wₖ) where {T} =
    sample_bz(transpose(fₖₙ), wₖ)  # Just want to align axis, `transpose` is enough.

DimensionalData.name(::Type{<:Wavevector}) = "Wavevector"
DimensionalData.name(::Type{<:Branch}) = "Branch"

end
