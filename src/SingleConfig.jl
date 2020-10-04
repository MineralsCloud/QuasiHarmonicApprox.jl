module SingleConfig

using DimensionalData:
    AbstractDimArray, AbstractDimMatrix, AbstractDimVector, DimArray, Dim, dims, swapdims
using Unitful: Temperature, Energy

import DimensionalData
import ..StatMech: ho_free_energy, ho_internal_energy, ho_entropy, ho_vol_sp_ht

export Wavevector, Branch, Temp, Vol, Press, FreeEnergy, sample_bz

const Wavevector = Dim{:Wavevector}  # TODO: Should I add more constraints?
const Branch = Dim{:Branch}
const Temp = Dim{:Temp}
const Vol = Dim{:Vol}
const Press = Dim{:Press}
const NormalModes = AbstractDimMatrix{
    T,
    <:Union{Tuple{Wavevector,Branch},Tuple{Branch,Wavevector}},
} where {T}
const TempIndependentNormalModes = Union{
    AbstractDimArray{
        T,
        3,
        <:Union{
            Tuple{Wavevector,Branch,Vol},
            Tuple{Branch,Wavevector,Vol},
            Tuple{Vol,Branch,Wavevector},
            Tuple{Vol,Wavevector,Branch},
            Tuple{Wavevector,Vol,Branch},
            Tuple{Branch,Vol,Wavevector},
        },
    },
    AbstractDimVector{<:NormalModes,<:Tuple{Vol}},
} where {T}
const TempDependentNormalModes = Union{
    AbstractDimArray{
        T,
        4,
        <:Union{
            Tuple{Wavevector,Branch,Vol,Temp},
            Tuple{Wavevector,Branch,Temp,Vol},
            Tuple{Wavevector,Vol,Branch,Temp},
            Tuple{Wavevector,Vol,Temp,Branch},
            Tuple{Wavevector,Temp,Branch,Vol},
            Tuple{Wavevector,Temp,Vol,Branch},
            Tuple{Branch,Wavevector,Vol,Temp},
            Tuple{Branch,Wavevector,Temp,Vol},
            Tuple{Branch,Vol,Wavevector,Temp},
            Tuple{Branch,Vol,Temp,Wavevector},
            Tuple{Branch,Temp,Wavevector,Vol},
            Tuple{Branch,Temp,Vol,Wavevector},
            Tuple{Vol,Wavevector,Branch,Temp},
            Tuple{Vol,Wavevector,Temp,Branch},
            Tuple{Vol,Branch,Wavevector,Temp},
            Tuple{Vol,Branch,Temp,Wavevector},
            Tuple{Vol,Temp,Wavevector,Branch},
            Tuple{Vol,Temp,Branch,Wavevector},
            Tuple{Temp,Wavevector,Branch,Vol},
            Tuple{Temp,Wavevector,Vol,Branch},
            Tuple{Temp,Branch,Wavevector,Vol},
            Tuple{Temp,Branch,Vol,Wavevector},
            Tuple{Temp,Vol,Wavevector,Branch},
            Tuple{Temp,Vol,Branch,Wavevector},
        },
    },
    AbstractDimMatrix{<:NormalModes,<:Union{Tuple{Temp,Vol},Tuple{Vol,Temp}}},
} where {T}

function testconverge(t, ωs, wₖs, N = 3)
    perm = sortperm(wₖs; by = length)
    fe = map(ωs[perm[(end-N+1):end]], wₖs[perm[(end-N+1):end]]) do ω, wₖ
        ho_free_energy(t, ω, wₖ)
    end
    return all(y / x < 1 for (x, y) in zip(fe, fe[2:end]))
end

const FreeEnergy = DimArray{<:Energy,2,<:Union{Tuple{Temp,Vol},Tuple{Vol,Temp}}}
const InternalEnergy = DimArray{<:Energy,2,<:Union{Tuple{Temp,Vol},Tuple{Vol,Temp}}}
const Entropy = DimArray{T,2,<:Union{Tuple{Temp,Vol},Tuple{Vol,Temp}}} where {T}
const VolSpHt = DimArray{T,2,<:Union{Tuple{Temp,Vol},Tuple{Vol,Temp}}} where {T}

for (T, f) in zip(
    (:FreeEnergy, :InternalEnergy, :Entropy, :VolSpHt),
    (:ho_free_energy, :ho_internal_energy, :ho_entropy, :ho_vol_sp_ht),
)
    expr = quote
        function $T(ω::TempIndependentNormalModes, wₖ, axes)
            t, v = dims(axes, (Temp, Vol))
            M, N = map(length, (t, v))
            arr = [sample_bz(x -> $f(t[i], x), ω[Vol(j)], wₖ) for i in 1:M, j in 1:N]  # Slower than `eachslice(ω; dims = Vol)`
            return swapdims(DimArray(reshape(arr, (M, N)), (t, v)), axes)
        end
        function $T(ω::TempDependentNormalModes, wₖ, axes)  # FIXME: fix `ω[Vol(i)]` method error
            t, v = dims(axes, (Temp, Vol))
            M, N = map(length, (t, v))
            arr = [
                sample_bz(x -> $f(t[i], x), ω[Temp(i), Vol(j)], wₖ) for i in 1:M, j in 1:N
            ]  # `eachslice` is not easy to use here
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
