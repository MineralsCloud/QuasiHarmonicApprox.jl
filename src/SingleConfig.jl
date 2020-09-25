module SingleConfig

using DimensionalData:
    AbstractDimArray, AbstractDimMatrix, AbstractDimVector, DimArray, Dim, dims
using Unitful: Temperature, Energy

import DimensionalData
import ..StatMech: ho_free_energy, ho_internal_energy, ho_entropy, ho_vol_sp_ht

export Wavevector, Branch, Temp, Vol, Press

const Wavevector = Dim{:Wavevector}
const Branch = Dim{:Branch}
const Temp = Dim{:Temp}
const Vol = Dim{:Vol}
const Press = Dim{:Press}
const NormalModes = AbstractDimMatrix{
    T,
    <:Union{Tuple{Wavevector,Branch},Tuple{Branch,Wavevector}},
} where {T}
const TempIndependentNormalModes = AbstractDimArray{
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
} where {T}
const TempDependentNormalModes = AbstractDimArray{
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
} where {T}

function testconverge(t, ωs, wₖs, N = 3)
    perm = sortperm(wₖs; by = length)
    fe = map(ωs[perm[(end-N+1):end]], wₖs[perm[(end-N+1):end]]) do ω, wₖ
        ho_free_energy(t, ω, wₖ)
    end
    return all(y / x < 1 for (x, y) in zip(fe, fe[2:end]))
end

for f in (:ho_free_energy, :ho_internal_energy, :ho_entropy, :ho_vol_sp_ht)
    eval(
        quote
            function $f(t::Temperature, ω::NormalModes, wₖ)
                wₖ = wₖ ./ sum(wₖ)  # Normalize weights
                fₙₖ = $f.(t, ω)  # Physical property on each harmonic oscillator
                return sample_bz(fₙₖ, wₖ)  # Scalar
            end
            function $f(t, ω::TempIndependentNormalModes, wₖ)
                return DimArray(
                    [$f(t₀, ωᵥ, wₖ) for t₀ in t, ωᵥ in eachslice(ω; dims = Vol)],
                    (Temp(t), dims(ω, Vol)),
                )
            end
            function $f(t, ω::TempDependentNormalModes, wₖ)
                adims = dims(ω, (Temp, Vol))
                arr = map(eachslice(ω; dims = Temp), t) do ωₜ₀ᵥ, t₀
                    map(eachslice(ωₜ₀ᵥ; dims = Vol)) do ωₜᵥ
                        $f(t₀, ωₜᵥ, wₖ)
                    end
                end
                arr = reshape(collect(Iterators.flatten(arr))', length.(adims))
                return DimArray(arr, adims)
            end
        end,
    )
end


# Relax the constraint on wₖ, it can even be a 2×1 matrix!
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
