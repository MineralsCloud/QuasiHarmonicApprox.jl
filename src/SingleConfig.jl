module SingleConfig

using DimensionalData:
    AbstractDimMatrix, AbstractDimArray, DimArray, Dim, dims, val, refdims, swapdims, dimnum
using EquationsOfStateOfSolids.Collections: BirchMurnaghan3rd, EnergyEOS, PressureEOS
using EquationsOfStateOfSolids.Fitting: eosfit
using EquationsOfStateOfSolids.Volume
using Unitful: Temperature, Frequency, Energy, Wavenumber

import DimensionalData
import ..StatMech: ho_free_energy

export Wavevector, Branch, Temp, Vol, v2p

const Wavevector = Dim{:Wavevector}
const Branch = Dim{:Branch}
const Temp = Dim{:Temp}
const Vol = Dim{:Vol}
const Press = Dim{:Press}
const NormalMode = Union{Tuple{Wavevector,Branch},Tuple{Branch,Wavevector}}

function ho_free_energy(
    t,
    ω::AbstractDimMatrix{T,<:NormalMode,<:AbstractMatrix{T}},
    wₖ,
) where {T}
    wₖ = wₖ ./ sum(wₖ)  # Normalize weights
    fₕₒ = ho_free_energy.(t, ω)  # free energy on each harmonic oscillator
    return sum(sample_bz(fₕₒ, wₖ))  # Scalar
end

# Relax the constraint on wₖ, it can even be a 2×1 matrix!
function sample_bz(ω::AbstractDimMatrix{T,<:Tuple{Branch,Wavevector}}, wₖ) where {T}
    if any(wₖ .<= zero(eltype(wₖ)))  # Must hold, or else wₖ is already wrong
        throw(DomainError("All the values of the weights should be greater than 0!"))
    end
    return ω * collect(wₖ)  # Allow wₖ to be a tuple
end
sample_bz(ω::AbstractDimMatrix{T,<:Tuple{Wavevector,Branch}}, wₖ) where {T} =
    sample_bz(transpose(ω), wₖ)  # Just want to align axis, `transpose` is enough.

function v2p(
    fₜᵥ,
    initparam = BirchMurnaghan3rd(
        minimum(dims(fₜᵥ, Vol)),
        zero(eltype(fₜᵥ)) / minimum(dims(fₜᵥ, Vol)),
        4,
    ),
)
    t, v = dims(fₜᵥ, (Temp, Vol))
    volumes = val(v)
    arr = map(eachslice(fₜᵥ; dims = Temp)) do fₜ₀ᵥ
        eosparam = eosfit(EnergyEOS(initparam), volumes, fₜ₀ᵥ)
        p = map(PressureEOS(eosparam), volumes)
        fₜ₀ᵥ = if dimnum(fₜᵥ, Temp) == 1
            DimArray(reshape(fₜ₀ᵥ, 1, :), (Temp([val(refdims(fₜ₀ᵥ))]), v))
        else
            DimArray(reshape(fₜ₀ᵥ, :, 1), (v, Temp([val(refdims(fₜ₀ᵥ))])))
        end
        replacedim(fₜ₀ᵥ, Vol => Press(p))
    end
    return DimArray(arr, (t,))
end

DimensionalData.name(::Type{<:Wavevector}) = "Wavevector"
DimensionalData.name(::Type{<:Branch}) = "Branch"

DimensionalData.shortname(::Type{<:Wavevector}) = "𝐪"
DimensionalData.shortname(::Type{<:Branch}) = "𝑛"

function replacedim(A::AbstractDimArray, dimensions::Pair...)
    for (olddim, newdim) in dimensions
        alldims = Any[dims(A)...]
        alldims[findfirst(x -> x isa olddim, alldims)] = newdim
        A = swapdims(A, Tuple(alldims))
    end
    return A
end

end
