module SingleConfig

using DimensionalData: AbstractDimMatrix, AbstractDimArray, DimArray, Dim, dims, val
using EquationsOfStateOfSolids.Collections
using EquationsOfStateOfSolids.Fitting
using EquationsOfStateOfSolids.Volume
using OptionalArgChecks: @argcheck
using Unitful: Temperature, Frequency, Energy, Wavenumber

import DimensionalData
import ..StatMech: free_energy

export Wavevector, Branch, Temp, Vol, v2p

const Wavevector = Dim{:Wavevector}
const Branch = Dim{:Branch}
const Temp = Dim{:Temp}
const Vol = Dim{:Vol}
const Press = Dim{:Press}
const FreqAxes = Union{Tuple{Wavevector,Branch},Tuple{Branch,Wavevector}}
const Freq = AbstractDimMatrix{<:Union{Frequency,Energy,Wavenumber},<:FreqAxes}

function free_energy(t::Temperature, ω::Freq, wₖ)
    wₖ /= sum(wₖ)  # Normalize weights
    fₕₒ = free_energy.(t, ω)  # free energy on each harmonic oscillator
    return sum(sample_bz(fₕₒ, wₖ))  # Scalar
end


function free_energy(
    t::AbstractVector{<:Temperature},
    ω::AbstractVector{<:Frequencies},
    e0::AbstractVector{<:Energy} = zeros(length(ω)) * 0u"eV",
)  # For T-independent frequencies
    length(ω) == length(e0) ||
        throw(DimensionMismatch("ω and e0 should be the same length!"))
    return [free_energy(tt, ww, wₖ, e00) for tt in t, (ww, e00) in zip(ω, e0)]
# Relax the constraint on wₖ, it can even be a 2×1 matrix!
function sample_bz(ω::AbstractDimMatrix{T,<:Tuple{Branch,Wavevector}}, wₖ) where {T}
    @argcheck all(wₖ .> zero(eltype(wₖ)))
    return ω * wₖ
end
sample_bz(ω::AbstractDimMatrix{T,<:Tuple{Wavevector,Branch}}, wₖ) where {T} =
    sample_bz(ω', wₖ)

function v2p(fₜᵥ, eosparam)
    t, v = dims(fₜᵥ, (Temp, Vol))
    volumes = val(v)
    map(eachslice(fₜᵥ; dims = Temp), t) do fₜ₀ᵥ, t0
        eos = eosfit(EnergyEOS(eosparam), volumes, fₜ₀ᵥ)
        DimArray(fₜ₀ᵥ, (Temp(t0), Press(map(PressureEOS(eos), volumes))))
    end
end

function interpolate_f_v(f, t0, ω, wk, e0, p, eos, volumes)
    v = v_from_p(t0, ω, wk, e0, p, eos)
    return interpolate(f, volumes)(v)
end

DimensionalData.name(::Type{<:Wavevector}) = "Wavevector"
DimensionalData.name(::Type{<:Branch}) = "Branch"

DimensionalData.shortname(::Type{<:Wavevector}) = "𝐪"
DimensionalData.shortname(::Type{<:Branch}) = "𝑛"

end
