module SingleConfig

using DimensionalData: AbstractDimMatrix, AbstractDimArray, DimArray, @dim
using EquationsOfStateOfSolids.Collections
using EquationsOfStateOfSolids.Fitting
using EquationsOfStateOfSolids.Volume
using OptionalArgChecks: @argcheck
using Unitful: Temperature, Frequency, Energy, Wavenumber, @u_str

import ..StatMech: free_energy

export WaveVector, Branch, Frequencies, Weights

@dim WaveVector "WaveVector"
@dim Branch "Branch"

const FrequenciesAxes = Union{Tuple{WaveVector,Branch},Tuple{Branch,WaveVector}}

const Frequencies = DimArray{<:Union{Frequency,Energy,Wavenumber},2,<:FrequenciesAxes}
Frequencies(data, axes::FrequenciesAxes) = DimArray(data, axes)


function free_energy(t::Temperature, ω::Freq, wₖ, e0::Energy = 0u"eV")
    @argcheck all(w >= zero(w) for w in wₖ)
    wₖ = wₖ / sum(wₖ)
    f = map(Base.Fix1(free_energy, t), ω)
    sum(times(f, wₖ)) + e0
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
function sample_bz(ω::AbstractDimMatrix{T,<:Tuple{Branch,WaveVector}}, wₖ) where {T}
    @argcheck all(wₖ .> zero(eltype(wₖ)))
    return ω * wₖ
end
sample_bz(ω::AbstractDimMatrix{T,<:Tuple{WaveVector,Branch}}, wₖ) where {T} =
    sample_bz(ω', wₖ)


function v_from_p(t0, ω, wk, e0, p, eosparam)
    f_t0v = free_energy(t0, ω, wk, e0)
    eos = eosfit(EnergyEOS(eosparam), volumes, f_t0v)
    return mustfindvolume(PressureEOS(eosparam), p)
end

function interpolate_f_v(f, t0, ω, wk, e0, p, eos, volumes)
    v = v_from_p(t0, ω, wk, e0, p, eos)
    return interpolate(f, volumes)(v)
end

end
