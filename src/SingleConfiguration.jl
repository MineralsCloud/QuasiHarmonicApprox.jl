module SingleConfiguration

using DimensionalData
using EquationsOfState.Collections
using EquationsOfState.NonlinearFitting
using EquationsOfState.Find
using OptionalArgChecks: @argcheck
using Unitful: Temperature, Frequency, Energy, @u_str

import ..StatMech: free_energy

export WaveVector, Branch, Frequencies, Weights, times

@dim WaveVector "WaveVector"
@dim Branch "Branch"

Frequencies = DimensionalArray{
    T,
    2,
    <:Union{Tuple{WaveVector,Branch},Tuple{Branch,WaveVector}},
} where {T}

Weights = DimensionalArray{T,1,<:Tuple{WaveVector}} where {T}

function free_energy(t::Temperature, ω::Frequencies, wₖ::Weights, e0::Energy = 0u"eV")
    @argcheck all(w >= zero(w) for w in wₖ)
    wₖ = wₖ / sum(wₖ)
    f = map(Base.Fix1(free_energy, t), ω)
    sum(times(f, wₖ)) + e0
end

times(ω::DimensionalArray{T,2,<:Tuple{WaveVector,Branch}}, wₖ::Weights) where {T} = ω' * wₖ
times(ω::DimensionalArray{T,2,<:Tuple{Branch,WaveVector}}, wₖ::Weights) where {T} = ω * wₖ
function v_from_p(t0, ω, wk, e0, p, eos)
    f_t0v = free_energy(t0, ω, wk, e0)
    eos = lsqfit(eos(Collections.Energy()), volumes, f_t0v)
    return findvolume(eos(Pressure()), p, (0.5, 1.3) .* eos.v0)
end

function interpolate_f_v(f, t0, ω, wk, e0, p, eos, volumes)
    v = v_from_p(t0, ω, wk, e0, p, eos)
    return interpolate(f, volumes)(v)
end

end
