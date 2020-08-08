module SingleConfiguration

using EquationsOfState.Collections
using EquationsOfState.NonlinearFitting
using EquationsOfState.Find
using OptionalArgChecks: @argcheck
using Unitful: Temperature, Frequency, Energy, @u_str

import ..StatMech: free_energy

function free_energy(
    t::Temperature,
    ω::AbstractArray{<:Frequency,3},
    wₖ::AbstractVector,
    e0::AbstractVector{<:Energy} = zeros(size(ω, 3)) * u"eV",
)
    @argcheck all(w >= zero(w) for w in wₖ)
    wₖ = wₖ / sum(wₖ)
    f = map(Base.Fix1(free_energy, t), ω)
    g = dropdims(sum(f; dims = 1); dims = 1)
    return g' * wₖ + e0
end
function free_energy(
    t::AbstractVector{<:Temperature},
    ω::AbstractArray{<:Frequency,4},
    wₖ::AbstractVector,
    e0::AbstractMatrix{<:Energy} = zeros(size(ω)[3:4]) * u"eV",
)
    f = map(t, eachslice(ω; dims = 3), eachrow(e0)) do t1, ω1, e1
        free_energy(t1, ω1, wₖ, e1)
    end
    reshape(Iterators.flatten(f) |> collect, size(ω)[3:4])
end

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
