module SingleConfiguration

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

end
