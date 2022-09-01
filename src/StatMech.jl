module StatMech

using Functors: fmap, @functor
using Unitful: Frequency, Energy, Wavenumber, NoUnits, ħ, k, c0

export HarmonicOscillator
export bose_einstein_dist,
    partition_function, free_energy, internal_energy, entropy, volumetric_heat_capacity

struct HarmonicOscillator{T<:Frequency}
    ω::T
    function HarmonicOscillator(ω)
        @assert ω >= zero(ω)
        return new(ω)
    end
end
@functor HarmonicOscillator

bose_einstein_dist(ho::HarmonicOscillator, t) = 1 / expm1(ħ * ho.ω / (k * t))

function partition_function(ho::HarmonicOscillator, t)
    return iszero(ho.ω) ? 1 : csch(ħ * ho.ω / (2k * t)) / 2
end

function free_energy(ho::HarmonicOscillator, t)
    if iszero(ho.ω)
        return zero(k * t)
    else
        ħω, kt = ħ * ho.ω, k * t
        return ħω / 2 + kt * log(-expm1(-ħω / kt))
    end
end

function internal_energy(ho::HarmonicOscillator, t)
    if iszero(ho.ω)
        return k * t
    else
        ħω = ħ * ho.ω / 2
        return ħω * coth(ħω / (k * t))
    end
end

function entropy(ho::HarmonicOscillator, t)
    if iszero(t) || iszero(ho.ω)
        return zero(k)
    else
        n = bose_einstein_dist(t, ho.ω)
        return k * ((1 + n) * log1p(n) - n * log(n))
    end
end

function volumetric_heat_capacity(ho::HarmonicOscillator, t)
    if iszero(t)
        return zero(k)
    else
        if iszero(ho.ω)
            return k
        else
            x = NoUnits(ħ * ho.ω / (k * 2t))
            return k * (x * csch(x))^2
        end
    end
end

foreach((
    :bose_einstein_dist,
    :partition_function,
    :free_energy,
    :internal_energy,
    :entropy,
    :volumetric_heat_capacity,
)) do func
    # See https://docs.julialang.org/en/v1/manual/metaprogramming/#Code-Generation
    @eval $func(ho::HarmonicOscillator, t) = $func(fmap(to_ω, ho.ω), t)
end

end
