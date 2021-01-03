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
import Unitful

import ..StatMech: ho_free_energy, ho_internal_energy, ho_entropy, ho_vol_sp_ht
import DimensionalData: name

export Wv,
    Br,
    Temp,
    Vol,
    Press,
    HoFreeEnergy,
    HoInternalEnergy,
    HoEntropy,
    HoVolSpHt,
    collectmodes,
    sample_bz

const Wv = const Wavevector = Dim{:Wavevector}  # TODO: Should I add more constraints?
const Br = const Branch = Dim{:Branch}
const Temp = const Temperature = Dim{:Temperature}
const Vol = const Volume = Dim{:Volume}
const Press = const Pressure = Dim{:Pressure}
const NormalModes = AbstractDimMatrix{T,<:Union{Tuple{Wv,Br},Tuple{Br,Wv}}} where {T}
const TempVolOrVolTemp = Union{Tuple{Temp,Vol},Tuple{Vol,Temp}}
const TempIndepNormalModes = AbstractDimVector{<:NormalModes,<:Tuple{Vol}}
const TempDepNormalModes = AbstractDimMatrix{<:NormalModes,<:TempVolOrVolTemp}

function testconverge(t, ωs, wₖs, N = 3)
    perm = sortperm(wₖs; by = length)
    fe = map(ωs[perm[(end-N+1):end]], wₖs[perm[(end-N+1):end]]) do ω, wₖ
        ho_free_energy(t, ω, wₖ)
    end
    return all(y / x < 1 for (x, y) in zip(fe, fe[2:end]))
end

function collectmodes(ω::AbstractDimArray{T,3})::TempIndepNormalModes where {T}
    @argcheck all(hasdim(ω, (Branch, Wavevector, Vol)))
    return DimArray([ωᵥ for ωᵥ in eachslice(ω; dims = Vol)], dims(ω, Vol))
end
function collectmodes(ω::AbstractDimArray{T,4})::TempDepNormalModes where {T}
    @argcheck all(hasdim(ω, (Branch, Wavevector, Vol, Temp)))
    M, N = map(Base.Fix1(size, ω), (Temp, Vol))
    return DimArray([ω[Temp(i), Vol(j)] for i in 1:M, j in 1:N], dims(ω, (Temp, Vol)))
end

foreach((:ho_free_energy, :ho_internal_energy, :ho_entropy, :ho_vol_sp_ht)) do f
    quote
        # Relax the constraint on wₖ, it can even be a 2×1 matrix!
        function $f(t::Unitful.Temperature, ω::NormalModes, wₖ)  # Scalar
            if any(wₖ .<= zero(eltype(wₖ)))  # Must hold, or else wₖ is already wrong
                throw(DomainError("all weights should be greater than 0!"))
            end
            wₖ = wₖ ./ sum(wₖ)  # Normalize weights
            fₙₖ = map(Base.Fix1($f, t), ω)  # Physical property on each harmonic oscillator
            return _sample_bz(fₙₖ, wₖ)  # Scalar
        end
    end |> eval
end

for (T, f) in zip(
    (:HoFreeEnergy, :HoInternalEnergy, :HoEntropy, :HoVolSpHt),
    (:ho_free_energy, :ho_internal_energy, :ho_entropy, :ho_vol_sp_ht),
)
    expr = quote
        function $T(ω::TempIndepNormalModes, wₖ, ax::TempVolOrVolTemp)
            t, v = dims(ax, (Temp, Vol))
            arr = [sample_bz(x -> $f(t₀, x), ωᵥ, wₖ) for t₀ in t, ωᵥ in ω]  # Slower than `eachslice(ω; dims = Vol)`
            return swapdims(DimArray(arr, (t, v)), map(typeof, ax))
        end
        $T(ω::TempIndepNormalModes, wₖ, t::Union{Temp,Tuple{<:Temp}}) =
            $T(ω, wₖ, (t, dims(ω, Vol)...))
        function $T(ω::TempDepNormalModes, wₖ, ax::TempVolOrVolTemp = dims(ω, (Temp, Vol)))
            t, v = dims(axes, (Temp, Vol))
            M, N = size(ω)
            arr = [sample_bz(x -> $f(t[i], x), ω[i, j], wₖ) for i in 1:M, j in 1:N]  # `eachslice` is not easy to use here
            return swapdims(DimArray(arr, (t, v)), map(typeof, ax))
        end
    end
    eval(expr)
end

_sample_bz(fₙₖ::AbstractDimMatrix{T,<:Tuple{Br,Wv}}, wₖ) where {T} = sum(fₙₖ * collect(wₖ))  # `collect` allows wₖ to be a tuple
_sample_bz(fₖₙ::AbstractDimMatrix{T,<:Tuple{Wv,Br}}, wₖ) where {T} =
    _sample_bz(transpose(fₖₙ), wₖ)  # Just want to align axis, `transpose` is enough.

name(::Type{<:Wv}) = "Wavevector"
name(::Type{<:Br}) = "Branch"
name(::Type{<:Vol}) = "Volume"
name(::Type{<:Temp}) = "Temperature"
name(::Type{<:Press}) = "Pressure"

end
