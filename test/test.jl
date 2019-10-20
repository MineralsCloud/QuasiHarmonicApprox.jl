#=
test.jl:
- Julia version:
- Author: qz
- Date: 2019-01-26
=#
using QuasiHarmonicApproximation.StatisticalMechanics
using QuasiHarmonicApproximation.CoreDataStructures
using QuasiHarmonicApproximation.Sampling
using QuasiHarmonicApproximation.Tools
using QuasiHarmonicApproximation.Interpolation
using QuasiHarmonicApproximation.Thermodynamics
using EquationsOfState
using Dierckx
# using CSV
# using DataFrames
using Setfield: set

temperature = collect(1000.0:1000:10000)
volume = collect(100.0:100:500)
qaxes = (rand(10), rand(6))
w = QSpaceField{:q,:s}(qaxes..., rand(10, 6))
freq = rand(10, 6) * 100 .+ 400

function generate(f)
    ret = zeros(10, 5)
    qs = QSpaceField{:q,:s}(qaxes..., zeros(10, 6))
    for (i, (t, v)) in enumerate(Iterators.product(temperature, volume))
        qs = set(qs, Abstractions.DATALENS, f.(t, freq))
        ret[i] = sample_brillouin_zone(w, qs)
    end
    ret
end

free_energy = ThermodynamicField{:T,:V}(
    temperature,
    volume,
    generate(subsystem_free_energy),
)
internal_energy = ThermodynamicField{:T,:V}(
    temperature,
    volume,
    generate(subsystem_internal_energy),
)
entropy = ThermodynamicField{:T,:V}(temperature, volume, generate(subsystem_entropy))
cv = ThermodynamicField{:T,:V}(
    temperature,
    volume,
    generate(subsystem_volumetric_specific_heat),
)

fits = []
ps = []
for row in eachrow(free_energy)
    fit = fit_energy(BirchMurnaghan3rd(volume[1], 200.0, 4.0), volume, row)
    push!(fits, fit)
    push!(ps, map(eval_pressure(BirchMurnaghan3rd(fit.param[1:end-1]...)), volume))
end
ps = convert(Matrix{Float64}, deepflatten(convert(Vector{Vector}, ps)))

# f = legendre_transformation(free_energy, NaturalVariable{:P}(pressures))
# interpolator = NDInterpolator{1}(Spline1D)
# f(interpolator)

free_energy + ThermodynamicField{:T,:V}(temperature, volume, ps .* reshape(volume, 1,size(volume)...))
