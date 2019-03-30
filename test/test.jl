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
using CSV

df = CSV.read("f.csv")
temperature = df[:, 1]

# temperature = NaturalVariable{:T}(collect(1000.0:1000:10000))
# volume = NaturalVariable{:V}(collect(100.0:100:500))
# w = NormalMode{:q}(rand(10))
# freq = rand(10, 6) * 1e14

# function generate(func)
#     ret = zeros(5, 10)
#     for (i, t) in enumerate(temperature)
#         for (j, v) in enumerate(volume)
#             qs = QSpaceField(NormalMode{:q}(1:10 |> collect), NormalMode{:s}(1:6 |> collect), func.(t, freq))
#             ret[j, i] = sample_brillouin_zone(w, qs)
#         end
#     end
#     ret
# end

# free_energy = ThermodynamicField(volume, temperature, rand(5, 10))
# internal_energy = ThermodynamicField(volume, temperature, rand(5, 10))
# entropy = ThermodynamicField(volume, temperature, rand(5, 10))
# cv = ThermodynamicField(volume, temperature, rand(5, 10))

fit = fit_energy(BirchMurnaghan3rd(volume[1], 200.0, 4.0), volume, free_energy[:,1])
pressures = eval_pressure(BirchMurnaghan3rd(fit.param[1:end-1]...)).(volume)

f = legendre_transformation(free_energy, NaturalVariable{:P}(pressures))
interpolator = NDInterpolator{1}(Spline1D)
f(interpolator)
