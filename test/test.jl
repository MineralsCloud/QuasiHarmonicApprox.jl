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
using QuasiHarmonicApproximation.Thermodynamics
using EquationsOfState

freq = rand(10, 6, 5, 10) * 1e14
temperature = collect(1000.0:1000:10000)
volume = collect(100.0:100:500)
w = NormalMode{:q}(rand(10))

function generate(f, func)
    ret = zeros(5, 10)
    for (i, t) in enumerate(temperature)
        for (j, v) in enumerate(volume)
            ret[j, i] = sample_brillouin_zone(w, QSpaceField{:q, :s}(1:10 |> collect, 1:6 |> collect, func.(t, f[:, :, j, i])))
        end
    end
    ret
end

free_energy = ThermodynamicField{:V, :T}(volume, temperature, generate(freq, subsystem_free_energy))
internal_energy = ThermodynamicField{:V, :T}(volume, temperature, generate(freq, subsystem_internal_energy))
entropy = ThermodynamicField{:V, :T}(volume, temperature, generate(freq, subsystem_entropy))
cv = ThermodynamicField{:V, :T}(volume, temperature, generate(freq, subsystem_volumetric_specific_heat))

fit = fit_energy(BirchMurnaghan3rd(volume[1], 200.0, 4.0), volume, free_energy.values[:,1])
pressures = eval_pressure(BirchMurnaghan3rd(fit.param[1:end-1]...)).(volume)

legendre_transformation(free_energy, :T)
