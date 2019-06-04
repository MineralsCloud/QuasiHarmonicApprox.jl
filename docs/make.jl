using Documenter, QuasiHarmonicApproximation

makedocs(;
    modules=[QuasiHarmonicApproximation],
    format=Documenter.HTML(),
    pages=[
        "Home" => "index.md",
    ],
    repo="https://github.com/singularitti/QuasiHarmonicApproximation.jl/blob/{commit}{path}#L{line}",
    sitename="QuasiHarmonicApproximation.jl",
    authors="Qi Zhang <singularitti@outlook.com>",
    assets=String[],
)

deploydocs(;
    repo="github.com/singularitti/QuasiHarmonicApproximation.jl",
)
