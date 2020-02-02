using QuasiHarmonicApproximation
using Documenter

makedocs(;
    modules=[QuasiHarmonicApproximation],
    authors="Qi Zhang <singularitti@outlook.com>",
    repo="https://github.com/MineralsCloud/QuasiHarmonicApproximation.jl/blob/{commit}{path}#L{line}",
    sitename="QuasiHarmonicApproximation.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://MineralsCloud.github.io/QuasiHarmonicApproximation.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/MineralsCloud/QuasiHarmonicApproximation.jl",
)
