using QuasiHarmonicApprox
using Documenter

makedocs(;
    modules=[QuasiHarmonicApprox],
    authors="Qi Zhang <singularitti@outlook.com>",
    repo="https://github.com/MineralsCloud/QuasiHarmonicApprox.jl/blob/{commit}{path}#L{line}",
    sitename="QuasiHarmonicApprox.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://MineralsCloud.github.io/QuasiHarmonicApprox.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/MineralsCloud/QuasiHarmonicApprox.jl",
)
