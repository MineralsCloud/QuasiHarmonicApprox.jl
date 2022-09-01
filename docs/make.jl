using QuasiHarmonicApprox
using Documenter

DocMeta.setdocmeta!(QuasiHarmonicApprox, :DocTestSetup, :(using QuasiHarmonicApprox); recursive=true)

makedocs(;
    modules=[QuasiHarmonicApprox],
    authors="singularitti <singularitti@outlook.com>",
    repo="https://github.com/MineralsCloud/QuasiHarmonicApprox.jl/blob/{commit}{path}#{line}",
    sitename="QuasiHarmonicApprox.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://MineralsCloud.github.io/QuasiHarmonicApprox.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
        "Manual" => ["Installation guide" => "installation.md"],
        "API Reference" => [],
        "Developer Docs" => [
            "Contributing" => "developers/contributing.md",
            "Style Guide" => "developers/style.md",
        ],
        "Troubleshooting" => "troubleshooting.md",
    ],
)

deploydocs(;
    repo="github.com/MineralsCloud/QuasiHarmonicApprox.jl",
)
