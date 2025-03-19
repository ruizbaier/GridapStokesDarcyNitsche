using GridapStokesDarcyNitsche
using Documenter

DocMeta.setdocmeta!(GridapStokesDarcyNitsche, :DocTestSetup, :(using GridapStokesDarcyNitsche); recursive=true)

makedocs(;
    modules=[GridapStokesDarcyNitsche],
    authors="Ricardo Ruiz-Baier <ricardo.ruizbaier@monash.edu>",
    sitename="GridapStokesDarcyNitsche.jl",
    format=Documenter.HTML(;
        canonical="https://ruizbaier.github.io/GridapStokesDarcyNitsche.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/ruizbaier/GridapStokesDarcyNitsche.jl",
    devbranch="main",
)
