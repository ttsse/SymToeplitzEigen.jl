using ToeplitzCUDA_Eigen_Refinement
using Documenter

DocMeta.setdocmeta!(ToeplitzCUDA_Eigen_Refinement, :DocTestSetup, :(using ToeplitzCUDA_Eigen_Refinement); recursive=true)

makedocs(;
    modules=[ToeplitzCUDA_Eigen_Refinement],
    authors="David Meadon",
    repo="https://github.com/Dedavond/ToeplitzCUDA_Eigen_Refinement.jl/blob/{commit}{path}#{line}",
    sitename="ToeplitzCUDA_Eigen_Refinement.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://Dedavond.github.io/ToeplitzCUDA_Eigen_Refinement.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/Dedavond/ToeplitzCUDA_Eigen_Refinement.jl",
    devbranch="main",
)
