using Documenter, FastaIO

makedocs(
    modules  = [FastaIO],
    format = Documenter.HTML(prettyurls = "--local" ∉ ARGS),
    sitename = "FastaIO.jl",
    pages    = Any[
        "Home" => "index.md",
       ]
    )

deploydocs(
    repo   = "github.com/carlobaldassi/FastaIO.jl.git",
)
