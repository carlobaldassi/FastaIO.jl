using Documenter, FastaIO

makedocs(
    modules  = [FastaIO],
    format   = :html,
    sitename = "FastaIO.jl",
    pages    = Any[
        "Home" => "index.md",
       ]
    )

deploydocs(
    repo   = "github.com/carlobaldassi/FastaIO.jl.git",
    target = "build",
    deps = nothing,
    make = nothing,
    julia  = "0.7"
)
