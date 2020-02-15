# example 4: read with iterators,
#            custom return type

module FastaIOExample4
using FastaIO

function read_fasta_file(filename::String)
    # We set the reader so that it returns a Vector{Char};
    # If you have BioSeq package, you can use more
    # specialized uints, e.g. Vector{AminoAcid}
    fr = FastaReader{Vector{Char}}(filename)

    for (desc, seq) in fr
        println("num=$(fr.num_parsed) desc=$desc seq=$seq")
    end
    println("read $(fr.num_parsed) entries")
end

read_fasta_file(joinpath(@__DIR__, "example.fasta.gz"))

end # module
