# example 1: read with iterators,
#            default return type (String)

module FastaIOExample1
using FastaIO

function read_fasta_file(filename::String)
    # note: FastaReader can take a readable IO descriptor as
    #       argument instead of a filename; it also supports
    #       do-notation:
    #
    #       FastaReader(filename) do fr
    #         # do something
    #       end
    fr = FastaReader(filename)

    for (desc, seq) in fr
        println("num=$(fr.num_parsed) desc=$desc seq=$seq")
    end
    println("read $(fr.num_parsed) entries")
end

read_fasta_file(joinpath(@__DIR__, "example.fasta.gz"))

end # module
