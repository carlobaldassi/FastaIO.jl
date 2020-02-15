# example 2: read with readentry
#            default return type (String)

module FastaIOExample2
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

    while !eof(fr)
        desc, seq = readentry(fr)
        println("num=$(fr.num_parsed) desc=$desc seq=$seq")
    end
    println("read $(fr.num_parsed) entries")
end

read_fasta_file(joinpath(@__DIR__, "example.fasta.gz"))

end # module
