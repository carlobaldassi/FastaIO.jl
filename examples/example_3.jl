# example 3: read all at once with readfasta
#            default return type (String)

module FastaIOExample3
using FastaIO

function read_fasta_file(filename::String)
    out = readfasta(filename)

    n = 0
    for (desc, seq) in out
        n += 1
        println("num=$n desc=$desc seq=$seq")
    end
    # here of course n == fr.num_parsed == length(out)
    println("read $n entries")
end

read_fasta_file(joinpath(@__DIR__, "example.fasta.gz"))

end # module
