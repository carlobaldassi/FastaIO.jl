# example 3: read all at once with readall (or readfasta)
#            default return type (String)

module FastaIOExample3
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

    out = readall(fr)

    # alternatively:
    # out = readfasta(filename)

    n = 0
    for (desc, seq) in out
        n += 1
        println("num=$n desc=$desc seq=$seq")
    end
    # here of course n == fr.num_parsed == length(out)
    println("read $n entries")
end

examples_dir = dirname(Base.source_path())
read_fasta_file(joinpath(examples_dir, "example.fasta.gz"))
end
