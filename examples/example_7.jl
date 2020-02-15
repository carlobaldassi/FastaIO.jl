# example 7: stream in chunks with write

module FastaIOExample7
using FastaIO
using GZip

# notes about the data:
#   1. the sequence type could be
#      different, e.g. Vector{Uint8},
#      Vector{AminoAcid} etc.
#   2. the only thing which marks descriptions
#      is the '>' sign at the beginning; the
#      sequence data can be split arbitrarily;
#      whitespace is ignored
#
fastadata = [
    ">ITEM1",
    "-----------STVELTKEN-F--D-Q",
    "E--F-V--LI--------D-----F--",
    "P----C-----R----Q---F",
    ">ITEM2",
    "----------MVKEITDAT-F--E-Q-",
    "--L-V--LT--------D-----F--W",
    "-------M-------D-----",
    ">ITEM3",
    "----------NLESVEQFD--------",
    "-------G--K-S--VF--------M-",
    "-D------L------A-----"]

function write_fasta_data(filename::String, data)
    println("Writing to file $filename")

    # note: FastaWriter can take a writable IO descriptor as
    #       argument instead of a filename; also, it has
    #       an optional mode argument which defaults to "w"
    #       (use "a" for appending to an exisitng file).
    #       The do-notation ensures closing of the stream
    #       when finished, which is crucial.
    FastaWriter(filename) do fw
        for item in data
            write(fw, item)
        end
    end

    println("This is the content of the file:")
    gzopen(filename) do f
        println(read(f, String))
    end
end

# note: gzip compression is used because the filename ends with .gz
write_fasta_data(joinpath(@__DIR__, "example_out.fasta.gz"), fastadata)

# remove the file
rm(joinpath(@__DIR__, "example_out.fasta.gz"))

end # module
