# example 6: write one entry at a time with writeentry

module FastaIOExample6
using FastaIO
using GZip

# note about the data:
#   the sequence type could be
#   different, e.g. Vector{Uint8},
#   Vector{AminoAcid} etc.
fastadata = [
    ("ITEM1","""
     -----------STVELTKEN-F--D-Q
     E--F-V--LI--------D-----F--
     P----C-----R----Q---F
     """),
    ("ITEM2","""
     ----------MVKEITDAT-F--E-Q-
     --L-V--LT--------D-----F--W
     -------M-------D-----
     """),
    ("ITEM3","""
     ----------NLESVEQFD--------
     -------G--K-S--VF--------M-
     -D------L------A-----
     """)]

function write_fasta_data(filename::String, data)
    println("Writing $(length(data)) entries to file $filename")

    # note: FastaWriter can take a writable IO descriptor as
    #       argument instead of a filename; also, it has
    #       an optional mode argument which defaults to "w"
    #       (use "a" for appending to an exisitng file).
    #       The do-notation ensures closing of the stream
    #       when finished, which is crucial.
    FastaWriter(filename) do fw
        for (desc, seq) in data
            writeentry(fw, desc, seq)
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
