# example 8: stream one char at a time with write

module FastaIOExample8
using FastaIO
using GZip

# notes about the data:
#   1. the streamed data can be anything as long as
#      it's convertible to an ASCII Char
#   2. descriptions start with '>' and end with a
#      newline; sequence data can be split across
#      multiple lines (they will be reformatted
#      anyways) and end with a newline; whitespace
#      is ignored
#
fastadata = """
    >ITEM1
    -----------STVELTKEN-F--D-Q
    E--F-V--LI--------D-----F--
    P----C-----R----Q---F
    >ITEM2
    ----------MVKEITDAT-F--E-Q-
    --L-V--LT--------D-----F--W
    -------M-------D-----
    >ITEM3
    ----------NLESVEQFD--------
    -------G--K-S--VF--------M-
    -D------L------A-----
    """

function write_fasta_data(filename::String, data)
    println("Writing to file $filename")

    # note: FastaWriter can take a writable IO descriptor as
    #       argument instead of a filename; also, it has
    #       an optional mode argument which defaults to "w"
    #       (use "a" for appending to an exisitng file).
    #       The do-notation ensures closing of the stream
    #       when finished, which is crucial.
    FastaWriter(filename) do fw
        for char in data
            write(fw, char)
        end
    end

    println("This is the content of the file:")
    gzopen(filename) do f
        println(readall(f))
    end
end

examples_dir = dirname(Base.source_path())
# note: gzip compression is used because the filename ends with .gz
write_fasta_data(joinpath(examples_dir, "example_out.fasta.gz"), fastadata)

# remove the file
rm(joinpath(examples_dir, "example_out.fasta.gz"))

end
