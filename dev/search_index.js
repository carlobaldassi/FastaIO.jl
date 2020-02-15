var documenterSearchIndex = {"docs":
[{"location":"#FastaIO.jl-—-FASTA-file-reader-and-writer-module-1","page":"Home","title":"FastaIO.jl — FASTA file reader and writer module","text":"","category":"section"},{"location":"#","page":"Home","title":"Home","text":"CurrentModule = FastaIO","category":"page"},{"location":"#","page":"Home","title":"Home","text":"This module provides ways to parse and write files in FASTA format in Julia. It is designed to be lightweight and fast; the parsing method is inspired by kseq.h. It can read and write files on the fly, keeping only one entry at a time in memory, and it can read and write gzip-compressed files.","category":"page"},{"location":"#","page":"Home","title":"Home","text":"Here is a quick example for reading a file:","category":"page"},{"location":"#","page":"Home","title":"Home","text":"julia> using FastaIO\n\njulia> FastaReader(\"somefile.fasta\") do fr\n           for (desc, seq) in fr\n               println(\"$desc : $seq\")\n           end\n       end","category":"page"},{"location":"#","page":"Home","title":"Home","text":"And for writing:","category":"page"},{"location":"#","page":"Home","title":"Home","text":"julia> using FastaIO\n\njulia> FastaWriter(\"somefile.fasta\") do fw\n           for s in [\">GENE1\", \"GCATT\", \">GENE2\", \"ATTAGC\"]\n               write(fw, s)\n           end\n       end","category":"page"},{"location":"#Installation-and-usage-1","page":"Home","title":"Installation and usage","text":"","category":"section"},{"location":"#","page":"Home","title":"Home","text":"To install the module, use Julia's package manager: start pkg mode by pressing <kbd>]</kbd> and then enter:","category":"page"},{"location":"#","page":"Home","title":"Home","text":"(v1.3) pkg> add FastaIO","category":"page"},{"location":"#","page":"Home","title":"Home","text":"Dependencies will be installed automatically. The module can then be loaded like any other Julia module:","category":"page"},{"location":"#","page":"Home","title":"Home","text":"julia> using FastaIO","category":"page"},{"location":"#Introductory-notes-1","page":"Home","title":"Introductory notes","text":"","category":"section"},{"location":"#","page":"Home","title":"Home","text":"For both reading and writing, there are quick methods to read/write all the data at once: readfasta and writefasta. These, however, require all the data to be stored in memory at once, which may be impossible or undesirable for very large files. Therefore, for both reading and writing, the preferred way is actually to use specialized types, FastaReader and FastaWriter, which have the ability to process one entry (description + sequence data) at a time (the writer can actually process one char at a time); however, note that these two object types are not symmetric: the reader acts as an iterable object, while the writer behaves similarly to an IO stream.","category":"page"},{"location":"#The-FASTA-format-1","page":"Home","title":"The FASTA format","text":"","category":"section"},{"location":"#","page":"Home","title":"Home","text":"The FASTA format which is assumed by this module is as follows:","category":"page"},{"location":"#","page":"Home","title":"Home","text":"description lines must start with a > character, and cannot be empty\nonly one description line per entry is allowed\nall characters must be ASCII\nwhitespace is not allowed within sequence data (except for newlines) and at the beginning or end of the description\nEmpty lines are ignored (note however that lines containing whitespace will still trigger an error)","category":"page"},{"location":"#","page":"Home","title":"Home","text":"When writing, description lines longer than 80 characters will trigger a warning message; sequence data is formatted in lines of 80 characters each; extra whitespace is silently discarded. No other restriction is put on the content of the sequence data, except that the > character is forbidden.","category":"page"},{"location":"#","page":"Home","title":"Home","text":"When reading, almost no explicit checks are performed to test that the data actually conforms to these specifications.","category":"page"},{"location":"#The-sequence-storage-type-1","page":"Home","title":"The sequence storage type","text":"","category":"section"},{"location":"#","page":"Home","title":"Home","text":"When reading FASTA files, the container type used to store the sequence data can be chosen (as an optional argument to readfasta or as a parametric type of FastaReader). The default is String, which is the most memory-efficient and the fastest; another performance-optimal option is Vector{UInt8}, which is a less friendly representation, but has the advantage of being mutable. Any other container T for which convert(::Type{T}, ::Vector{UInt8}) is defined can be used (e.g. Vector{Char}, or a more specialized Vector{AminoAcid} if you use BioSeq, but the conversion will generally slightly reduce the performance.","category":"page"},{"location":"#Reading-files-1","page":"Home","title":"Reading files","text":"","category":"section"},{"location":"#","page":"Home","title":"Home","text":"readfasta","category":"page"},{"location":"#FastaIO.readfasta","page":"Home","title":"FastaIO.readfasta","text":"readfasta(file::Union{String,IO}, [sequence_type::Type = String])\n\nThis function parses a whole FASTA file at once and stores it into memory. The result is a Vector{Any} whose elements are tuples consisting of (description, sequence), where description is a String and sequence contains the sequence data, stored in a container type defined by the sequence_type optional argument (see The sequence storage type section for more information).\n\n\n\n\n\n","category":"function"},{"location":"#","page":"Home","title":"Home","text":"FastaReader(::Union{AbstractString,IO})","category":"page"},{"location":"#FastaIO.FastaReader-Tuple{Union{AbstractString, IO}}","page":"Home","title":"FastaIO.FastaReader","text":"FastaReader{T}(file::Union{AbstractString,IO})\n\nThis creates an object which is able to parse FASTA files, one entry at a time. file can be a plain text file or a gzip-compressed file (it will be autodetected from the content). The type T determines the output type of the sequences (see The sequence storage type section for more information) and it defaults to String.\n\nThe data can be read out by iterating the FastaReader object:\n\nfor (name, seq) in FastaReader(\"somefile.fasta\")\n    # do something with name and seq\nend\n\nAs shown, the iterator returns a tuple containing the description (always a String) and the data (whose type is set when creating the FastaReader object (e.g. FastaReader{Vector{UInt8}}(filename)).\n\nThe FastaReader type has a field num_parsed which contains the number of entries parsed so far.\n\nOther ways to read out the data are via the readentry and readfasta functions.\n\n\n\n\n\n","category":"method"},{"location":"#","page":"Home","title":"Home","text":"FastaReader(::Function, ::AbstractString)","category":"page"},{"location":"#FastaIO.FastaReader-Tuple{Function,AbstractString}","page":"Home","title":"FastaIO.FastaReader","text":"FastaReader(f::Function, filename::AbstractString, [sequence_type::Type = String])\n\nThis format of the constructor is useful for do-notation, i.e.:\n\nFastaReader(filename) do fr\n    # read out the data from fr, e.g.\n    for (name, seq) in fr\n        # do something with name and seq\n    end\nend\n\nwhich ensures that the close function is called and is thus recommended (otherwise the file is closed by the garbage collector when the FastaReader object goes out of scope).\n\n\n\n\n\n","category":"method"},{"location":"#","page":"Home","title":"Home","text":"readentry","category":"page"},{"location":"#FastaIO.readentry","page":"Home","title":"FastaIO.readentry","text":"readentry(fr::FastaReader)\n\nThis function can be used to read entries one at a time:\n\nfr = FastaReader(\"somefile.fasta\")\nname, seq = readentry(fr)\n\nSee also the eof function.\n\n\n\n\n\n","category":"function"},{"location":"#","page":"Home","title":"Home","text":"rewind(::FastaReader)","category":"page"},{"location":"#FastaIO.rewind-Tuple{FastaReader}","page":"Home","title":"FastaIO.rewind","text":"rewind(fr::FastaReader)\n\nThis function rewinds the reader, so that it can restart the parsing again without closing and re-opening it. It also resets the value of the num_parsed field.\n\n\n\n\n\n","category":"method"},{"location":"#","page":"Home","title":"Home","text":"eof(::FastaReader)","category":"page"},{"location":"#Base.eof-Tuple{FastaReader}","page":"Home","title":"Base.eof","text":"eof(fr::FastaReader)\n\nThis function extends Base.eof and tests for end-of-file condition; it is useful when using readentry:\n\nfr = FastaReader(\"somefile.fasta\")\nwhile !eof(fr)\n   name, seq = readentry(fr)\n   # do something\nend\nclose(fr)\n\n\n\n\n\n","category":"method"},{"location":"#","page":"Home","title":"Home","text":"close(::FastaReader)","category":"page"},{"location":"#Base.close-Tuple{FastaReader}","page":"Home","title":"Base.close","text":"close(fr::FastaReader)\n\nThis function extends Base.close and closes the stream associated with the FastaReader; the reader must not be used any more after this function is called.\n\n\n\n\n\n","category":"method"},{"location":"#Writing-files-1","page":"Home","title":"Writing files","text":"","category":"section"},{"location":"#","page":"Home","title":"Home","text":"writefasta(::AbstractString, ::Any)","category":"page"},{"location":"#FastaIO.writefasta-Tuple{AbstractString,Any}","page":"Home","title":"FastaIO.writefasta","text":"writefasta(filename::String, data, [mode::String = \"w\"])\n\nThis function dumps data to a FASTA file, auto-formatting it so to follow the specifications detailed in the section titled The FASTA format. The data can be anything which is iterable and which produces (description, sequence) tuples upon iteration, where the description must be convertible to a String and the sequence can be any iterable object which yields elements convertible to ASCII characters (e.g. a String, a Vector{UInt8} etc.).\n\nExamples:\n\nwritefasta(\"somefile.fasta\", [(\"GENE1\", \"GCATT\"), (\"GENE2\", \"ATTAGC\")])\nwritefasta(\"somefile.fasta\", [\"GENE1\" => \"GCATT\", \"GENE2\" => \"ATTAGC\"])\n\nIf the filename ends with .gz, the result will be a gzip-compressed file.\n\nThe mode flag determines how the filename is open; use \"a\" to append the data to an existing file.\n\n\n\n\n\n","category":"method"},{"location":"#","page":"Home","title":"Home","text":"writefasta(::IO, ::Any)","category":"page"},{"location":"#FastaIO.writefasta-Tuple{IO,Any}","page":"Home","title":"FastaIO.writefasta","text":"writefasta([io::IO = STDOUT], data)\n\nThis version of the function writes to an already opened IO stream, defaulting to STDOUT.\n\n\n\n\n\n","category":"method"},{"location":"#","page":"Home","title":"Home","text":"FastaWriter","category":"page"},{"location":"#FastaIO.FastaWriter","page":"Home","title":"FastaIO.FastaWriter","text":"FastaWriter(filename::AbstractString, [mode::String = \"w\"])\nFastaWriter([io::IO = STDOUT])\nFastaWriter(f::Function, args...)\n\nThis creates an object which is able to write formatted FASTA files which conform to the specifications detailed in the section titled The FASTA format, via the write and writeentry functions.\n\nThe third form allows to use do-notation:\n\nFastaWriter(\"somefile.fasta\") do fw\n    # write the file\nend\n\nwhich is strongly recommended since it ensures that the close function is called at the end of writing: this is crucial, as failing to do so may result in incomplete files (this is done by the finalizer, so it will still happen automatically if the FastaWriter object goes out of scope and is garbage-collected, but there is no guarantee that this will happen if Julia exits).\n\nIf the filename ends with .gz, the result will be gzip-compressed.\n\nThe mode flag can be used to set the opening mode of the file; use \"a\" to append to an existing file.\n\nThe FastaWriter object has an entry::Int field which stores the number of the entry which is currently being written.\n\n\n\n\n\n","category":"type"},{"location":"#","page":"Home","title":"Home","text":"writeentry","category":"page"},{"location":"#FastaIO.writeentry","page":"Home","title":"FastaIO.writeentry","text":"writeentry(fw::FastaWriter, description::AbstractString, sequence)\n\nThis function writes one entry to the FASTA file, following the specifications detailed in the section titled The FASTA format. The description is without the initial '>' character. The sequence can be any iterable object whose elements are convertible to ASCII characters.\n\nExample:\n\nFastaWriter(\"somefile.fasta\") do fw\n    for (desc,seq) in [(\"GENE1\", \"GCATT\"), (\"GENE2\", \"ATTAGC\")]\n        writeentry(fw, desc, seq)\n    end\nend\n\n\n\n\n\n","category":"function"},{"location":"#","page":"Home","title":"Home","text":"write(::FastaWriter, ::Any)","category":"page"},{"location":"#Base.write-Tuple{FastaWriter,Any}","page":"Home","title":"Base.write","text":"write(fw::FastaWriter, item)\n\nThis function extends Base.write and streams items to a FASTA file, which will be formatted according to the specifications detailed in the section titled The FASTA format.\n\nWhen using this method, description lines are marked by the fact that they begin with a '>' character; anything else is assumed to be part of the sequence data.\n\nIf item is a Vector, write will be called iteratively over it; if it is a String, a newline will be appended to it and it will be dumped. For example the following code:\n\nFastaWriter(\"somefile.fasta\") do fw\n    for s in [\">GENE1\", \"GCA\", \"TTT\", \">GENE2\", \"ATTAGC\"]\n        write(fw, s)\n    end\nend\n\nwill result in the file:\n\n>GENE1\nGCATTT\n>GENE2\nATTAGC\n\nIf item is not a Vector nor a String, it must be convertible to an ASCII character, and it will be piped into the file. For example the following code:\n\ndata = \"\"\"\n  >GENE1\n  GCA\n  TTT\n  >GENE2\n  ATT\n  AGC\n  \"\"\"\n\nFastaWriter(\"somefile.fasta\") do fw\n    for ch in data\n        write(fw, ch)\n    end\nend\n\nwill result in the same file as above.\n\n\n\n\n\n","category":"method"},{"location":"#","page":"Home","title":"Home","text":"close(::FastaWriter)","category":"page"},{"location":"#Base.close-Tuple{FastaWriter}","page":"Home","title":"Base.close","text":"close(fw::FastaWriter)\n\nThis function extends Base.close and it should always be explicitly used for finalizing the FastaWriter once the writing has finished, unless the do-notation is used when creating it.\n\n\n\n\n\n","category":"method"}]
}
