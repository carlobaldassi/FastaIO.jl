"""
This module provides ways to parse and write files in
[FASTA format](http://en.wikipedia.org/wiki/FASTA_format) in Julia.
It is designed to be lightweight and fast; the parsing method is inspired by
[kseq.h](http://lh3lh3.users.sourceforge.net/kseq.shtml). It can read and write
files on the fly, keeping only one entry at a time in memory, and it can read and
write gzip-compressed files.

See [`FastaReader`](@ref), [`FastaWriter`](@ref), [`readfasta`](@ref),
[`writefasta`](@ref), [`readentry`](@ref) and [`writeentry`](@ref).
"""
module FastaIO

using GZip

export
    FastaReader,
    readentry,
    rewind,
    readfasta,
    FastaWriter,
    writeentry,
    writefasta

import Base: close, show, eof, write, iterate

const fasta_buffer_size = 4096

mutable struct FastaReader{T}
    # public but read-only
    num_parsed::Int        # number of parsed entries so far
    # private
    f::IO
    is_eof::Bool           # did we reach end of file?
    rbuffer::Vector{UInt8} # read buffer
    rbuf_sz::Int           # read buffer size
    rbuf_pos::Int          # read buffer cursor
    lbuffer::Vector{UInt8} # line buffer
    lbuf_sz::Int           # line buffer size
    mbuffer::Vector{UInt8} # multi-line buffer
    mbuf_sz::Int           # multi-line buffer size
    own_f::Bool
    function FastaReader{T}(filename::AbstractString) where {T}
        fr = new{T}(0, gzopen(filename), false,
                    Array{UInt8}(undef, fasta_buffer_size), 0, 0,
                    Array{UInt8}(undef, fasta_buffer_size), 0,
                    Array{UInt8}(undef, fasta_buffer_size), 0,
                    true)
        finalizer(close, fr)
        return fr
    end
    function FastaReader{T}(io::IO) where {T}
        new{T}(0, io, false,
               Array{UInt8}(undef, fasta_buffer_size), 0, 0,
               Array{UInt8}(undef, fasta_buffer_size), 0,
               Array{UInt8}(undef, fasta_buffer_size), 0,
               false)
    end
end

"""
    FastaReader{T}(file::Union{AbstractString,IO})

This creates an object which is able to parse FASTA files, one entry at a time. `file` can be a plain text
file or a gzip-compressed file (it will be autodetected from the content). The type `T` determines the output
type of the sequences (see [The sequence storage type](@ref) section for more information) and it defaults to
`String`.

The data can be read out by iterating the `FastaReader` object:

```julia
for (name, seq) in FastaReader("somefile.fasta")
    # do something with name and seq
end
```

As shown, the iterator returns a tuple containing the description (always a `String`) and the data (whose
type is set when creating the `FastaReader` object (e.g. `FastaReader{Vector{UInt8}}(filename)`).

The `FastaReader` type has a field `num_parsed` which contains the number of entries parsed so far.

Other ways to read out the data are via the [`readentry`](@ref) and [`readfasta`](@ref) functions.
"""
FastaReader(file::Union{AbstractString,IO}) = FastaReader{String}(file)


"""
    FastaReader(f::Function, filename::AbstractString, [sequence_type::Type = String])

This format of the constructor is useful for do-notation, i.e.:

```julia
FastaReader(filename) do fr
    # read out the data from fr, e.g.
    for (name, seq) in fr
        # do something with name and seq
    end
end
```

which ensures that the [`close`](@ref) function is called and is thus recommended (otherwise
the file is closed by the garbage collector when the `FastaReader` object goes out of scope).
"""
function FastaReader(f::Function, filename::AbstractString, T::Type=String)
    fr = FastaReader{T}(filename)
    try
        f(fr)
    finally
        close(fr)
    end
end

"""
    close(fr::FastaReader)

This function extends `Base.close` and closes the stream associated with the `FastaReader`; the reader must not
be used any more after this function is called.
"""
close(fr::FastaReader) = fr.own_f && close(fr.f)

"""
    rewind(fr::FastaReader)

This function rewinds the reader, so that it can restart the parsing again without closing and re-opening it.
It also resets the value of the `num_parsed` field.
"""
function rewind(fr::FastaReader)
    seek(fr.f, 0)
    fr.is_eof = false
    fr.num_parsed = 0
    fr.rbuf_sz = 0
    fr.rbuf_pos = 0
    fr.lbuf_sz = 0
    fr.mbuf_sz = 0
    return
end

read_chunk_ll(fr::FastaReader, s::GZipStream) = gzread(s, pointer(fr.rbuffer), fasta_buffer_size)
read_chunk_ll(fr::FastaReader, s::IOStream) =
    ccall(:ios_readall, UInt, (Ptr{Nothing}, Ptr{Nothing}, UInt), fr.f.ios, fr.rbuffer, fasta_buffer_size)
function read_chunk_ll(fr::FastaReader, s::IO)
    ret = 0
    while !eof(fr.f) && ret < fasta_buffer_size
        ret += 1
        fr.rbuffer[ret] = read(fr.f, UInt8)
    end
    return ret
end

function read_chunk(fr::FastaReader)
    fr.is_eof && return
    ret = read_chunk_ll(fr, fr.f)
    ret == -1 && error("read failure")
    fr.rbuf_sz = ret
    fr.rbuf_pos = 1
    ret == 0 && (fr.is_eof = true)
    return
end

function readline(fr::FastaReader)
    fr.lbuf_sz = 0
    found = false
    while !fr.is_eof
        fr.rbuf_pos == 0 && read_chunk(fr::FastaReader)
        i = fr.rbuf_pos
        cr = false
        while i <= fr.rbuf_sz
            c = fr.rbuffer[i]
            if c == UInt8('\n')
                found = true
                break
            else
                cr = (c == UInt8('\r'))
            end
            i += 1
        end
        i -= 1 + cr
        chunk_len = i - fr.rbuf_pos + 1
        free_sbuf = length(fr.lbuffer) - fr.lbuf_sz
        gap = chunk_len - free_sbuf
        gap > 0 && resize!(fr.lbuffer, length(fr.lbuffer) + gap)

        copyto!(fr.lbuffer, fr.lbuf_sz + 1, fr.rbuffer, fr.rbuf_pos, chunk_len)
        fr.lbuf_sz += chunk_len

        i += 2 + cr
        i > fr.rbuf_sz && (i = 0)
        fr.rbuf_pos = i
        found && break
    end
    # if we found an empty line, ignore it and run again
    !fr.is_eof && fr.lbuf_sz == 0 && return readline(fr)
    return
end

# gets the ID from lbuffer and reads FASTA entry into mbuffer
function _next_step(fr::FastaReader)
    if fr.lbuffer[1] != UInt8('>')
        error("invalid FASTA file: description does not start with '>'")
    end
    if fr.lbuf_sz == 1
        error("invalid FASTA file: empty description")
    end
    name = String(fr.lbuffer[2:fr.lbuf_sz])
    isascii(name) || error("invalid non-ASCII description in FASTA file")
    fr.mbuf_sz = 0
    while true
        readline(fr)
        if fr.lbuf_sz == 0 || fr.lbuffer[1] == UInt8('>')
            break
        end
        gap = fr.lbuf_sz - (length(fr.mbuffer) - fr.mbuf_sz)
        gap > 0 && resize!(fr.mbuffer, length(fr.mbuffer) + gap)
        copyto!(fr.mbuffer, fr.mbuf_sz + 1, fr.lbuffer, 1, fr.lbuf_sz)
        fr.mbuf_sz += fr.lbuf_sz
    end
    return name
end

# extracts the sequence of the FASTA entry from the mbuffer, the result is of type T
_next_seq(fr::FastaReader{Vector{UInt8}}) = fr.mbuffer[1:fr.mbuf_sz]
_next_seq(fr::FastaReader{String}) = ccall(:jl_pchar_to_string, Ref{String}, (Ptr{UInt8},Int), fr.mbuffer, fr.mbuf_sz)
_next_seq(fr::FastaReader{T}) where T = T(fr.mbuffer[1:fr.mbuf_sz])

function _next(fr::FastaReader{T}) where T
    name = _next_step(fr)
    seq = _next_seq(fr)::T
    fr.num_parsed += 1
    return (name, seq)
end

Base.eltype(fr::FastaReader{T}) where T = Tuple{String, T}
Base.IteratorSize(fr::FastaReader) = Base.SizeUnknown()

function iterate(fr::FastaReader)
    rewind(fr)
    readline(fr)
    fr.lbuf_sz == 0 && error("empty FASTA file")
    return fr.is_eof ? nothing : (_next(fr), nothing)
end
function iterate(fr::FastaReader, x::Nothing)
    fr.is_eof && return nothing
    return _next(fr), nothing
end

"""
    readentry(fr::FastaReader)

This function can be used to read entries one at a time:

```julia
fr = FastaReader("somefile.fasta")
name, seq = readentry(fr)
```

See also the [`eof`](@ref) function.
"""
function readentry(fr::FastaReader)
    fr.is_eof && throw(EOFError())
    if fr.num_parsed == 0
        readline(fr)
        fr.lbuf_sz == 0 && error("empty FASTA file")
    end
    item = _next(fr)
    return item
end

"""
    eof(fr::FastaReader)

This function extends `Base.eof` and tests for end-of-file condition; it is useful when using
[`readentry`](@ref):

```julia
fr = FastaReader("somefile.fasta")
while !eof(fr)
   name, seq = readentry(fr)
   # do something
end
close(fr)
```
"""
eof(fr::FastaReader) = fr.is_eof

show(io::IO, fr::FastaReader{T}) where {T} =
    print(io, "FastaReader(input=\"$(fr.f)\", out_type=$T, num_parsed=$(fr.num_parsed), eof=$(fr.is_eof))")

"""
    readfasta(file::Union{String,IO}, [sequence_type::Type = String])

This function parses a whole FASTA file at once and stores it into memory. The result is a `Vector{Any}`
whose elements are tuples consisting of `(description, sequence)`, where `description` is a
`String` and `sequence` contains the sequence data, stored in a container type defined by
the `sequence_type` optional argument (see [The sequence storage type](@ref) section for more information).
"""
readfasta(filename::AbstractString, ::Type{T}=String) where T =
    gzopen(io -> readfasta(io, T), filename)
readfasta(io::IO, ::Type{T}=String) where T = collect(FastaReader{T}(io))

mutable struct FastaWriter
    f::IO
    in_seq::Bool
    entry_chars::Int
    desc_chars::Int
    parsed_nl::Bool
    pos::Int
    entry::Int
    own_f::Bool
    at_start::Bool
    function FastaWriter(io::IO)
        fw = new(io, false, 0, 0, false, 0, 1, false, true)
        finalizer(close, fw)
        return fw
    end
    function FastaWriter(filename::AbstractString, mode::AbstractString = "w")
        fopen = endswith(filename, ".gz") ? gzopen : open
        fw = new(fopen(filename, mode), false, 0, 0, false, 0, 1, true, true)
        finalizer(close, fw)
        return fw
    end
end

FastaWriter() = FastaWriter(STDOUT)

"""
    FastaWriter(filename::AbstractString, [mode::String = "w"])
    FastaWriter([io::IO = STDOUT])
    FastaWriter(f::Function, args...)

This creates an object which is able to write formatted FASTA files which conform to the specifications
detailed in the section titled [The FASTA format](@ref), via the [`write`](@ref) and [`writeentry`](@ref)
functions.

The third form allows to use do-notation:

```julia
FastaWriter("somefile.fasta") do fw
    # write the file
end
```

which is strongly recommended since it ensures that the [`close`](@ref) function is called at the end of
writing: this is crucial, as failing to do so may result in incomplete files (this is done by the
finalizer, so it will still happen automatically if the `FastaWriter` object goes out of scope and is
garbage-collected, but there is no guarantee that this will happen if Julia exits).

If the `filename` ends with `.gz`, the result will be gzip-compressed.

The `mode` flag can be used to set the opening mode of the file; use `"a"` to append to an existing
file.

The `FastaWriter` object has an `entry::Int` field which stores the number of the entry which is
currently being written.
"""
function FastaWriter(f::Function, args...)
    fw = FastaWriter(args...)
    try
        f(fw)
    finally
        close(fw)
    end
end

"""
    write(fw::FastaWriter, item)

This function extends `Base.write` and streams items to a FASTA file, which will be formatted
according to the specifications detailed in the section titled [The FASTA format](@ref).

When using this method, description lines are marked by the fact that they begin with a `'>'`
character; anything else is assumed to be part of the sequence data.

If `item` is a `Vector`, `write` will be called iteratively over it; if it is a `String`,
a newline will be appended to it and it will be dumped. For example the following code:

```julia
FastaWriter("somefile.fasta") do fw
    for s in [">GENE1", "GCA", "TTT", ">GENE2", "ATTAGC"]
        write(fw, s)
    end
end
```

will result in the file:

```text
>GENE1
GCATTT
>GENE2
ATTAGC
```

If `item` is not a `Vector` nor a `String`, it must be convertible to an ASCII character, and
it will be piped into the file. For example the following code:

```julia
data = \"\"\"
  >GENE1
  GCA
  TTT
  >GENE2
  ATT
  AGC
  \"\"\"

FastaWriter("somefile.fasta") do fw
    for ch in data
        write(fw, ch)
    end
end
```

will result in the same file as above.
"""
function write(fw::FastaWriter, c)
    ch = convert(Char, c)
    isascii(ch) || error("invalid (non-ASCII) character: $c (entry $(fw.entry) of FASTA input)")
    if ch == '\n' && !fw.at_start
        fw.parsed_nl = true
        if !fw.in_seq
            fw.desc_chars == 1 && error("empty description (entry $(fw.entry) of FASTA input")
            write(fw.f, '\n')
            fw.pos = 0
            fw.in_seq = true
        end
    end
    isspace(ch) && (fw.at_start || fw.in_seq || fw.desc_chars <= 1) && return
    fw.at_start && ch != '>' && error("no desctiption given (entry $(fw.entry) of FASTA input")
    fw.at_start = false
    if fw.parsed_nl
        @assert fw.in_seq
        if ch == '>'
            fw.entry_chars > 0 || error("description must span a single line (entry $(fw.entry) of FASTA input)")
            write(fw.f, '\n')
            fw.in_seq = false
            fw.pos = 0
            fw.entry += 1
            fw.entry_chars = 0
            fw.desc_chars = 0
        end
    elseif fw.in_seq && ch == '>'
        error("character '>' not allowed in sequence data (entry $(fw.entry) of FASTA input)")
    end
    if fw.pos == 80
        if !fw.in_seq
            @warn("description line longer than 80 characters (entry $(fw.entry) of FASTA input)")
        else
            write(fw.f, '\n')
            fw.pos = 0
        end
    end
    write(fw.f, ch)
    fw.pos += 1
    if fw.in_seq
        fw.entry_chars += 1
    else
        fw.desc_chars += 1
    end
    fw.parsed_nl = false
    return
end

function write(fw::FastaWriter, s::Vector)
    for c in s
        write(fw, c)
    end
end

function write(fw::FastaWriter, s::AbstractString)
    for c in s
        write(fw, c)
    end
    write(fw, '\n')
end

"""
    writeentry(fw::FastaWriter, description::AbstractString, sequence)

This function writes one entry to the FASTA file, following the specifications detailed in
the section titled [The FASTA format](@ref). The `description` is without the initial `'>'` character.
The `sequence` can be any iterable object whose elements are convertible to ASCII characters.

Example:

```julia
FastaWriter("somefile.fasta") do fw
    for (desc,seq) in [("GENE1", "GCATT"), ("GENE2", "ATTAGC")]
        writeentry(fw, desc, seq)
    end
end
```
"""
function writeentry(fw::FastaWriter, desc::AbstractString, seq)
    !fw.at_start && write(fw, '\n')
    desc = strip(String(desc))
    isascii(desc) || error("description must be ASCII (entry $(fw.entry+1) of FASTA input)")
    if findfirst(==('\n'), desc) ≢ nothing
        error("newlines are not allowed within description (entry $(fw.entry+1) of FASTA input)")
    end
    write(fw, '>')
    write(fw, desc)
    write(fw, '\n')
    fw.entry_chars = writefastaseq(fw.f, seq, fw.entry, false)
    fw.in_seq = true
    fw.parsed_nl = false
    fw.pos = 0
    fw.entry_chars > 0 || error("empty sequence data (entry $(fw.entry) of FASTA input)")
    return
end

"""
    close(fw::FastaWriter)

This function extends `Base.close` and it should always be explicitly used for finalizing
the `FastaWriter` once the writing has finished, unless the do-notation is used when creating
it.
"""
function close(fw::FastaWriter)
    try
        write(fw.f, '\n')
        flush(fw.f)
    catch err
        err isa EOFError || rethrow(err)
    end
    fw.pos = 0
    fw.parsed_nl = true
    fw.own_f && close(fw.f)
    return
end

function show(io::IO, fw::FastaWriter)
    print(io, "FastaWriter(input=\"$(fw.f)\", entry=$(fw.entry)")
end

function writefastaseq(io::IO, seq, entry::Int, nl::Bool = true)
    i = 0
    entry_chars = 0
    for c in seq
        if i == 80
            write(io, '\n')
            i = 0
        end
        ch = convert(Char, c)
        isascii(ch) || error("invalid (non-ASCII) character: $c (entry $entry of FASTA input)")
        isspace(ch) && continue
        ch == '>' && error("character '>' not allowed in sequence data (entry $entry of FASTA input)")
        write(io, ch)
        i += 1
        entry_chars += 1
    end
    nl && write(io, '\n')
    return entry_chars
end

"""
    writefasta([io::IO = STDOUT], data)

This version of the function writes to an already opened `IO` stream, defaulting to `STDOUT`.
"""
function writefasta(io::IO, data)
    entry = 0
    for (desc, seq) in data
        entry += 1
        desc = strip(String(desc))
        isascii(desc) || error("description must be ASCII (entry $entry of FASTA input)")
        isempty(desc) && error("empty description (entry $entry of FASTA input")
        if findfirst(==('\n'), desc) ≢ nothing
            error("newlines are not allowed within description (entry $entry of FASTA input)")
        end
        if length(desc) > 79
            @warn("description line longer than 80 characters (entry $entry of FASTA input)")
        end
        println(io, ">", desc)
        entry_chars = writefastaseq(io, seq, entry)
        entry_chars > 0 || error("empty sequence data (entry $entry of FASTA input)")
    end
end
writefasta(data) = writefasta(STDOUT, data)

"""
    writefasta(filename::String, data, [mode::String = "w"])

This function dumps data to a FASTA file, auto-formatting it so to follow the specifications detailed in
the section titled [The FASTA format](@ref). The `data` can be anything which is iterable and which produces
`(description, sequence)` tuples upon iteration, where the `description` must be convertible to
a `String` and the `sequence` can be any iterable object which yields elements convertible
to ASCII characters (e.g. a `String`, a `Vector{UInt8}` etc.).

Examples:

```julia
writefasta("somefile.fasta", [("GENE1", "GCATT"), ("GENE2", "ATTAGC")])
writefasta("somefile.fasta", ["GENE1" => "GCATT", "GENE2" => "ATTAGC"])
```

If the `filename` ends with `.gz`, the result will be a gzip-compressed file.

The `mode` flag determines how the `filename` is open; use `"a"` to append the data to an existing
file.
"""
function writefasta(filename::AbstractString, data, mode::AbstractString = "w")
    fopen = endswith(filename, ".gz") ? gzopen : open
    fopen(filename, mode) do f
        writefasta(f, data)
    end
end

end
