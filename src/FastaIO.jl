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

import Base.start, Base.done, Base.next, Base.readall,
       Base.close, Base.show, Base.eof, Base.write

const fasta_buffer_size = 4096

type FastaReader{T}
    # public but read-only
    num_parsed::Int        # number of parsed entries so far
    # private
    f::IO
    is_eof::Bool           # did we reach end of file?
    rbuffer::Vector{Uint8} # read buffer
    rbuf_sz::Int           # read buffer size
    rbuf_pos::Int          # read buffer cursor
    lbuffer::Vector{Uint8} # line buffer
    lbuf_sz::Int           # line buffer size
    mbuffer::Vector{Uint8} # multi-line buffer
    mbuf_sz::Int           # multi-line buffer size
    own_f::Bool
    function FastaReader(filename::String)
        fr = new(0, gzopen(filename), false,
            Array(Uint8, fasta_buffer_size), 0, 0,
            Array(Uint8, fasta_buffer_size), 0,
            Array(Uint8, fasta_buffer_size), 0,
            true)
        finalizer(fr, close)
        return fr
    end
    function FastaReader(io::IO)
        new(0, io, false,
            Array(Uint8, fasta_buffer_size), 0, 0,
            Array(Uint8, fasta_buffer_size), 0,
            Array(Uint8, fasta_buffer_size), 0,
            false)
    end
end

FastaReader(filename::String) = FastaReader{ASCIIString}(filename)
FastaReader(io::IO) = FastaReader{ASCIIString}(io)

function FastaReader(f::Function, filename::String, T::Type=ASCIIString)
    fr = FastaReader{T}(filename)
    try
        f(fr)
    finally
        close(fr)
    end
end

close(fr::FastaReader) = fr.own_f && close(fr.f)

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
    ccall(:ios_readall, Uint, (Ptr{Void}, Ptr{Void}, Uint), fr.f.ios, fr.rbuffer, fasta_buffer_size)
function read_chunk_ll(fr::FastaReader, s::IO)
    ret = 0
    while !eof(fr.f) && ret < fasta_buffer_size
        ret += 1
        fr.rbuffer[ret] = read(fr.f, Uint8)
    end
    return ret
end

function read_chunk(fr::FastaReader)
    if fr.is_eof
        return
    end
    ret = read_chunk_ll(fr, fr.f)
    ret == -1 && error("read failure")
    fr.rbuf_sz = ret
    fr.rbuf_pos = 1
    if ret == 0
        fr.is_eof = true
    end
    return
end

function readline(fr::FastaReader)
    fr.lbuf_sz = 0
    found = false
    while !fr.is_eof
        if fr.rbuf_pos == 0
            read_chunk(fr::FastaReader)
        end
        i = fr.rbuf_pos
        cr = false
        while i <= fr.rbuf_sz
            c = fr.rbuffer[i]
            if c == '\n'
                found = true
                break
            else
                cr = (c == '\r')
            end
            i += 1
        end
        i -= 1 + cr
        chunk_len = i - fr.rbuf_pos + 1
        free_sbuf = length(fr.lbuffer) - fr.lbuf_sz
        gap = chunk_len - free_sbuf
        if gap > 0
            resize!(fr.lbuffer, length(fr.lbuffer) + gap)
        end

        #fr.lbuffer[fr.lbuf_sz + (1:chunk_len)] = fr.rbuffer[fr.rbuf_pos:i]
        copy!(fr.lbuffer, fr.lbuf_sz + 1, fr.rbuffer, fr.rbuf_pos, chunk_len)
        fr.lbuf_sz += chunk_len

        i += 2 + cr
        if i > fr.rbuf_sz
            i = 0
        end
        fr.rbuf_pos = i
        if found
            break
        end
    end
    return
end

function start(fr::FastaReader)
    rewind(fr)
    readline(fr)
    if fr.lbuf_sz == 0
        error("empty FASTA file")
    end
    return
end
done(fr::FastaReader, x::Nothing) = fr.is_eof
function _next_step(fr::FastaReader)
    if fr.lbuffer[1] != '>'
        error("invalid FASTA file: description does not start with '>'")
    end
    if fr.lbuf_sz == 1
        error("invalid FASTA file: empty description")
    end
    name = ascii(fr.lbuffer[2:fr.lbuf_sz])
    fr.mbuf_sz = 0
    while true
        readline(fr)
        if fr.is_eof || fr.lbuffer[1] == '>'
            break
        end
        gap = fr.lbuf_sz - (length(fr.mbuffer) - fr.mbuf_sz)
        if gap > 0
            resize!(fr.mbuffer, length(fr.mbuffer) + gap)
        end
        #fr.mbuffer[fr.mbuf_sz + (1:fr.lbuf_sz)] = fr.lbuffer[1:fr.lbuf_sz]
        copy!(fr.mbuffer, fr.mbuf_sz + 1, fr.lbuffer, 1, fr.lbuf_sz)
        fr.mbuf_sz += fr.lbuf_sz
    end
    return name
end
function _next(fr::FastaReader{Vector{Uint8}})
    name = _next_step(fr)
    fr.num_parsed += 1
    return (name, fr.mbuffer[1:fr.mbuf_sz])
end
function _next(fr::FastaReader{ASCIIString})
    name = _next_step(fr)
    out_str = ccall(:jl_pchar_to_string, ByteString, (Ptr{Uint8},Int), fr.mbuffer, fr.mbuf_sz)
    fr.num_parsed += 1
    return (name, out_str)
end
function _next{T}(fr::FastaReader{T})
    name = _next_step(fr)
    fr.num_parsed += 1
    return (name, convert(T, fr.mbuffer[1:fr.mbuf_sz]))
end

next(fr::FastaReader, x::Nothing) = (_next(fr), nothing)

function readall(fr::FastaReader)
    ret = Any[]
    for item in fr
        push!(ret, item)
    end
    return ret
end

function readentry(fr::FastaReader)
    fr.is_eof && throw(EOFError())
    if fr.num_parsed == 0
        readline(fr)
        if fr.lbuf_sz == 0
            error("empty FASTA file")
        end
    end
    item, _ = next(fr, nothing)
    return item
end

eof(fr::FastaReader) = fr.is_eof

function show{T}(io::IO, fr::FastaReader{T})
    print(io, "FastaReader(input=\"$(fr.f)\", out_type=$T, num_parsed=$(fr.num_parsed), eof=$(fr.is_eof))")
end

function readfasta(filename::String, T::Type=ASCIIString)
    FastaReader(filename, T) do fr
        readall(fr)
    end
end
readfasta(io::IO, T::Type=ASCIIString) = readall(FastaReader{T}(io))

type FastaWriter
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
        finalizer(fw, close)
        return fw
    end
    function FastaWriter(filename::String, mode::String = "w")
        if endswith(filename, ".gz")
            of = gzopen
        else
            of = open
        end
        fw = new(of(filename, mode), false, 0, 0, false, 0, 1, true, true)
        finalizer(fw, close)
        return fw
    end
end

FastaWriter() = FastaWriter(STDOUT)

function FastaWriter(f::Function, args...)
    fw = FastaWriter(args...)
    try
        f(fw)
    finally
        close(fw)
    end
end

function write(fw::FastaWriter, c)
    ch = char(c)
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
            warn("description line longer than 80 characters (entry $(fw.entry) of FASTA input)")
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

function write(fw::FastaWriter, s::String)
    for c in s
        write(fw, c)
    end
    write(fw, '\n')
end

function writeentry(fw::FastaWriter, desc::String, seq)
    !fw.at_start && write(fw, '\n')
    desc = strip(ascii(desc))
    if search(desc, '\n') != 0
        error("newlines are not allowed within description (entry $(fw.entry+1) of FASTA input)")
    end
    write(fw, '>')
    write(fw, strip(desc))
    write(fw, '\n')
    #write(fw, seq)
    #write(fw, '\n')
    fw.entry_chars = writefastaseq(fw.f, seq, fw.entry, false)
    fw.in_seq = true
    fw.parsed_nl = false
    fw.pos = 0
    fw.entry_chars > 0 || error("empty sequence data (entry $(fw.entry) of FASTA input)")
    return
end

function close(fw::FastaWriter)
    try
        write(fw.f, '\n')
        flush(fw.f)
    catch err
        isa(err, EOFError) || rethrow(err)
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
        ch = char(c)
        isascii(ch) || error("invalid (non-ASCII) character: $c (entry $entry of FASTA input)")
        isspace(ch) && continue
        ch != '>' || error("character '>' not allowed in sequence data (entry $entry of FASTA input)")
        write(io, ch)
        i += 1
        entry_chars += 1
    end
    nl && write(io, '\n')
    return entry_chars
end

function writefasta(io::IO, data)
    entry = 0
    for (desc, seq) in data
        entry += 1
        desc = strip(ascii(desc))
        if isempty(desc)
            error("empty description (entry $entry of FASTA input")
        end
        if search(desc, '\n') != 0
            error("newlines are not allowed within description (entry $entry of FASTA input)")
        end
        if length(desc) > 79
            warn("description line longer than 80 characters (entry $entry of FASTA input)")
        end
        println(io, ">", desc)
        entry_chars = writefastaseq(io, seq, entry)
        entry_chars > 0 || error("empty sequence data (entry $entry of FASTA input)")
    end
end
writefasta(data) = writefasta(STDOUT, data)

function writefasta(filename::String, data, mode::String = "w")
    if endswith(filename, ".gz")
        of = gzopen
    else
        of = open
    end
    of(filename, mode) do f
        writefasta(f, data)
    end
end

end
