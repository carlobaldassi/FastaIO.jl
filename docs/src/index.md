# FastaIO.jl — FASTA file reader and writer module

```@meta
CurrentModule = FastaIO
```

This module provides ways to parse and write files in
[FASTA format](http://en.wikipedia.org/wiki/FASTA_format) in [Julia](http://julialang.org).
It is designed to be lightweight and fast; the parsing method is inspired by
[kseq.h](http://lh3lh3.users.sourceforge.net/kseq.shtml). It can read and write
files on the fly, keeping only one entry at a time in memory, and it can read and
write gzip-compressed files.

Here is a quick example for reading a file:

```text
julia> using FastaIO

julia> FastaReader("somefile.fasta") do fr
           for (desc, seq) in fr
               println("$desc : $seq")
           end
       end
```

And for writing:

```text
julia> using FastaIO

julia> FastaWriter("somefile.fasta") do fw
           for s in [">GENE1", "GCATT", ">GENE2", "ATTAGC"]
               write(fw, s)
           end
       end
```

## Installation and usage


To install the module, use Julia's package manager: start pkg mode by pressing `]` and then enter:

```
(v1.3) pkg> add FastaIO
```

Dependencies will be installed automatically.
The module can then be loaded like any other Julia module:

```
julia> using FastaIO
```

## Introductory notes

For both reading and writing, there are quick methods to read/write all the data at once: [`readfasta`](@ref) and
[`writefasta`](@ref). These, however, require all the data to be stored in memory at once, which may be impossible
or undesirable for very large files. Therefore, for both reading and writing, the preferred way is actually to
use specialized types, [`FastaReader`](@ref) and [`FastaWriter`](@ref), which have the ability to process one entry
(description + sequence data) at a time (the writer can actually process one char at a time); however, note
that these two object types are not symmetric: the reader acts as an iterable object, while the writer behaves
similarly to an `IO` stream.

## The FASTA format

The FASTA format which is assumed by this module is as follows:

1. description lines must start with a `>` character, and cannot be empty
2. only one description line per entry is allowed
3. all characters must be ASCII
4. whitespace is not allowed within sequence data (except for newlines) and
   at the beginning or end of the description
5. Empty lines are ignored (note however that lines containing whitespace will still trigger an error)

When writing, description lines longer than 80 characters will trigger a warning message; sequence data is
formatted in lines of 80 characters each; extra whitespace is silently discarded.
No other restriction is put on the content of the sequence data, except that the `>` character is
forbidden.

When reading, almost no explicit checks are performed to test that the data actually conforms to these
specifications.

## The sequence storage type

When reading FASTA files, the container type used to store the sequence data can be chosen (as an optional
argument to [`readfasta`](@ref) or as a parametric type of [`FastaReader`](@ref)). The default is
`String`, which is the most memory-efficient and the fastest; another performance-optimal option is
`Vector{UInt8}`, which is a less friendly representation, but has the advantage of being mutable. Any
other container `T` for which `convert(::Type{T}, ::Vector{UInt8})` is defined can be used (e.g.
`Vector{Char}`, or a more specialized `Vector{AminoAcid}` if you use
[BioSeq](https://github.com/diegozea/BioSeq.jl), but the conversion will generally slightly
reduce the performance.

## Reading files

```@docs
readfasta
```

```@docs
FastaReader(::Union{AbstractString,IO})
```

```@docs
FastaReader(::Function, ::AbstractString)
```

```@docs
readentry
```

```@docs
rewind(::FastaReader)
```

```@docs
eof(::FastaReader)
```

```@docs
close(::FastaReader)
```

## Writing files

```@docs
writefasta(::AbstractString, ::Any)
```

```@docs
writefasta(::IO, ::Any)
```

```@docs
FastaWriter
```

```@docs
writeentry
```

```@docs
write(::FastaWriter, ::Any)
```

```@docs
close(::FastaWriter)
```
