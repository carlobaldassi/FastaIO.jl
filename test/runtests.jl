module FastaTests

using FastaIO
using GZip
using Test

const fastadata_ascii = Any[
    ("A0ADS9_STRAM/3-104",
     "-------------------------------------------STVELTKEN-F--D-Q-" *
     "-T-V----T-D---------N-----------E--F-V--LI--------D-----F--W" *
     "----A----S-W---------C---G------P----C-----R----Q---F------A" *
     "------P---V----Y-------E---K---A---A--E---A---NP------------" *
     "----------------DL---V--F--G----K-------V-------D-----------" *
     "T-----------E--------A------Q-----------------P------E------" *
     "------L----A--------Q-----A------F------G--------I------Q---" *
     "-----S------I------P-T---L--M------I---V-------R---D------Q-" *
     "---V----A-VF-----------AQP-GA----LP-EE-ALTDVI---GQA---------" *
     "------------------------------------------------------------" *
     "---------"),
    ("A0AHX4_LISW6/1-102",
     "-------------------------------------------MVKEITDAT-F--E-Q-" *
     "-E-T----S-E---------------------G--L-V--LT--------D-----F--W" *
     "----A----T-W---------C---G------P----C-----R----M---V------A" *
     "------P---V----L-------E---E---I---Q--E---E---RGE-----------" *
     "----------------AL---K--I--V----K-------M-------D-----------" *
     "V-----------D--------E------N-----------------P------E------" *
     "------T----P--------G-----S------F------G--------V------M---" *
     "-----S------I------P-T---L--L------I---K-------K---D------G-" *
     "---E----V-VE-----------TII-GY----HP-KE-ELDEVINK-Y-----------" *
     "------------------------------------------------------------" *
     "---------"),
    ("A0AJ61_LISW6/3-103",
     "-----------------------------------NLESVEQFD----------------" *
     "---V----I-K-------S-E-----------G--K-S--VF--------M-----F--S" *
     "----A----D-W---------C---G------D----C-----K----Y---I------E" *
     "------P---V----M-------P---E---I---E--A---E---NE------------" *
     "----------------DF---T--F--Y----H-------V-------D-----------" *
     "R-----------D--------E------F-----------------I------D------" *
     "------L----C--------A-----D------L------A--------I------F---" *
     "-----G------I------P-S---F--L------V---F-------E---D------G-" *
     "---E----E-VG-----------RFV-SKD--RKT-KE-EINDFLA--AI----------" *
     "------------------------------------------------------------" *
     "---------"),
    ("A0AK08_LISW6/41-151",
     "-----------------------------------------FLN--TISTKD-F--K-Q-" *
     "-Q-M----A-D---------K-----------T--T-G--FV--------Y-----V--G" *
     "----R----P-T---------C---E------D----C-----Q----A---F------Q" *
     "------P---I----L-------K---K---E---L--K---E---RKLN----------" *
     "---------------QNM---N--Y--Y----N-------T-------D-----------" *
     "K-----------A--------S------E-----------------K------SRDD---" *
     "---MIAL----L--------K-----K------M------D--------I------D---" *
     "-----S------V------P-T---M--V------Y---L-------K---D------G-" *
     "---K----V-AS-----------TYA-AT----DE-PE-KLTHWMNK-V-----------" *
     "------------------------------------------------------------" *
     "---------")]

const fastadata_uint8 = map(x->(x[1],Vector{UInt8}(codeunits(x[2]))), fastadata_ascii)
const fastadata_char = map(x->(x[1],Vector{Char}(x[2])), fastadata_ascii)

function test_fastaread(T::Type, infile, fastadata)
    FastaReader(infile, T) do fr
        @test readentries(fr) == fastadata

        rewind(fr)

        for (desc, seq) in fr
            @test desc == fastadata[fr.num_parsed][1]
            @test seq == fastadata[fr.num_parsed][2]
        end

        rewind(fr)
        while !eof(fr)
            desc, seq = readentry(fr)
            @test desc == fastadata[fr.num_parsed][1]
            @test seq == fastadata[fr.num_parsed][2]
        end
    end
end

function test_fastawrite(infile, outfile, fastadata)
    # char-by-char mode
    FastaWriter(outfile) do fw
        gzopen(infile) do f
            while !eof(f)
                write(fw, read(f, UInt8))
            end
        end
    end
    @test readfasta(outfile) == fastadata

    # writeentry mode
    FastaWriter(outfile) do fw
        for (desc, seq) in fastadata
            writeentry(fw, desc, seq)
        end
    end
    @test readfasta(outfile) == fastadata

    # one element at a time mode
    FastaWriter(outfile) do fw
        fd_plain = Any[]
        for (desc, seq) in fastadata
            push!(fd_plain, ">" * desc)
            i = 0
            while i < length(seq)
                j = rand(i+1:length(seq))
                push!(fd_plain, seq[i+1:j])
                i = j
            end
            #push!(fd_plain, seq)
        end
        write(fw, fd_plain)
    end
    @test readfasta(outfile) == fastadata

    # writefasta (vector)
    writefasta(outfile, fastadata)
    @test readfasta(outfile) == fastadata

    # writefasta (dict)
    ## how to do this in a backward-compatible way?
    ##   fd_dict = Dict(desc=>seq for (desc, seq) in fastadata)
    fd_dict = Dict(); for (desc, seq) in fastadata; fd_dict[desc]=seq; end
    writefasta(outfile, fd_dict)
    @test sort!(readfasta(outfile)) == sort(fastadata)
    return
end

for suffix in ["", ".gz", ".win", ".win.gz", ".no_eof", ".no_eof.gz", ".blanklines", ".blanklines.gz"]
    infile = joinpath(dirname(Base.source_path()), "test.fasta" * suffix)
    outfile = joinpath(dirname(Base.source_path()), "test_out.fasta" * suffix)

    tests = [(Vector{UInt8}, fastadata_uint8),
             (Vector{Char}, fastadata_char),
             (String, fastadata_ascii)]

    for (T, fastadata) in tests
        test_fastaread(T, infile, fastadata)
    end

    try
        test_fastawrite(infile, outfile, fastadata_ascii)
    finally
        isfile(outfile) && rm(outfile)
    end
end

for i in 1:4
    infile = joinpath(dirname(Base.source_path()), "invalid_test_$i.fasta.gz")

    @test_throws Exception readfasta(infile)
    @test_throws Exception FastaReader(infile) do fr
        while !eof(fr)
            desc, seq = readentry(fr)
        end
    end
end

outfile = joinpath(dirname(Base.source_path()), "invalid_test_out.fasta.gz")

for i in 2:6
    infile = joinpath(dirname(Base.source_path()), "invalid_test_$i.fasta.gz")

    try
        @test_throws ErrorException FastaWriter(outfile) do fw
            gzopen(infile) do f
                while !eof(f)
                    write(fw, read(f, UInt8))
                end
            end
        end
    finally
        isfile(outfile) && rm(outfile)
    end
end

const invalid_fastadata_entries = [("DE\nSC", "DATA" ),
                                   (""      , "DATA" ),
                                   ("DESC"  , ""     ),
                                   ("ΔΕΣΞ"  , "DATA" ),
                                   ("DESC"  , "ΔΑΤΑ" ),
                                   ("DESC"  , "DA>TA")]

for (desc, data) in invalid_fastadata_entries
    try
        @test_throws Exception FastaWriter(outfile) do fw
            writeentry(fw, desc, data)
        end
    finally
        isfile(outfile) && rm(outfile)
    end
end

for invalid_fastadata in [[ife] for ife in invalid_fastadata_entries]
    try
        @test_throws Exception writefasta(outfile, invalid_fastadata)
    finally
        isfile(outfile) && rm(outfile)
    end
end

end
