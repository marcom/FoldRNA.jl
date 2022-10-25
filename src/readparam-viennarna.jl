function readparam_viennarna_to_dict(filename::AbstractString;
                                     INF = 10_000_000,
                                     DEF = -50)
    INT_TABLES = [
        "stack", "stack_enthalpies",
        "mismatch_hairpin", "mismatch_hairpin_enthalpies",
        "mismatch_interior", "mismatch_interior_enthalpies",
        "mismatch_interior_1n", "mismatch_interior_1n_enthalpies",
        "mismatch_interior_23", "mismatch_interior_23_enthalpies",
        "mismatch_multi", "mismatch_multi_enthalpies",
        "mismatch_exterior", "mismatch_exterior_enthalpies",
        "dangle5", "dangle5_enthalpies",
        "dangle3", "dangle3_enthalpies",
        "int11", "int11_enthalpies",
        "int21", "int21_enthalpies",
        "int22", "int22_enthalpies",
        "hairpin", "hairpin_enthalpies",
        "bulge", "bulge_enthalpies",
        "interior", "interior_enthalpies",
    ]
    EXTRA_TABLES = ["NINIO", "ML_params", "Misc"]
    LOOP_TABLES = ["Triloops", "Tetraloops", "Hexaloops"]
    paramdict = Dict{String,Any}()
    lines = readlines(filename)
    headeridx = findall(x -> occursin(r"^#", x), lines)
    if lines[headeridx[1]] != "## RNAfold parameter file v2.0"
        error("couldn't find header line: \"## RNAfold parameter file v2.0\"")
    end
    if lines[headeridx[end]][2:end] |> strip != "END"
        error("couldn't find END marker")
    end
    for i in 2:length(headeridx)-1
        startidx = headeridx[i] + 1
        endidx = headeridx[i+1] - 1
        header = lines[headeridx[i]][2:end] |> strip
        datalines = (lines[startidx:endidx]
                     .|> x -> replace(x, r"/\*.*\*/" => "") # remove comments
                     .|> x -> replace(x, "DEF" => DEF)
                     .|> x -> replace(x, "INF" => INF)
                    ) |> x -> filter(a -> ! isempty(strip(a)), x) # ignore empty lines
        datatxt = join(datalines, "\n")
        if header in INT_TABLES
            paramdict[header] = map(x -> parse(Int, x), split(datatxt))
        elseif header in EXTRA_TABLES
            if header == "NINIO"
                a = map(x -> parse(Int, x), split(datatxt))
                if length(a) != 3
                    error("expected exactly 3 values in NINIO table")
                end
                paramdict[header] = Dict{String, Int}(
                    "m" => a[1], "m_dH" => a[2], "max" => a[3]
                )
            elseif header == "ML_params"
                a = map(x -> parse(Int, x), split(datatxt))
                length(a) == 6 || error("expected exactly 6 values in ML_params table")
                paramdict[header] = Dict{String, Int}(
                    "cu" => a[1], "cu_dH" => a[2],
                    "cc" => a[3], "cc_dH" => a[4],
                    "ci" => a[5], "ci_dH" => a[6]
                )
            elseif header == "Misc"
                a = split(datatxt)
                if length(a) != 4 && length(a) != 6
                    error("expected either 4 or 6 values in Misc table")
                end
                paramdict[header] = Dict{String, Int}(
                    "DuplexInit" => parse(Int, a[1]), "DuplexInit_dH" => parse(Int, a[2]),
                    "TerminalAU" => parse(Int, a[3]), "TerminalAU_dH" => parse(Int, a[4])
                )
                if length(a) == 6
                    paramdict["Misc_LXC"] = Dict{String, Float64}(
                        "LXC" => parse(Float64, a[5]), "LXC_dH" => parse(Float64, a[6])
                    )
                end
            else
                error("unhandled header $header")
            end
        elseif header in LOOP_TABLES
            looptable = Dict{String, Tuple{Int, Int}}()
            for d in datalines
                seq, a, b = split(d)
                looptable[String(seq)] = (parse(Int, a), parse(Int, b))
            end
            paramdict[header] = looptable
        else
            error("unhandled header $header")
        end
    end
    return paramdict
end

"""
    readparam_viennarna(filepath; options...)

Read nearest-neighbour energy parameter file in ViennaRNA RNAfold v2.0
format.

Some implicit defaults for the ViennaRNA parameter file format:
- energy values are in units of 0.01 * kcal/mol
- parameters are assumed to be at 37 degrees Celsius
These can be configured with keyword arguments.

Note: for the intloop22 tables (named `int22` and `int22_enthalpies`
in the ViennaRNA parameter file), the tables are assumed to not
consider bases or basepairs that contain an `'N'`.
"""
function readparam_viennarna(filepath::AbstractString;
                             # ViennaRNA file format defaults
                             NUMTYPE = Int,
                             SEQTYPE = Int,
                             MAXLOOPS = 30,
                             unit = 0.01u"kcal/mol",
                             DEFAULT_LXC = 107.856,
                             DEFAULT_LXC_DH = 0.0,
                             hpmin = 3,
                             RT = RT37,
                             # file_bases, file_basepairs: base and basepair ordering as
                             # used in the ViennaRNA RNAfold v2.0 file format
                             file_bases = ['N', 'A', 'C', 'G', 'U'],
                             file_basepairs = [('C','G'), ('G','C'), ('G','U'), ('U','G'), ('A','U'), ('U','A'), ('N','N')],
                             # file_wildcards: bases that are wildcard bases and not real bases
                             file_wildcards = ['N'],
                             )
    isgcbp(bp::Tuple{Char,Char}) = bp == ('G','C') || bp == ('C','G')
    reshape_rowmajor(A, dims) = permutedims(reshape(A, reverse(dims)), length(dims):-1:1)
    nb = length(file_bases)
    nbp = length(file_basepairs)

    # outputs
    #
    # Note: output base ordering is different, as the intloop22
    # parameters in ViennaRNA don't consider N bases or basepairs
    # containing N (probably to conserve space).
    out_bases = ['A', 'C', 'G', 'U', 'N']
    out_basepairs = [('C','G'), ('G','C'), ('G','U'), ('U','G'), ('A','U'), ('U','A'), ('N','N')]
    alphabet = Alphabet("RNA", out_bases; wildcard_chars = file_wildcards)
    deltaG = LoopModel{NUMTYPE,SEQTYPE,nb,nbp,MAXLOOPS}(
        alphabet = alphabet,
        unit = unit,
        name = basename(filepath),
        hpmin = hpmin,
        RT = RT
    )
    deltaH = LoopModel{NUMTYPE,SEQTYPE,nb,nbp,MAXLOOPS}(
        alphabet = alphabet,
        unit = unit,
        name = basename(filepath),
        hpmin = hpmin,
        RT = RT
    )
    # set bptype
    for (i, (a,b)) in enumerate(out_basepairs)
        ka = first(encode(alphabet, a))
        kb = first(encode(alphabet, b))
        deltaG.bptype[ka, kb] = i
        deltaH.bptype[ka, kb] = i
    end

    p = readparam_viennarna_to_dict(filepath)

    # remapping of indices between our ordering and ViennaRNA's
    # for bases and basepairs, e.g. we might have A = 3 and
    # ViennaRNA might have A = 1
    # IbN, IbpN: the same arrays but excluding base 'N'
    Ib   = indexin(out_bases, file_bases)
    IbN  = indexin(filter(x -> x != 'N', out_bases), filter(x -> x != 'N', file_bases))
    Ibp  = indexin(out_basepairs, file_basepairs)
    IbpN = indexin(filter(x -> x != ('N','N'), out_basepairs), filter(x -> x != ('N','N'), file_basepairs))

    # stack, stack_enthalpies
    deltaG.stack = reshape_rowmajor(p["stack"], (nbp,nbp))[Ibp,Ibp]
    deltaH.stack = reshape_rowmajor(p["stack_enthalpies"], (nbp,nbp))[Ibp,Ibp]

    # hairpin, bulge, intloop
    # TODO: type of fields is hardcoded here
    deltaG.hairpin_init  = OffsetArray(p["hairpin"], 0:MAXLOOPS)
    deltaG.bulge_init    = OffsetArray(p["bulge"], 0:MAXLOOPS)
    deltaG.intloop_init  = OffsetArray(p["interior"], 0:MAXLOOPS)

    # intloop11, intloop12, intloop22
    deltaG.intloop11  = reshape_rowmajor(p["int11"], (nbp,nbp,nb,nb))[Ibp,Ibp,Ib,Ib]
    deltaH.intloop11  = reshape_rowmajor(p["int11_enthalpies"], (nbp,nbp,nb,nb))[Ibp,Ibp,Ib,Ib]
    deltaG.intloop12  = reshape_rowmajor(p["int21"], (nbp,nbp,nb,nb,nb))[Ibp,Ibp,Ib,Ib,Ib]
    deltaH.intloop12  = reshape_rowmajor(p["int21_enthalpies"], (nbp,nbp,nb,nb,nb))[Ibp,Ibp,Ib,Ib,Ib]
    # int22 table doesn't have entries for NN basepairs and N bases
    deltaG.intloop22  = reshape_rowmajor(p["int22"], (nbp-1,nbp-1,nb-1,nb-1,nb-1,nb-1))[IbpN,IbpN,IbN,IbN,IbN,IbN]
    deltaH.intloop22  = reshape_rowmajor(p["int22_enthalpies"], (nbp-1,nbp-1,nb-1,nb-1,nb-1,nb-1))[IbpN,IbpN,IbN,IbN,IbN,IbN]

    # dangle5, dangle3
    deltaG.dangle5 = reshape_rowmajor(p["dangle5"], (nbp,nb))[Ibp,Ib]
    deltaH.dangle5 = reshape_rowmajor(p["dangle5_enthalpies"], (nbp,nb))[Ibp,Ib]
    deltaG.dangle3 = reshape_rowmajor(p["dangle3"], (nbp,nb))[Ibp,Ib]
    deltaH.dangle3 = reshape_rowmajor(p["dangle3_enthalpies"], (nbp,nb))[Ibp,Ib]

    # mismatch_{hairpin,intloop,intloop1n,intloop23,multiloop,extloop}
    deltaG.mismatch_hairpin   = reshape_rowmajor(p["mismatch_hairpin"], (nbp,nb,nb))[Ibp,Ib,Ib]
    deltaH.mismatch_hairpin   = reshape_rowmajor(p["mismatch_hairpin_enthalpies"], (nbp,nb,nb))[Ibp,Ib,Ib]
    deltaG.mismatch_intloop   = reshape_rowmajor(p["mismatch_interior"], (nbp,nb,nb))[Ibp,Ib,Ib]
    deltaH.mismatch_intloop   = reshape_rowmajor(p["mismatch_interior_enthalpies"], (nbp,nb,nb))[Ibp,Ib,Ib]
    deltaG.mismatch_intloop1n = reshape_rowmajor(p["mismatch_interior_1n"], (nbp,nb,nb))[Ibp,Ib,Ib]
    deltaH.mismatch_intloop1n = reshape_rowmajor(p["mismatch_interior_1n_enthalpies"], (nbp,nb,nb))[Ibp,Ib,Ib]
    deltaG.mismatch_intloop23 = reshape_rowmajor(p["mismatch_interior_23"], (nbp,nb,nb))[Ibp,Ib,Ib]
    deltaH.mismatch_intloop23 = reshape_rowmajor(p["mismatch_interior_23_enthalpies"], (nbp,nb,nb))[Ibp,Ib,Ib]
    deltaG.mismatch_multiloop = reshape_rowmajor(p["mismatch_multi"], (nbp,nb,nb))[Ibp,Ib,Ib]
    deltaH.mismatch_multiloop = reshape_rowmajor(p["mismatch_multi_enthalpies"], (nbp,nb,nb))[Ibp,Ib,Ib]
    deltaG.mismatch_extloop   = reshape_rowmajor(p["mismatch_exterior"], (nbp,nb,nb))[Ibp,Ib,Ib]
    deltaH.mismatch_extloop   = reshape_rowmajor(p["mismatch_exterior_enthalpies"], (nbp,nb,nb))[Ibp,Ib,Ib]

    # multiloop_{stem,unpaired,init}
    # TODO: check where ci, cu, cc belong
    deltaG.multiloop_branch     = p["ML_params"]["ci"]
    deltaH.multiloop_branch     = p["ML_params"]["ci_dH"]
    deltaG.multiloop_unpaired = p["ML_params"]["cu"]
    deltaH.multiloop_unpaired = p["ML_params"]["cu_dH"]
    deltaG.multiloop_init     = p["ML_params"]["cc"]
    deltaH.multiloop_init     = p["ML_params"]["cc_dH"]

    # TODO: extloop_unpaired isn't set, check to confirm

    # specialhairpins
    shp = Dict(p["Triloops"]..., p["Tetraloops"]..., p["Hexaloops"]...)
    deltaG.specialhairpins = Dict(encode(alphabet, k) => shp[k][1] for k in keys(shp))
    deltaH.specialhairpins = Dict(encode(alphabet, k) => shp[k][2] for k in keys(shp))

    # ninio_{m,max}
    deltaG.ninio_m = p["NINIO"]["m"]
    deltaH.ninio_m = p["NINIO"]["m_dH"]
    deltaG.ninio_max = p["NINIO"]["max"]
    deltaH.ninio_max = p["NINIO"]["max"]

    # duplex_init, terminal_AU, lxc
    deltaG.duplex_init = p["Misc"]["DuplexInit"]
    deltaH.duplex_init = p["Misc"]["DuplexInit_dH"]
    deltaG.terminal_nonGC = p["Misc"]["TerminalAU"]
    deltaH.terminal_nonGC = p["Misc"]["TerminalAU_dH"]

    deltaG.terminal_nonGC_bp = [isgcbp(b) ? zero(NUMTYPE) : deltaG.terminal_nonGC for b in out_basepairs]
    deltaH.terminal_nonGC_bp = [isgcbp(b) ? zero(NUMTYPE) : deltaH.terminal_nonGC for b in out_basepairs]

    if haskey(p, "Misc_LXC")
        deltaG.lxc = p["Misc_LXC"]["LXC"]
        deltaH.lxc = p["Misc_LXC"]["LXC_dH"]
    else
        deltaG.lxc = DEFAULT_LXC
        deltaH.lxc = DEFAULT_LXC_DH
    end

    return deltaG, deltaH
end
