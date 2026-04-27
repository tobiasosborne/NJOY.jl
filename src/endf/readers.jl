# readers.jl -- Additional ENDF section readers for module orchestration
#
# These readers extract data needed by heatr and thermr that isn't
# available through the existing reconr/MF2 readers.

# =========================================================================
# MF1/MT451: Material descriptive header (for acer ZAID / incident particle)
# =========================================================================

"""
    read_mf1_mt451_header(filename, mat) -> NamedTuple

Parse the first three records of MF1/MT451 and return:
  (za, awr, lrp, lfi, nlib, nmod, elis, sta, lis, liso, nfor,
   awi, emax, lrel, nsub, nver)

Ref: ENDF-6 manual section 1.1 + njoy-reference/src/endf.f90 mt451 read;
also acefc.f90:280-360 where izai is derived from NSUB/10.

NSUB encodes the incident particle:
  10     → incident neutrons
  0      → photon
  10010  → protons                (IZA = 1001)
  10020  → deuterons              (IZA = 1002)
  10030  → tritons                (IZA = 1003)
  20030  → He-3                   (IZA = 2003)
  20040  → alphas                 (IZA = 2004)
  113    → electro-nuclear
"""
function read_mf1_mt451_header(filename::AbstractString, mat::Integer)
    open(filename, "r") do io
        found = find_section(io, 1, 451; target_mat=Int(mat))
        found || error("MF1/MT451 not found for MAT=$mat in $filename")
        h1 = read_cont(io)
        h2 = read_cont(io)
        h3 = read_cont(io)
        return (
            za   = Int(round(Float64(h1.C1))),
            awr  = Float64(h1.C2),
            lrp  = Int(h1.L1),
            lfi  = Int(h1.L2),
            nlib = Int(h1.N1),
            nmod = Int(h1.N2),
            elis = Float64(h2.C1),
            sta  = Float64(h2.C2),
            lis  = Int(h2.L1),
            liso = Int(h2.L2),
            nfor = Int(h2.N2),
            awi  = Float64(h3.C1),
            emax = Float64(h3.C2),
            lrel = Int(h3.L1),
            nsub = Int(h3.N1),
            nver = Int(h3.N2),
        )
    end
end

"""
    read_mf6_incident_energies(filename, mat, mt) -> Vector{Float64}

Extract the incident-neutron/particle energy grid at which MF6/MT angular
(+ energy) distributions are tabulated. Returns eV values from the first
subsection's TAB2 record — all subsections share the same incident grid
for standard evaluations (NJOY's `topfil` verifies this).

Used by `acer` to union MF6 energies into the ESZ grid (acefc.f90
unionx adds these so each MF6 tabulation point lands on a grid point).
Returns empty vector if MF6/MT is absent.

Ref: ENDF-6 manual §6, plus njoy-reference/src/acefc.f90:1538-1702
(unionx loop over MF6 grids).
"""
function read_mf6_incident_energies(filename::AbstractString, mat::Integer, mt::Integer)
    energies = Float64[]
    open(filename, "r") do io
        found = find_section(io, 6, mt; target_mat=Int(mat))
        found || return energies
        # HEAD: ZA, AWR, JP, LCT, NK, 0
        head = read_cont(io)
        nk = Int(head.N1)
        nk <= 0 && return energies
        # First subsection: yield TAB1 whose HEAD carries LAW in L2 (ENDF-6 §6).
        # read_tab1 consumes the whole TAB1 (CONT + interp + data). After the
        # TAB1, what follows depends on LAW.
        yield = read_tab1(io)
        law = Int(yield.L2)
        if law in (1, 5, 6, 7)
            # TAB2: NR=yield.L1? no, TAB2 has its own interp range. NZ field is
            # the number of incident energies (NE) in the energy distribution.
            tab2 = read_tab2(io)
            ne = Int(tab2.NZ)
            # Each of NE records starts with a CONT header whose C2 = incident E.
            for _ in 1:ne
                rec = read_cont(io)
                push!(energies, Float64(rec.C2))
                nw = Int(rec.N1)
                nlines = cld(nw, 6)
                for _ in 1:nlines; readline(io); end
            end
        elseif law == 2
            # LAW=2 embeds incident energies in a different structure — port later.
        end
    end
    energies
end

"""
    MF6Law5Subsection

One incident-energy block of an MF6/MT LAW=5 (charged-particle elastic)
distribution. ENDF-6 §6.2.5: LIST record with header
(SPI, E, LTP, LIDP, NW, NL), followed by NW=2*NL real numbers in
(μ_k, p_k) pairs.

For LTP=12 (the form T50 uses), `data` holds tabulated cosines and
probabilities directly. For LTP<12, `data` holds Legendre amplitude
coefficients that must be expanded via `ptlegc` to (μ, p) pairs.
"""
struct MF6Law5Subsection
    spi::Float64
    e::Float64       # incident energy [eV]
    ltp::Int
    lidp::Int
    nl::Int          # number of (μ,p) pairs (or Legendre order, depending on LTP)
    mu::Vector{Float64}    # length nl
    prob::Vector{Float64}  # length nl
end

"""
    read_mf6_law5(filename, mat, mt) -> (header::NamedTuple, subs::Vector{MF6Law5Subsection})

Read every NE-incident-energy LIST record from MF6/MT LAW=5 (identical-
particle Coulomb+nuclear elastic). Header carries top-of-section CONT
fields (ZA, AWR, JP, LCT, NK) and the yield TAB1's LAW; subs is one
entry per incident energy.

Returns `(nothing, [])` if MF6/MT or LAW=5 is absent.

Ref: ENDF-6 §6.2.5; njoy-reference/src/acefc.f90:6529-6539 (acecpe LIST
read), 6566-6625 (μ,p extraction).
"""
function read_mf6_law5(filename::AbstractString, mat::Integer, mt::Integer)
    subs = MF6Law5Subsection[]
    header = (za=0.0, awr=0.0, jp=0, lct=0, nk=0, law=0, ne=0,
              spi=0.0, lidp=0)
    open(filename, "r") do io
        find_section(io, 6, mt; target_mat=Int(mat)) || return
        h = read_cont(io)
        za, awr, jp, lct, nk = Float64(h.C1), Float64(h.C2), Int(h.L1), Int(h.L2), Int(h.N1)
        nk <= 0 && return
        yield = read_tab1(io)
        law = Int(yield.L2)
        law == 5 || return  # only LAW=5 here
        # TAB2 for LAW=5 carries the section-level SPI (C1) and LIDP (L1).
        # Ref: njoy-reference/src/acefc.f90:5830-5831 and ENDF-6 §6.2.5.
        tab2 = read_tab2(io)
        ne   = Int(tab2.NZ)
        spi  = Float64(tab2.C1)
        lidp = Int(tab2.L1)
        header = (; za, awr, jp, lct, nk, law, ne, spi, lidp)
        for _ in 1:ne
            rec = read_list(io)
            e   = Float64(rec.C2)
            ltp = Int(rec.L1)
            # rec.N1 = NW, rec.N2 = NL. NL pairs of (μ,p) starting at data[1].
            nl  = Int(rec.N2)
            mu   = Vector{Float64}(undef, nl)
            prob = Vector{Float64}(undef, nl)
            for k in 1:nl
                mu[k]   = rec.data[2k-1]
                prob[k] = rec.data[2k]
            end
            # spi and lidp inherited from the TAB2 (per Fortran acefc.f90:5830).
            push!(subs, MF6Law5Subsection(spi, e, ltp, lidp, nl, mu, prob))
        end
    end
    header, subs
end

"""
    acer_incident_letter(nsub) -> Char

Map ENDF-6 NSUB (MF1/MT451 third CONT N1 field) to the single-letter ACE
suffix tag for the incident particle. Matches Fortran acer.f90:391-407.
Unknown NSUB → 'c' (neutron) with a warning.
"""
function acer_incident_letter(nsub::Integer)
    nsub == 10    ? 'c' :
    nsub == 0     ? 'p' :  # photons (photoatomic)
    nsub == 3     ? 'e' :  # electrons — acer uses 'u' for electroatomic in practice; keep 'e'
    nsub == 113   ? 'n' :  # electro-nuclear
    nsub == 10010 ? 'h' :  # protons
    nsub == 10020 ? 'o' :  # deuterons — NJOY MCNPX mapping uses 'o' for deuterons
    nsub == 10030 ? 'r' :  # tritons
    nsub == 20030 ? 's' :  # He-3
    nsub == 20040 ? 'a' :  # alphas
    (@warn("acer: unknown NSUB=$nsub — defaulting to 'c'"); 'c')
end

# =========================================================================
# MF12: Photon multiplicity data (for heatr gamma recoil)
# =========================================================================

"""
    read_mf12_gammas(filename::AbstractString, mat::Integer; mt::Integer=102)
        -> Vector{Tuple{Float64, Float64}}

Read MF12 photon yield data for a given material and MT. Returns a vector
of (E_gamma, yield) tuples, one per discrete gamma line.

The first TAB1 subsection (total yield) is skipped; per-gamma subsections
follow, each starting with E_gamma in the TAB1 C1 field.

Replaces the manual 60-line parsing block in t01_pipeline.jl.
"""
function read_mf12_gammas(filename::AbstractString, mat::Integer; mt::Integer=102)
    gammas = Tuple{Float64, Float64}[]

    open(filename, "r") do io
        found = find_section(io, 12, mt; target_mat=Int(mat))
        found || return gammas

        # HEAD record: ZA, AWR, LO, 0, NK, 0
        head = read_cont(io)
        nk = Int(head.N1)
        nk <= 0 && return gammas

        # First subsection: total yield TAB1 — skip it
        _discard_tab1(io)

        # NK per-gamma subsections follow the total
        for k in 1:nk
            tab1 = read_tab1(io)
            eg = Float64(tab1.C1)              # gamma energy [eV]
            yield_val = isempty(tab1.y) ? 0.0 : Float64(tab1.y[1])  # yield at first energy
            push!(gammas, (eg, yield_val))
        end
    end

    gammas
end

"""Skip one TAB1 record (interp params + data pairs) from an IO stream."""
function _discard_tab1(io::IO)
    read_tab1(io)  # just read and discard
    nothing
end

# =========================================================================
# MF5: Energy distribution data (for heatr evaporation)
# =========================================================================

"""
    read_mf5_evaporation(filename::AbstractString, mat::Integer; mt::Integer=91)
        -> Union{Nothing, NamedTuple{(:u, :theta, :lf), Tuple{Float64, Float64, Int}}}

Read MF5 evaporation spectrum parameters for a given MT.
Returns (u=threshold, theta=nuclear_temperature, lf=distribution_law)
or `nothing` if MF5/MT not found.

Only handles LF=9 (simple evaporation spectrum) and LF=5 (general evaporation).
For LF=9, theta is the constant nuclear temperature.
"""
function read_mf5_evaporation(filename::AbstractString, mat::Integer; mt::Integer=91)
    open(filename, "r") do io
        found = find_section(io, 5, mt; target_mat=Int(mat))
        found || return nothing

        # HEAD: ZA, AWR, 0, 0, NK, 0
        head = read_cont(io)
        nk = Int(head.N1)
        nk <= 0 && return nothing

        # First subsection: probability TAB1 with LF in L2
        prob_tab1 = read_tab1(io)
        u = Float64(prob_tab1.C1)     # threshold energy [eV]
        lf = Int(prob_tab1.L2)        # distribution law

        if lf == 9  # simple fission / evaporation spectrum
            # Next TAB1: theta(E)
            theta_tab1 = read_tab1(io)
            theta = isempty(theta_tab1.y) ? 0.0 : Float64(theta_tab1.y[1])
            return (u=u, theta=theta, lf=lf)
        elseif lf == 5  # general evaporation
            # Next TAB1: theta(E), same format
            theta_tab1 = read_tab1(io)
            theta = isempty(theta_tab1.y) ? 0.0 : Float64(theta_tab1.y[1])
            return (u=u, theta=theta, lf=lf)
        end

        # Other LF values: return what we have
        return (u=u, theta=0.0, lf=lf)
    end
end
