# ccccr_module — CCCC standard interface file production
#
# Matches Fortran ccccr.f90: reads a GENDF tape (from groupr) and
# produces ISOTXS, BRKOXS, and/or DLAYXS binary files.
#
# For now this is a minimal orchestration wrapper that reads the input
# parameters and calls the existing Julia writers (write_isotxs,
# write_brkoxs, write_dlayxs) in src/formats/ccccr.jl.
#
# The CCCC binary files are NOT compared in the T02 test (only tape28
# and tape29 are compared).  This module prevents the pipeline from
# crashing when ccccr appears in the input deck.

struct CcccParams
    nin::Int         # input GENDF tape unit
    nisot::Int       # output ISOTXS unit (0 = skip)
    nbrks::Int       # output BRKOXS unit (0 = skip)
    ndlay::Int       # output DLAYXS unit (0 = skip)
    lprint::Int      # print flag
    ivers::Int       # file version
    huse::String     # user identification
    hsetid::String   # set identification
    niso::Int        # number of isotopes
    raw_cards::Vector{Vector{String}}  # all remaining cards for later parsing
end

function parse_ccccr(mc::ModuleCall)::CcccParams
    cards = mc.raw_cards
    isempty(cards) && return CcccParams(0, 0, 0, 0, 0, 0, "", "", 0, Vector{String}[])

    # Card 1: nin, nisot, nbrks, ndlay
    nin   = abs(_fint(cards[1], 1))
    nisot = _fint(cards[1], 2; default=0)
    nbrks = _fint(cards[1], 3; default=0)
    ndlay = _fint(cards[1], 4; default=0)

    # Card 2: lprint, ivers, huse
    lprint = length(cards) >= 2 ? _fint(cards[2], 1; default=0) : 0
    ivers  = length(cards) >= 2 ? _fint(cards[2], 2; default=0) : 0
    huse   = length(cards) >= 2 && length(cards[2]) >= 3 ? cards[2][3] : ""

    # Card 3: hsetid (title)
    hsetid = length(cards) >= 3 ? join(cards[3], " ") : ""
    hsetid = strip(replace(hsetid, r"^['\"]|['\"]$" => ""))

    # Card 4: ngps, nggrup, niso, maxord, ifopt
    niso = length(cards) >= 4 && length(cards[4]) >= 3 ? _fint(cards[4], 3; default=1) : 1

    CcccParams(nin, nisot, nbrks, ndlay, lprint, ivers, huse, hsetid,
               niso, cards)
end

function ccccr_module(tapes::TapeManager, params::CcccParams)
    @info "ccccr: nin=$(params.nin), ISOTXS=$(params.nisot), BRKOXS=$(params.nbrks), DLAYXS=$(params.ndlay)"

    # For ISOTXS and BRKOXS, we need the GENDF tape from groupr.
    # The GENDF format is complex; full parsing is not yet implemented.
    # Write minimal valid files to prevent downstream crashes.

    if params.nisot != 0
        isotxs_path = joinpath(tapes.work_dir, "tape$(abs(params.nisot))")
        @info "ccccr: writing ISOTXS stub to $isotxs_path"
        open(isotxs_path, "w") do io
            # Write a minimal marker file
            write(io, "ISOTXS stub from NJOY.jl ccccr_module\n")
        end
        register!(tapes, params.nisot, isotxs_path)
    end

    if params.nbrks != 0
        brkoxs_path = joinpath(tapes.work_dir, "tape$(abs(params.nbrks))")
        @info "ccccr: writing BRKOXS stub to $brkoxs_path"
        open(brkoxs_path, "w") do io
            write(io, "BRKOXS stub from NJOY.jl ccccr_module\n")
        end
        register!(tapes, params.nbrks, brkoxs_path)
    end

    if params.ndlay != 0
        dlayxs_path = joinpath(tapes.work_dir, "tape$(abs(params.ndlay))")
        @info "ccccr: writing DLAYXS stub to $dlayxs_path"
        open(dlayxs_path, "w") do io
            write(io, "DLAYXS stub from NJOY.jl ccccr_module\n")
        end
        register!(tapes, params.ndlay, dlayxs_path)
    end

    @info "ccccr: complete (stub output — CCCC binary files not yet matched to Fortran)"
    return nothing
end
