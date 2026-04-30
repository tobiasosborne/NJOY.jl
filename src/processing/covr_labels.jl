# covr_labels.jl -- Plot-label generation for covr.
#
# Ports the label/name-table helpers from NJOY2016 covr.f90:
#   elem    (covr.f90:1753-1779) -- Z → element name
#   mtno    (covr.f90:1781-1890) -- MT → reaction-name strings
#   matmes  (covr.f90:1705-1751) -- (iza, mat, mt) → "<isotope><reaction>"
#   smilab  (covr.f90:1649-1703) -- semi-log plot title for one (mat,mt) frame
#
# All character-table ordering and padding matches the Fortran source
# verbatim; covr emits these labels into the plot tape, where viewr
# consumes the embedded #H markup to render PostScript text. Any deviation
# breaks the byte-for-byte plot-tape diff.

# Element symbol table — matches covr.f90:1761-1773 verbatim. Index = Z (1..103).
# The viewr typesetter renders `<x>` lowercase as a small-caps glyph (so
# `<h>` → "H", `<u>` → "U", etc.); spaces keep one-letter symbols 4-wide.
const _COVR_ELEM_NAMES = String[
    "<h> ","<h>e","<l>i","<b>e","<b> ","<c> ","<n> ","<o> ","<f> ",
    "<n>e","<n>a","<m>g","<a>l","<s>i","<p> ","<s> ","<c>l","<a>r",
    "<k> ","<c>a","<s>c","<t>i","<v> ","<c>r","<m>n","<f>e","<c>o",
    "<n>i","<c>u","<z>n","<g>a","<g>e","<a>s","<s>e","<b>r","<k>r",
    "<r>b","<s>r","<y> ","<z>r","<n>b","<m>o","<t>c","<r>u","<r>h",
    "<p>d","<a>g","<c>d","<i>n","<s>n","<s>b","<t>e","<i> ","<x>e",
    "<c>s","<b>a","<l>a","<c>e","<p>r","<n>d","<p>m","<s>m","<e>u",
    "<g>d","<t>b","<d>y","<h>o","<e>r","<t>m","<y>b","<l>u","<h>f",
    "<t>a","<w> ","<r>e","<o>s","<i>r","<p>t","<a>u","<h>g","<t>l",
    "<p>b","<b>i","<p>o","<a>t","<r>n","<f>r","<r>a","<a>c","<t>h",
    "<p>a","<u> ","<n>p","<p>u","<a>m","<c>m","<b>k","<c>f","<e>s",
    "<f>m","<m>d","<n>o","<l>w",
]

"""
    covr_elem(z::Integer) -> (name::String, len::Int)

Element name for Z (1..103) in viewr markup; returns the 4-char name and
its trimmed length (3 if the trailing char is a space, else 4). Matches
covr.f90:1753-1779 (`elem`).
"""
function covr_elem(z::Integer)
    1 <= z <= length(_COVR_ELEM_NAMES) ||
        error("covr_elem: Z=$z out of supported range 1..103 (covr.f90:1761-1773)")
    name = _COVR_ELEM_NAMES[z]
    len  = name[4] == ' ' ? 3 : 4
    (name, len)
end

# MT → reaction-name lookup tables. Match covr.f90:1795-1807 verbatim.
# `_MTNO_IRA1[k]` is a probe MT, `_MTNO_HIRA[k]` its right-hand name string,
# `_MTNO_IRA2[k]` the trimmed length of that string. mtno walks the table
# until ira1[k]==mt, falling back to "(mt   ) cont." when no entry hits.
const _MTNO_IRA1 = Int[
    1,2,3,4,16,17,18,22,28,37,51,91,102,103,104,105,106,107,111,
    207,780,781,452,455,456,25,24,32,33,41,112,115,116,
]
const _MTNO_IRA2 = Int[
    5,4,7,6,3,3,2,5,3,3,1,6,4,2,2,2,6,4,3,5,5,5,5,5,5,6,6,3,3,4,5,3,3,
]
const _MTNO_HIRA = String[
    "tot.)   ","el.)    ","nonel.) ","inel.)  ","2n)     ",
    "3n)     ","f)      ","n]a<)   ","np)     ","4n)     ",
    ")       ","cont.)  ","]g<)    ","p)      ","d)      ",
    "t)      ","<h>e3)  ","]a<)    ","2p      ","<h>e)   ",
    "]a<0)   ","]a<1)   "," ]n<)   "," ]n<)   "," ]n<)   ",
    "3n]a<)  ","2n]a<)  ","nd)     ","nt)     ","2np)    ",
    "p]a<)   ","pd)     ","pt)     ",
]
# Row 13 (`]g<)`) is the (n,γ) reaction marker (covr.f90:1807).

"""
    covr_mtno(mt::Integer; iverf::Integer=6) -> (lnamel, inamel, ivl, lnamer, inamer)

Convert an MT-number to reaction-name parts (left + level + right). Mirrors
covr.f90:1781-1890 (`mtno`). `iverf` is the ENDF version (5 or 6) read from
the MF1/MT451 directory; for ENDF-6 the proton/d/t/3He/α MT-ranges 600-849
are recognised, ENDF-5 they are not.

Returns:
- `lnamel` : 8-char left part, e.g. "(n,n    " or "(mt 052 "
- `inamel` : useful length of `lnamel`
- `ivl`    : discrete level number (-1 if not a level reaction; 0..99 else)
- `lnamer` : 8-char right part from the hira table (e.g. "tot.)   ")
- `inamer` : useful length of `lnamer`
"""
function covr_mtno(mt::Integer; iverf::Integer=6)
    blank   = "        "
    nmea1   = "(n,     "
    nmeh1   = "(mt     "
    nmed1   = "(total  "
    nmeb1   = "(n,n    "
    nmep1   = "(n,p    "
    nmez1   = "(n,d    "
    nmet1   = "(n,t    "
    nmex1   = "(n,h    "
    nmey1   = "(n,a    "
    nmee1   = "(delayed"
    nmef1   = "(prompt "
    nmeg1   = "(spectr."

    # Pick left name + length (covr.f90:1822-1860).
    if 50 < mt < 92
        lnamel = nmeb1; inamel = 4
    elseif 600 <= mt < 650 && iverf > 5
        lnamel = nmep1; inamel = 4
    elseif 650 <= mt < 700 && iverf > 5
        lnamel = nmez1; inamel = 4
    elseif 700 <= mt < 750 && iverf > 5
        lnamel = nmet1; inamel = 4
    elseif 750 <= mt < 800 && iverf > 5
        lnamel = nmex1; inamel = 4
    elseif 800 <= mt < 850 && iverf > 5
        lnamel = nmey1; inamel = 4
    elseif mt == 452
        lnamel = nmed1; inamel = 6
    elseif mt == 455
        lnamel = nmee1; inamel = 8
    elseif mt == 456
        lnamel = nmef1; inamel = 7
    elseif mt == 251
        lnamel = nmeh1; inamel = 3
    elseif mt == 261
        lnamel = nmeg1; inamel = 8
    else
        lnamel = nmea1; inamel = 3
    end

    # Discrete level number (covr.f90:1864-1870).
    ivl = -1
    if 50 < mt < 91
        ivl = mt - 50
    elseif 600 <= mt < 650 && iverf > 5
        ivl = mt - 600
    elseif 650 <= mt < 700 && iverf > 5
        ivl = mt - 650
    elseif 700 <= mt < 750 && iverf > 5
        ivl = mt - 700
    elseif 750 <= mt < 800 && iverf > 5
        ivl = mt - 750
    elseif 800 <= mt < 850 && iverf > 5
        ivl = mt - 800
    end

    # Right-hand name (covr.f90:1874-1888).  Walk ira1 until either an
    # exact MT match or, for level reactions (ivl in 0..99), the row whose
    # ira1 == 51 (which carries the inelastic level template).
    jloc = 0
    found = false
    for j in 1:33
        if _MTNO_IRA1[j] == mt
            jloc = j; found = true; break
        elseif ivl >= 0 && ivl <= 99 && _MTNO_IRA1[j] == 51
            jloc = j; found = true; break
        end
    end
    if !found
        # Fortran covr.f90:1881-1886: writes mt into cols 4..6 of nmeh1
        # (`(mt     `) with `i3` format. Cols 1..3 stay `(mt`, then 3 digits
        # (right-justified, width 3), then two trailing spaces. inamel=6 so
        # only the first 6 chars are used, giving `(mt252` style. The naive
        # concatenation `"(mt " * 252` would shift digits right and lose the
        # last digit when truncated to 6.
        lnamel = "(mt" * lpad(string(Int(mt)), 3) * "  "
        inamel = 6
        # Use jloc=11 (the "cont.)" entry) as right-name index — covr.f90:1885.
        jloc = 11
    end
    lnamer = _MTNO_HIRA[jloc]
    inamer = _MTNO_IRA2[jloc]

    (lnamel, inamel, ivl, lnamer, inamer)
end

"""
    covr_matmes(iza, mat, mt; iverf=6) -> (strng::String, lstrng::Int)

Build a viewer-markup reaction name like "#EH<235#HXEX<u>(n,f)".
Matches covr.f90:1705-1751 (`matmes`). `iza` is the ENDF ZA = 1000*Z + A.
"""
function covr_matmes(iza::Integer, mat::Integer, mt::Integer; iverf::Integer=6)
    iz = iza ÷ 1000
    ia = iza - iz * 1000
    lname, iname = covr_elem(iz)

    # Isotope mass-number markup (covr.f90:1720-1729). The #EH<…#HXEX
    # pair tells viewr to render the digits as a half-height super-script.
    if 0 < ia < 10
        liso = "#EH<" * string(ia) * "#HXEX"
        niso = 10
    elseif 10 <= ia < 100
        liso = "#EH<" * string(ia) * "#HXEX"
        niso = 11
    elseif ia >= 100
        liso = "#EH<" * string(ia) * "#HXEX"
        niso = 12
    else
        liso = ""
        niso = 0
    end

    if ia == 0
        lmat = lname[1:iname]
        nmat = iname
    else
        lmat = liso[1:niso] * lname[1:iname]
        nmat = niso + iname
    end

    lnamel, inamel, ivl, lnamer, inamer = covr_mtno(mt; iverf=iverf)

    # Discrete-level number markup (covr.f90:1738-1747).
    if ivl <= 0
        lnam = lnamer
        nnam = inamer
    elseif ivl < 10
        lnam = "#L.25H.75<" * string(ivl) * "#HXLX<)"
        nnam = 18
    else  # ivl < 100 (covr.f90 errors out if ≥ 100 — we let caller crash via guard upstream)
        lnam = "#L.25H.75<" * string(ivl) * "#HXLX<)"
        nnam = 19
    end

    strng = lmat[1:nmat] * lnamel[1:inamel] * lnam[1:nnam]
    lstrng = nmat + inamel + nnam
    (strng, lstrng)
end

"""
    covr_smilab(iza, mat, mt, izap, mtflg, mfflg, mf35, einc; iverf=6) -> String

Build the title string for one of the two semi-log frames in a covr plot
group. Mirrors covr.f90:1649-1703 (`smilab`). The five branches in `str1`
choose between "[d]s/s", "    ]s", "[d]n/n", "    ]n", etc.; `str2` picks
either "> vs. <e> for " or, for MF=5 cross-reaction frames, an inset
"(<e>#LH>in#HXLX>= …) ," form. The mfflg=-14 branches append ", MF40" or
", izap = NNN" markers (covr.f90:1693-1700).
"""
function covr_smilab(iza::Integer, mat::Integer, mt::Integer,
                     izap::Integer, mtflg::Integer,
                     mfflg::Integer, mf35::Integer, einc::Real;
                     iverf::Integer=6)
    # ---- left token (str1 in Fortran) ----
    l1 = 6
    if 452 <= mt <= 456 && mtflg == 0
        str1 = "[d]n/n"
    elseif 452 <= mt <= 456 && mtflg != 0
        str1 = "    ]n"
    elseif mt == 251 && mtflg != 0
        str1 = "    ]m"
    elseif mtflg != 0 && mf35 == 5
        str1 = "Grp-average ]f"; l1 = 14
    elseif mtflg == 0 && mf35 == 5
        str1 = "[d]f/f"
    elseif mtflg == 0 && mt == 251
        str1 = "[d]m/m"
    elseif mtflg != 0
        str1 = "    ]s"
    else
        str1 = "[d]s/s"
    end

    # ---- middle token (str2 in Fortran) ----
    if mtflg != 0 && mf35 == 5 && einc > 9.999e4
        str2 = ">(<e>#LH>in#HXLX>=" * @sprintf("%5.2f", einc/1.0e6) * " <m>e<v>), "
        l2 = 34
    elseif mtflg != 0 && mf35 == 5 && einc <= 9.999e4
        str2 = ">(<e>#LH>in#HXLX>= " * @sprintf("%9.2e", einc) * " e<v>), "
        l2 = 35
    else
        str2 = "> vs. <e> for "
        l2 = 14
    end

    # ---- reaction name (str3 in Fortran) ----
    str3, lstr3 = covr_matmes(iza, mat, mt; iverf=iverf)

    # ---- mfflg=-14 trailers (covr.f90:1693-1700) ----
    if mfflg != -14
        str1[1:l1] * str2[1:l2] * str3[1:lstr3]
    elseif mfflg == -14 && izap == 0
        str1[1:l1] * str2[1:l2] * str3[1:lstr3] * ", MF40"
    else
        str1[1:l1] * str2[1:l2] * str3[1:lstr3] * ", izap = " * lpad(string(Int(izap)), 7)
    end
end
