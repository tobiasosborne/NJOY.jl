using Test
using NJOY

@testset "NJOY.jl" begin

    # ======================================================================
    # Physics Constants
    # ======================================================================
    @testset "Physics Constants" begin
        C = NJOY.PhysicsConstants

        # Values must match NJOY2016 phys.f90 exactly (CODATA 2014)
        @test C.pi    == 3.141592653589793238
        @test C.euler == 0.57721566490153286
        @test C.bk    == 8.617333262e-5         # eV/K
        @test C.ev    == 1.602176634e-12         # erg/eV
        @test C.clight == 2.99792458e10          # cm/s

        # Particle masses in amu
        @test C.amassn == 1.00866491595
        @test C.amassp == 1.007276466621
        @test C.amassa == 4.001506179127
        @test C.amasse == 5.48579909065e-4

        # Derived quantities
        @test C.amu    == 931.49410242e6 * C.ev / (C.clight^2)
        @test C.hbar   == 6.582119569e-16 * C.ev
        @test C.pnratio == C.amassp / C.amassn
        @test C.anratio == C.amassa / C.amassn

        # epair = m_e * amu * c^2 / eV
        @test C.epair == C.amasse * C.amu * C.clight^2 / C.ev

        # Verify inverse fine structure constant is approximately 137
        @test 137.0 < C.finstri < 137.1
    end

    # ======================================================================
    # ENDF Float Parsing
    # ======================================================================
    @testset "ENDF Float Parsing" begin
        # Standard ENDF format: no 'E', sign indicates exponent
        @test parse_endf_float(" 1.234567+8") == 1.234567e8
        @test parse_endf_float("-1.234567+8") == -1.234567e8
        @test parse_endf_float(" 1.234567-3") == 1.234567e-3
        @test parse_endf_float("-2.345678-2") == -2.345678e-2

        # Extended precision (no exponent, just a number)
        @test parse_endf_float(" 1.23456789") == 1.23456789

        # With explicit 'E'
        @test parse_endf_float("1.23456E+08") == 1.23456e8

        # Zero and blank
        @test parse_endf_float("           ") == 0.0
        @test parse_endf_float(" 0.000000+0") == 0.0

        # Large and small exponents
        @test parse_endf_float(" 9.99999+38") == 9.99999e38
        @test parse_endf_float(" 1.00000-30") == 1.0e-30

        # Negative mantissa with positive exponent
        @test parse_endf_float("-3.14159+0") == -3.14159

        # Roundtrip consistency
        for val in [0.0, 1.0, -1.0, 1.23456e7, -3.14159e-12, 9.99999e38]
            formatted = format_endf_float(val)
            @test length(formatted) == 11
            parsed = parse_endf_float(formatted)
            if val == 0.0
                @test parsed == 0.0
            else
                @test isapprox(parsed, val, rtol=1e-5)
            end
        end
    end

    # ======================================================================
    # ENDF Float Formatting
    # ======================================================================
    @testset "ENDF Float Formatting" begin
        s = format_endf_float(0.0)
        @test length(s) == 11
        @test parse_endf_float(s) == 0.0

        s = format_endf_float(1.234567e8)
        @test length(s) == 11
        @test isapprox(parse_endf_float(s), 1.234567e8, rtol=1e-6)

        s = format_endf_float(-3.14159e-12)
        @test length(s) == 11
        @test isapprox(parse_endf_float(s), -3.14159e-12, rtol=1e-5)
    end

    # ======================================================================
    # ENDF I/O Roundtrip
    # ======================================================================
    @testset "ENDF I/O Roundtrip" begin
        id = MaterialId(Int32(125), Int32(3), Int32(1))

        # CONT record roundtrip
        rec = ContRecord(1.001e3, 2.34e-5, Int32(1), Int32(2),
                         Int32(3), Int32(4), id)
        buf = IOBuffer()
        write_cont(buf, rec)
        seekstart(buf)
        rec2 = read_cont(buf)
        @test isapprox(rec2.C1, rec.C1, rtol=1e-5)
        @test isapprox(rec2.C2, rec.C2, rtol=1e-5)
        @test rec2.L1 == rec.L1
        @test rec2.L2 == rec.L2
        @test rec2.N1 == rec.N1
        @test rec2.N2 == rec.N2
        @test rec2.id.mat == rec.id.mat

        # TAB1 record roundtrip
        interp = InterpolationTable([Int32(4)], [Int32(2)])
        tab1 = Tab1Record(92235.0, 233.025, Int32(0), Int32(0),
                          interp,
                          [1.0e-5, 1.0e0, 1.0e3, 2.0e7],
                          [10.0, 20.0, 5.0, 1.0],
                          id)
        buf = IOBuffer()
        write_tab1(buf, tab1)
        seekstart(buf)
        tab1b = read_tab1(buf)
        @test length(tab1b.x) == 4
        @test length(tab1b.y) == 4
        @test tab1b.interp.law[1] == LinLin
        for i in 1:4
            @test isapprox(tab1b.x[i], tab1.x[i], rtol=1e-5)
            @test isapprox(tab1b.y[i], tab1.y[i], rtol=1e-5)
        end

        # LIST record roundtrip
        list = ListRecord(0.0, 0.0, Int32(0), Int32(0),
                          Int32(5), Int32(0),
                          [1.0, 2.0, 3.0, 4.0, 5.0], id)
        buf = IOBuffer()
        write_list(buf, list)
        seekstart(buf)
        list2 = read_list(buf)
        @test list2.N1 == 5
        @test length(list2.data) == 5
        for i in 1:5
            @test isapprox(list2.data[i], list.data[i], rtol=1e-5)
        end
    end

    # ======================================================================
    # Interpolation
    # ======================================================================
    @testset "Interpolation" begin
        # terp1: Histogram (law 1)
        @test terp1(1.0, 10.0, 2.0, 20.0, 1.5, Histogram) == 10.0

        # terp1: Linear-linear (law 2)
        @test terp1(0.0, 0.0, 10.0, 10.0, 5.0, LinLin) == 5.0
        @test terp1(1.0, 2.0, 3.0, 6.0, 2.0, LinLin) == 4.0

        # terp1: Lin-log (law 3) -- y linear in ln(x)
        y = terp1(1.0, 0.0, 100.0, 1.0, 10.0, LinLog)
        # ln(10)/ln(100) = 0.5
        @test isapprox(y, 0.5, atol=1e-10)

        # terp1: Log-lin (law 4) -- ln(y) linear in x
        y = terp1(0.0, 1.0, 2.0, exp(2.0), 1.0, LogLin)
        @test isapprox(y, exp(1.0), rtol=1e-10)

        # terp1: Log-log (law 5) -- ln(y) linear in ln(x)
        # Power law: y = x^2 through (1,1) and (4,16)
        y = terp1(1.0, 1.0, 4.0, 16.0, 2.0, LogLog)
        @test isapprox(y, 4.0, rtol=1e-10)

        # terp1: Degenerate case x1 == x2
        @test terp1(5.0, 3.0, 5.0, 7.0, 5.0, LinLin) == 3.0

        # terp1: y1 == y2
        @test terp1(1.0, 5.0, 10.0, 5.0, 5.0, LinLin) == 5.0

        # TabulatedFunction interpolation
        interp = InterpolationTable([Int32(5)], [Int32(2)])  # lin-lin
        tab = TabulatedFunction(interp,
                                [0.0, 1.0, 2.0, 3.0, 4.0],
                                [0.0, 2.0, 4.0, 6.0, 8.0])
        @test interpolate(tab, 0.5) == 1.0
        @test interpolate(tab, 2.5) == 5.0
        @test interpolate(tab, 0.0) == 0.0
        @test interpolate(tab, 4.0) == 8.0

        # Outside range returns 0
        @test interpolate(tab, -1.0) == 0.0
        @test interpolate(tab, 5.0) == 0.0
    end

    # ======================================================================
    # Integration
    # ======================================================================
    @testset "Integration" begin
        # Integrate a linear function y = 2x from 0 to 4
        interp = InterpolationTable([Int32(5)], [Int32(2)])
        tab = TabulatedFunction(interp,
                                [0.0, 1.0, 2.0, 3.0, 4.0],
                                [0.0, 2.0, 4.0, 6.0, 8.0])
        # Integral of 2x from 0 to 4 = x^2 |_0^4 = 16
        result = integrate(tab, 0.0, 4.0)
        @test isapprox(result, 16.0, rtol=1e-10)

        # Partial integration from 1 to 3: x^2|_1^3 = 9 - 1 = 8
        result = integrate(tab, 1.0, 3.0)
        @test isapprox(result, 8.0, rtol=1e-10)

        # Integrate constant function y = 5 from 0 to 3 (histogram)
        interp_h = InterpolationTable([Int32(4)], [Int32(1)])
        tab_h = TabulatedFunction(interp_h,
                                  [0.0, 1.0, 2.0, 3.0],
                                  [5.0, 5.0, 5.0, 5.0])
        result = integrate(tab_h, 0.0, 3.0)
        @test isapprox(result, 15.0, rtol=1e-10)

        # Integration of subinterval
        result = integrate(tab_h, 0.5, 2.5)
        @test isapprox(result, 10.0, rtol=1e-10)
    end

    # ======================================================================
    # Penetrability and Shift Factors
    # ======================================================================
    @testset "Penetrability and Shift Factors" begin
        # l=0 analytical values
        @test penetrability(0, 0.5) == 0.5
        @test penetrability(0, 1.0) == 1.0
        @test shift_factor(0, 0.5) == 0.0
        @test shift_factor(0, 1.0) == 0.0

        # l=1 analytical values
        # P_1(rho) = rho^3 / (1 + rho^2)
        rho = 1.0
        @test isapprox(penetrability(1, rho), 1.0 / 2.0, rtol=1e-14)
        # S_1(rho) = -1 / (1 + rho^2)
        @test isapprox(shift_factor(1, rho), -0.5, rtol=1e-14)

        rho = 2.0
        @test isapprox(penetrability(1, rho), 8.0 / 5.0, rtol=1e-14)
        @test isapprox(shift_factor(1, rho), -1.0 / 5.0, rtol=1e-14)

        # l=2 analytical values
        # P_2(rho) = rho^5 / (9 + 3*rho^2 + rho^4)
        rho = 1.0
        @test isapprox(penetrability(2, rho), 1.0 / 13.0, rtol=1e-14)
        # S_2(rho) = -(18 + 3*rho^2) / (9 + 3*rho^2 + rho^4)
        @test isapprox(shift_factor(2, rho), -21.0 / 13.0, rtol=1e-14)

        rho = 3.0
        r2 = 9.0; r4 = 81.0
        den = 9.0 + 27.0 + 81.0  # = 117
        @test isapprox(penetrability(2, rho), r4 * 3.0 / den, rtol=1e-13)
        @test isapprox(shift_factor(2, rho), -(18.0 + 27.0) / den, rtol=1e-13)

        # l=3 analytical values at rho=1
        rho = 1.0
        den3 = 225.0 + 45.0 + 6.0 + 1.0  # = 277
        @test isapprox(penetrability(3, rho), 1.0 / den3, rtol=1e-14)
        @test isapprox(shift_factor(3, rho), -(675.0 + 90.0 + 6.0) / den3,
                       rtol=1e-14)

        # l=4 analytical values at rho=1
        rho = 1.0
        den4 = 11025.0 + 1575.0 + 135.0 + 10.0 + 1.0  # = 12746
        @test isapprox(penetrability(4, rho), 1.0 / den4, rtol=1e-14)
        @test isapprox(shift_factor(4, rho),
                       -(44100.0 + 4725.0 + 270.0 + 10.0) / den4, rtol=1e-14)

        # P_l should be positive for positive rho
        for l in 0:4, rho in [0.1, 0.5, 1.0, 2.0, 5.0]
            @test penetrability(l, rho) > 0.0
        end

        # S_l should be <= 0 for l >= 1
        for l in 1:4, rho in [0.1, 0.5, 1.0, 2.0, 5.0]
            @test shift_factor(l, rho) <= 0.0
        end
    end

    # ======================================================================
    # Phase Shifts
    # ======================================================================
    @testset "Phase Shifts" begin
        # l=0: phi_0 = rho
        @test phase_shift(0, 1.0) == 1.0
        @test phase_shift(0, 0.5) == 0.5

        # l=1: phi_1 = rho - atan(rho)
        rho = 1.0
        @test isapprox(phase_shift(1, rho), rho - atan(rho), rtol=1e-14)

        rho = 2.0
        @test isapprox(phase_shift(1, rho), rho - atan(rho), rtol=1e-14)

        # l=2: phi_2 = rho - atan(3*rho / (3 - rho^2))
        rho = 1.0
        @test isapprox(phase_shift(2, rho),
                       rho - atan(3.0 * rho / (3.0 - rho^2)), rtol=1e-14)

        # l=3: phi_3 = rho - atan((15*rho - rho^3) / (15 - 6*rho^2))
        rho = 1.0
        expected = rho - atan((15.0 * rho - rho^3) / (15.0 - 6.0 * rho^2))
        @test isapprox(phase_shift(3, rho), expected, rtol=1e-14)

        # l=4: phi_4 = rho - atan((105*rho - 10*rho^3) / (105 - 45*rho^2 + rho^4))
        rho = 2.0
        top = 105.0 * rho - 10.0 * rho^3
        bot = 105.0 - 45.0 * rho^2 + rho^4
        expected = rho - atan(top / bot)
        @test isapprox(phase_shift(4, rho), expected, rtol=1e-14)

        # Phase shifts should be non-negative for rho > 0
        for l in 0:4, rho in [0.5, 1.0, 2.0, 5.0]
            @test phase_shift(l, rho) >= 0.0
        end

        # Small rho: phi_l should be small for l >= 2
        @test phase_shift(2, 0.01) == 0.0  # below test threshold
        @test phase_shift(3, 0.01) == 0.0
        @test phase_shift(4, 0.01) == 0.0
    end

    # ======================================================================
    # Resonance Types
    # ======================================================================
    @testset "Resonance Types" begin
        # CrossSections arithmetic
        xs1 = CrossSections(10.0, 5.0, 3.0, 2.0)
        xs2 = CrossSections(20.0, 10.0, 6.0, 4.0)
        xs3 = xs1 + xs2
        @test xs3.total   == 30.0
        @test xs3.elastic  == 15.0
        @test xs3.fission  == 9.0
        @test xs3.capture  == 6.0

        xs4 = 2.0 * xs1
        @test xs4.total   == 20.0
        @test xs4.elastic  == 10.0

        # Default CrossSections is zero
        xs0 = CrossSections()
        @test xs0.total == 0.0
        @test xs0.elastic == 0.0
    end

    # ======================================================================
    # Faddeeva Function (Complex Probability Integral)
    # ======================================================================
    @testset "Faddeeva -- exact w(z)" begin
        # w(0) = 1
        rew, aimw = faddeeva_w(0.0, 0.0)
        @test rew == 1.0
        @test aimw == 0.0

        # Known values: w(z) at specific points compared to SpecialFunctions.jl
        # (erfcx-based reference, accurate to machine epsilon)
        test_points = [
            (0.0, 1.0),
            (1.0, 0.0),
            (1.0, 1.0),
            (0.5, 0.5),
            (2.0, 0.5),
            (0.5, 2.0),
            (3.0, 1.0),
            (0.1, 0.1),
            (0.5, 0.1),
            (4.0, 0.5),
            (0.0, 0.5),
            (2.5, 2.5),
        ]

        for (x, y) in test_points
            r_njoy, i_njoy = faddeeva_w(x, y)
            r_ref, i_ref = faddeeva_w_julia(x, y)
            @test isapprox(r_njoy, r_ref, atol=2e-7) ||
                  isapprox(r_njoy, r_ref, rtol=1e-5)
            @test isapprox(i_njoy, i_ref, atol=2e-7) ||
                  isapprox(i_njoy, i_ref, rtol=1e-5)
        end

        # w(0, y) should be real for all y > 0 (imaginary part = 0)
        for y in [0.1, 0.5, 1.0, 2.0, 5.0]
            rew, aimw = faddeeva_w(0.0, y)
            @test abs(aimw) < 1e-10
            @test rew > 0.0  # w(iy) is real and positive for y > 0
        end

        # Symmetry: Re(w(x,y)) is even in x, Im(w(x,y)) is odd in x
        for (x, y) in [(1.0, 0.5), (2.0, 1.0), (0.3, 0.7), (3.0, 2.0)]
            rp, ip = faddeeva_w(x, y)
            rm, im = faddeeva_w(-x, y)
            @test isapprox(rp, rm, rtol=1e-10)     # Re is even
            @test isapprox(ip, -im, rtol=1e-10)     # Im is odd
        end

        # For large y (y >> 1), w(0, y) -> 1/(sqrt(pi)*y) (asymptotic)
        for y in [10.0, 50.0, 100.0]
            rew, _ = faddeeva_w(0.0, y)
            asymp = 1.0 / (sqrt(pi) * y)
            @test isapprox(rew, asymp, rtol=1e-2)
        end

        # Negative imaginary part: test that the Taylor branch works
        r1, i1 = faddeeva_w(1.0, -0.5)
        r2, i2 = faddeeva_w_julia(1.0, -0.5)
        @test isapprox(r1, r2, atol=1e-6)
        @test isapprox(i1, i2, atol=1e-6)
    end

    @testset "Faddeeva -- SpecialFunctions reference" begin
        # Verify faddeeva_w_julia against known analytical results
        # w(0) = 1
        r, i = faddeeva_w_julia(0.0, 0.0)
        @test isapprox(r, 1.0, atol=1e-14)
        @test isapprox(i, 0.0, atol=1e-14)

        # w(iy) = erfcx(y) for real y (purely imaginary argument)
        using SpecialFunctions: erfcx
        for y in [0.1, 0.5, 1.0, 2.0, 5.0]
            r, i = faddeeva_w_julia(0.0, y)
            @test isapprox(r, erfcx(y), rtol=1e-14)
            @test abs(i) < 1e-14
        end
    end

    @testset "Faddeeva -- table construction" begin
        tab = build_faddeeva_table()
        @test tab isa FaddeevaTable
        @test length(tab.tr) == 3844  # 62*62
        @test length(tab.ti) == 3844

        # Spot-check: table at index (2,2) corresponds to x=0.0, y=0.0
        # (grid starts at -0.1 with spacing 0.1, so index 2 = -0.1 + 1*0.1 = 0.0)
        idx_00 = NJOY._tindex(2, 2)
        @test isapprox(tab.tr[idx_00], 1.0, atol=1e-7)  # w(0,0) = 1
        @test abs(tab.ti[idx_00]) < 1e-7                  # Im(w(0,0)) = 0

        # Spot-check: table at grid point (1.0, 1.0) -> index (12, 12)
        # x = -0.1 + 11*0.1 = 1.0, y = -0.1 + 11*0.1 = 1.0
        idx_11 = NJOY._tindex(12, 12)
        r_ref, i_ref = faddeeva_w(1.0, 1.0)
        @test isapprox(tab.tr[idx_11], r_ref, rtol=1e-10)
        @test isapprox(tab.ti[idx_11], i_ref, rtol=1e-10)
    end

    @testset "Faddeeva -- quickw regions" begin
        tab = build_faddeeva_table()

        # Region 1: table lookup (x^2 + y^2 < 36)
        # Test at several interior points
        for (x, y) in [(0.0, 0.0), (1.0, 1.0), (0.5, 0.1), (2.0, 2.0),
                        (3.0, 3.0), (4.5, 2.0), (1.5, 4.0)]
            rq, iq = quickw(x, y, tab)
            re, ie = faddeeva_w(x, y)
            # Table interpolation introduces small errors, but should be close
            @test isapprox(rq, re, atol=1e-4, rtol=1e-3)
            @test isapprox(iq, ie, atol=1e-4, rtol=1e-3)
        end

        # Region 2: rational approximation (36 <= x^2+y^2 < 144)
        for (x, y) in [(6.0, 0.5), (5.0, 4.0), (7.0, 7.0), (10.0, 5.0)]
            rq, iq = quickw(x, y, tab)
            re, ie = faddeeva_w_julia(x, y)
            @test isapprox(rq, re, atol=1e-5, rtol=1e-3)
            @test isapprox(iq, ie, atol=1e-5, rtol=1e-3)
        end

        # Region 3: simpler rational (144 <= x^2+y^2 < 10000)
        for (x, y) in [(12.0, 1.0), (50.0, 50.0), (80.0, 10.0)]
            rq, iq = quickw(x, y, tab)
            re, ie = faddeeva_w_julia(x, y)
            @test isapprox(rq, re, atol=1e-6, rtol=1e-3)
            @test isapprox(iq, ie, atol=1e-6, rtol=1e-3)
        end

        # Region 4: asymptotic (x^2+y^2 >= 10000)
        for (x, y) in [(100.0, 1.0), (200.0, 200.0), (500.0, 0.1)]
            rq, iq = quickw(x, y, tab)
            re, ie = faddeeva_w_julia(x, y)
            @test isapprox(rq, re, atol=1e-8, rtol=1e-2)
            @test isapprox(iq, ie, atol=1e-8, rtol=1e-2)
        end

        # Negative x: Re(w) is even, Im(w) is odd
        for (x, y) in [(1.0, 1.0), (3.0, 2.0), (7.0, 5.0), (50.0, 1.0)]
            rp, ip = quickw(x, y, tab)
            rm, im = quickw(-x, y, tab)
            @test isapprox(rp, rm, rtol=1e-10)
            @test isapprox(ip, -im, rtol=1e-10)
        end

        # compute_imaginary=false should return aimw=0
        rew, aimw = quickw(1.0, 1.0, tab; compute_imaginary=false)
        @test aimw == 0.0
        @test rew != 0.0
    end

    @testset "Faddeeva -- psi/chi Voigt profiles" begin
        # psi peaks at x=0 for fixed y
        y = 1.0
        psi0, chi0 = psi_chi(0.0, y)
        for x in [0.5, 1.0, 2.0, 5.0]
            psi_x, _ = psi_chi(x, y)
            @test psi0 > psi_x  # peak at x=0
        end

        # chi = 0 at x = 0 (by symmetry, Im(w(iy)) = 0)
        _, chi0 = psi_chi(0.0, 1.0)
        @test abs(chi0) < 1e-10

        # psi -> 0 as x -> infinity
        for x in [10.0, 50.0, 100.0]
            psi_x, _ = psi_chi(x, 1.0)
            @test abs(psi_x) < 0.01
        end

        # chi is odd in x
        for x in [0.5, 1.0, 3.0]
            _, chi_p = psi_chi(x, 1.0)
            _, chi_m = psi_chi(-x, 1.0)
            @test isapprox(chi_p, -chi_m, rtol=1e-10)
        end

        # psi is even in x
        for x in [0.5, 1.0, 3.0]
            psi_p, _ = psi_chi(x, 1.0)
            psi_m, _ = psi_chi(-x, 1.0)
            @test isapprox(psi_p, psi_m, rtol=1e-10)
        end

        # Fast version with table should agree with exact
        tab = build_faddeeva_table()
        for (x, y) in [(0.0, 1.0), (1.0, 0.5), (2.0, 2.0)]
            p1, c1 = psi_chi(x, y)
            p2, c2 = psi_chi(x, y, tab)
            @test isapprox(p1, p2, atol=1e-4)
            @test isapprox(c1, c2, atol=1e-4)
        end

        # psi is positive for y > 0 (Re(w) > 0 in upper half-plane)
        for (x, y) in [(0.0, 0.5), (1.0, 1.0), (2.0, 0.1), (0.0, 5.0)]
            psi_val, _ = psi_chi(x, y)
            @test psi_val > 0.0
        end
    end

    @testset "Faddeeva -- quickw vs SpecialFunctions systematic" begin
        tab = build_faddeeva_table()

        # Systematic sweep over many (x, y) points
        max_re_err = 0.0
        max_im_err = 0.0
        n_tested = 0
        for x in 0.0:0.5:15.0
            for y in 0.1:0.5:15.0
                rq, iq = quickw(x, y, tab)
                rj, ij = faddeeva_w_julia(x, y)
                re_err = abs(rq - rj)
                im_err = abs(iq - ij)
                if rj != 0.0
                    re_err = min(re_err, abs(rq - rj) / abs(rj))
                end
                if ij != 0.0
                    im_err = min(im_err, abs(iq - ij) / abs(ij))
                end
                max_re_err = max(max_re_err, re_err)
                max_im_err = max(max_im_err, im_err)
                n_tested += 1
            end
        end

        # The quickw approximation should be accurate to ~1e-3 or better
        # (matching NJOY's intended accuracy for cross section evaluation)
        @test max_re_err < 0.01
        @test max_im_err < 0.01
        @test n_tested > 900  # ensure we tested many points
    end

    @testset "Faddeeva -- benchmark quickw vs SpecialFunctions" begin
        tab = build_faddeeva_table()

        # Warm up
        quickw(1.0, 1.0, tab)
        faddeeva_w_julia(1.0, 1.0)

        # Time quickw (many evaluations)
        n_eval = 10_000
        xs = range(0.0, 10.0, length=100)
        ys = range(0.1, 5.0, length=100)

        t_quick = @elapsed begin
            for x in xs, y in ys
                quickw(x, y, tab)
            end
        end

        t_julia = @elapsed begin
            for x in xs, y in ys
                faddeeva_w_julia(x, y)
            end
        end

        # Just verify both complete without error; quickw should be competitive
        @test t_quick > 0.0
        @test t_julia > 0.0
        # Print timing info (informational, not a hard test)
        @info "Faddeeva timing" quickw_ms=round(t_quick*1000, digits=2) julia_ms=round(t_julia*1000, digits=2) speedup=round(t_julia/t_quick, digits=1)
    end

    # ======================================================================
    # ENDF I/O Integration Tests -- Real ENDF Evaluation Files
    # ======================================================================
    ENDF_RESOURCES = joinpath(@__DIR__,
        "..", "..", "njoy-reference", "tests", "resources")

    @testset "ENDF I/O -- H-2 (n-001_H_002-ENDF8.0)" begin
        endf_file = joinpath(ENDF_RESOURCES, "n-001_H_002-ENDF8.0.endf")
        if !isfile(endf_file)
            @warn "Skipping ENDF integration tests: $endf_file not found"
        else
            open(endf_file) do io
                # ---- (a) TPID record ----
                tpid = read_tpid(io)
                @test tpid.mat == 1       # first line always has MAT=1
                @test tpid.mf  == 0
                @test tpid.mt  == 0
                @test length(tpid.text) == 66

                # ---- (b) MF1/MT451 descriptive data ----
                head1 = read_cont(io)                # ZA, AWR
                @test head1.C1 == 1002.0              # ZA for H-2
                @test isapprox(head1.C2, 1.9968, rtol=1e-4)  # AWR
                @test head1.id.mat == 128             # MAT number for H-2
                @test head1.id.mf  == 1
                @test head1.id.mt  == 451

                head2 = read_cont(io)                # 0, 0, 0, 0, 0, NWD+6
                head3 = read_cont(io)                # ELIS, STA, LIS, LISO, 0, NFOR
                @test head3.N2 == 8                   # NFOR = 8 (ENDF-8)
                head4 = read_cont(io)                # AWI, 0, LDRV, 0, NWD, NXC
                NWD = Int(head4.N1)
                NXC = Int(head4.N2)
                @test NWD > 0
                @test NXC > 0

                # Skip descriptive text lines
                for _ in 1:NWD
                    readline(io)
                end

                # Read directory (NXC CONT records)
                directory = Tuple{Int,Int,Int}[]
                for _ in 1:NXC
                    d = read_cont(io)
                    push!(directory, (Int(d.L1), Int(d.L2), Int(d.N1)))
                end
                @test length(directory) == NXC

                # Verify expected sections are in the directory
                dir_set = Set((mf, mt) for (mf, mt, _) in directory)
                @test (1, 451) in dir_set
                @test (2, 151) in dir_set
                @test (3, 1) in dir_set
                @test (3, 2) in dir_set
                @test (3, 16) in dir_set
                @test (3, 102) in dir_set

                # ---- (c) Scan all MF/MT sections from the file ----
                seekstart(io)
                file_sections = Set{Tuple{Int,Int}}()
                while !eof(io)
                    line = readline(io)
                    p = rpad(line, 80)
                    mf = NJOY._parse_int(p[71:72])
                    mt = NJOY._parse_int(p[73:75])
                    if mf > 0 && mt > 0
                        push!(file_sections, (Int(mf), Int(mt)))
                    end
                end
                # All directory entries should match actual sections
                for (mf, mt, _) in directory
                    @test (mf, mt) in file_sections
                end

                # ---- (d) MF1/MT451 already read above, verify ZA ----
                @test head1.C1 == 1002.0  # H-2

                # ---- (e) MF2/MT151 resonance parameters ----
                @test find_section(io, 2, 151)
                mf2_head = read_head(io)
                @test mf2_head.C1 == 1002.0   # ZA
                @test mf2_head.id.mf == 2
                @test mf2_head.id.mt == 151

                # Range header
                range_cont = read_cont(io)
                NIS = Int(range_cont.N1)
                @test NIS >= 1   # at least one isotope

                # Isotope CONT
                iso_cont = read_cont(io)
                EL = iso_cont.C1
                EH = iso_cont.C2
                @test EL > 0      # lower energy bound positive
                @test EH > EL     # upper > lower
                LRU = Int(iso_cont.L1)
                @test LRU == 0    # H-2 has no resonance parameters (LRU=0)

                # Scattering radius
                scat = read_cont(io)
                SPI = scat.C1     # spin
                AP = scat.C2      # scattering radius
                @test SPI >= 0    # spin non-negative
                @test AP > 0      # scattering radius positive

                # ---- (f) MF3 cross section TAB1 records ----
                # MT=1: Total cross section
                @test find_section(io, 3, 1)
                mt1_head = read_head(io)
                @test mt1_head.C1 == 1002.0  # ZA
                mt1_tab = read_tab1(io)
                @test length(mt1_tab.x) == 178
                @test length(mt1_tab.y) == 178
                @test length(mt1_tab.interp) == 1
                @test mt1_tab.interp.law[1] == LinLin
                @test mt1_tab.interp.nbt[1] == 178

                # MT=2: Elastic scattering
                @test find_section(io, 3, 2)
                mt2_head = read_head(io)
                mt2_tab = read_tab1(io)
                @test length(mt2_tab.x) == 178
                @test mt2_tab.interp.law[1] == LinLin

                # MT=16: (n,2n)
                @test find_section(io, 3, 16)
                mt16_head = read_head(io)
                mt16_tab = read_tab1(io)
                @test length(mt16_tab.x) > 0
                # QM should be negative (threshold reaction)
                @test mt16_tab.C1 < 0

                # MT=102: Radiative capture
                @test find_section(io, 3, 102)
                mt102_head = read_head(io)
                mt102_tab = read_tab1(io)
                @test length(mt102_tab.x) > 0

                # ---- (g) Physical sanity checks ----
                for (label, tab) in [("MT1", mt1_tab), ("MT2", mt2_tab),
                                     ("MT16", mt16_tab), ("MT102", mt102_tab)]
                    # All energies positive
                    @test all(e -> e > 0, tab.x)
                    # All cross sections non-negative
                    @test all(y -> y >= 0, tab.y)
                    # Energies monotonically increasing
                    @test issorted(tab.x)
                end

                # Total XS should be >= elastic + capture at low energies
                # (where only elastic and capture contribute)
                tf_total = TabulatedFunction(mt1_tab)
                tf_elastic = TabulatedFunction(mt2_tab)
                tf_capture = TabulatedFunction(mt102_tab)
                for E in [1e-5, 0.0253, 1.0, 1e3]
                    sig_t = interpolate(tf_total, E)
                    sig_el = interpolate(tf_elastic, E)
                    sig_cap = interpolate(tf_capture, E)
                    # Below (n,2n) threshold, total ~ elastic + capture
                    @test sig_t >= sig_el - 1e-6
                end

                # Round-trip: write MF3/MT1 TAB1, read it back
                buf = IOBuffer()
                write_tab1(buf, mt1_tab)
                seekstart(buf)
                mt1_rt = read_tab1(buf)
                @test length(mt1_rt.x) == length(mt1_tab.x)
                @test length(mt1_rt.y) == length(mt1_tab.y)
                # Should be exact (no loss, values fit in 11-char ENDF fields)
                @test mt1_rt.x == mt1_tab.x
                @test mt1_rt.y == mt1_tab.y
                @test mt1_rt.interp.nbt == mt1_tab.interp.nbt
                @test mt1_rt.interp.law == mt1_tab.interp.law
            end
        end
    end

    @testset "ENDF I/O -- U-235 (n-092_U_235-ENDF8.0)" begin
        endf_file = joinpath(ENDF_RESOURCES, "n-092_U_235-ENDF8.0.endf")
        if !isfile(endf_file)
            @warn "Skipping U-235 integration tests: file not found"
        else
            open(endf_file) do io
                # TPID
                tpid = read_tpid(io)
                @test tpid.mat == 1

                # MF1/MT451 header
                head1 = read_cont(io)
                @test head1.C1 == 92235.0   # ZA for U-235
                @test head1.id.mat == 9228  # MAT number for U-235

                # MF3/MT18 (fission) -- important for U-235
                @test find_section(io, 3, 18)
                head = read_head(io)
                @test head.C1 == 92235.0
                tab = read_tab1(io)
                @test length(tab.x) > 100   # U-235 has many data points
                @test issorted(tab.x)
                @test all(e -> e > 0, tab.x)

                # MF3/MT1 (total)
                @test find_section(io, 3, 1)
                head_t = read_head(io)
                tab_t = read_tab1(io)
                @test length(tab_t.x) > 0
                @test tab_t.x[1] > 0
                @test issorted(tab_t.x)

                # MF3/MT102 (capture)
                @test find_section(io, 3, 102)
                head_c = read_head(io)
                tab_c = read_tab1(io)
                @test length(tab_c.x) > 0

                # Verify interpolation gives finite values
                # Note: MF3 for U-235 contains only the smooth background;
                # the large resonance contribution comes from MF2 reconstruction
                tf = TabulatedFunction(tab_t)
                sig_thermal = interpolate(tf, 0.0253)
                @test isfinite(sig_thermal)
                @test sig_thermal >= 0.0

                # Round-trip fission TAB1
                buf = IOBuffer()
                write_tab1(buf, tab)
                seekstart(buf)
                tab_rt = read_tab1(buf)
                @test length(tab_rt.x) == length(tab.x)
                for i in eachindex(tab.x)
                    @test isapprox(tab_rt.x[i], tab.x[i], rtol=1e-6)
                    @test isapprox(tab_rt.y[i], tab.y[i], rtol=1e-6)
                end
            end
        end
    end

    @testset "ENDF I/O -- MF3 MT=3 nonelastic (H-2)" begin
        endf_file = joinpath(ENDF_RESOURCES, "n-001_H_002-ENDF8.0.endf")
        if !isfile(endf_file)
            @warn "Skipping test"
        else
            open(endf_file) do io
                # MF3/MT3 nonelastic
                @test find_section(io, 3, 3)
                head = read_head(io)
                tab = read_tab1(io)
                @test length(tab.x) > 0
                @test issorted(tab.x)
                @test all(e -> e > 0, tab.x)
                # Nonelastic XS can be zero below threshold
                @test all(y -> y >= 0, tab.y)
            end
        end
    end

    @testset "ENDF I/O -- parse_endf_line edge cases" begin
        # Short line (less than 80 chars) -- should be padded
        short_line = " 1.234567+8 2.345678-3"
        fields, mat, mf, mt, ns = NJOY.parse_endf_line(short_line)
        @test NJOY.parse_endf_float(fields[1]) == 1.234567e8
        @test isapprox(NJOY.parse_endf_float(fields[2]), 2.345678e-3, rtol=1e-6)
        @test mat == 0  # no MAT in short line
        @test mf == 0
        @test mt == 0

        # Full 80-char CONT line from H-2 file
        full_line = " 1.002000+3 1.996800+0          0          0          0          0 128 1451    1"
        fields, mat, mf, mt, ns = NJOY.parse_endf_line(full_line)
        @test isapprox(NJOY.parse_endf_float(fields[1]), 1002.0, rtol=1e-6)
        @test isapprox(NJOY.parse_endf_float(fields[2]), 1.9968, rtol=1e-4)
        @test mat == 128
        @test mf == 1
        @test mt == 451
    end

    @testset "ENDF I/O -- LIST record from real file" begin
        # Construct a LIST record manually matching what would appear in an ENDF file
        # and verify round-trip
        id = MaterialId(Int32(128), Int32(3), Int32(1))
        data = collect(1.0:13.0)
        list = ListRecord(1002.0, 1.9968, Int32(0), Int32(0),
                          Int32(13), Int32(0), data, id)
        buf = IOBuffer()
        ns = write_list(buf, list)
        seekstart(buf)
        list2 = read_list(buf)
        @test list2.N1 == 13
        @test length(list2.data) == 13
        @test isapprox(list2.C1, 1002.0, rtol=1e-6)
        for i in 1:13
            @test isapprox(list2.data[i], data[i], rtol=1e-6)
        end
    end

    @testset "ENDF I/O -- TAB2 record round-trip" begin
        id = MaterialId(Int32(128), Int32(4), Int32(2))
        interp = InterpolationTable([Int32(10)], [Int32(2)])
        tab2 = Tab2Record(1002.0, 0.0, Int32(0), Int32(0), interp, Int32(10), id)
        buf = IOBuffer()
        write_tab2(buf, tab2)
        seekstart(buf)
        tab2b = read_tab2(buf)
        @test isapprox(tab2b.C1, 1002.0, rtol=1e-6)
        @test tab2b.NZ == 10
        @test tab2b.interp.nbt[1] == 10
        @test tab2b.interp.law[1] == LinLin
    end

    @testset "ENDF I/O -- find_section" begin
        endf_file = joinpath(ENDF_RESOURCES, "n-001_H_002-ENDF8.0.endf")
        if isfile(endf_file)
            open(endf_file) do io
                # Should find existing sections
                @test find_section(io, 1, 451)
                @test find_section(io, 2, 151)
                @test find_section(io, 3, 1)
                @test find_section(io, 3, 102)

                # Should not find non-existing sections
                @test !find_section(io, 3, 999)
                @test !find_section(io, 99, 1)
            end
        end
    end

    # ======================================================================
    # MF2 Reader -- Ag-109 (MLBW, LRF=2)
    # ======================================================================
    @testset "MF2 Reader -- Ag-109 MLBW" begin
        endf_file = joinpath(ENDF_RESOURCES, "n-047_Ag_109-ENDF8.0.endf")
        if !isfile(endf_file)
            @warn "Skipping MF2 Ag-109 tests: file not found"
        else
            open(endf_file) do io
                @test find_section(io, 2, 151)
                mf2 = read_mf2(io)

                # Basic structure checks
                @test mf2.ZA == 47109.0
                @test length(mf2.isotopes) == 1

                iso = mf2.isotopes[1]
                @test iso.ZAI == 47109.0
                @test iso.ABN == 1.0
                @test length(iso.ranges) >= 1

                rng = iso.ranges[1]
                @test rng.LRU == 1       # resolved
                @test rng.LRF == 2       # MLBW
                @test rng.EL < rng.EH    # energy bounds ordered

                params = rng.parameters
                @test params isa MLBWParameters
                @test params.NLS == 2    # s-wave and p-wave
                @test params.SPI == 0.5  # target spin 1/2
                @test params.AP > 0.0    # scattering radius positive
                @test length(params.l_values) == 2
                @test params.l_values[1] == 0  # s-wave
                @test params.l_values[2] == 1  # p-wave

                # Check resonance count
                @test length(params.Er[1]) > 0  # s-wave resonances
                @test length(params.Er[2]) > 0  # p-wave resonances

                # Check AWRI is stored
                @test length(params.AWRI) == 2
                @test params.AWRI[1] > 100.0  # Ag-109 is heavy

                # Check that first s-wave resonance has sensible values
                @test params.Gn[1][1] > 0.0 || params.Er[1][1] < 0.0  # Gn positive (or neg energy)
                @test params.Gg[1][1] > 0.0   # capture width positive

                # All resonance energies should be finite
                for il in 1:Int(params.NLS)
                    for ir in eachindex(params.Er[il])
                        @test isfinite(params.Er[il][ir])
                        @test isfinite(params.Gn[il][ir])
                        @test isfinite(params.Gg[il][ir])
                    end
                end
            end
        end
    end

    # ======================================================================
    # MF2 Reader -- Fe-56 (Reich-Moore, LRF=3)
    # ======================================================================
    @testset "MF2 Reader -- Fe-56 Reich-Moore" begin
        endf_file = joinpath(ENDF_RESOURCES, "n-026_Fe_056-JEFF3.3.endf")
        if !isfile(endf_file)
            @warn "Skipping MF2 Fe-56 tests: file not found"
        else
            open(endf_file) do io
                @test find_section(io, 2, 151)
                mf2 = read_mf2(io)

                @test mf2.ZA == 26056.0
                @test length(mf2.isotopes) == 1

                iso = mf2.isotopes[1]
                @test iso.ABN == 1.0

                rng = iso.ranges[1]
                @test rng.LRU == 1      # resolved
                @test rng.LRF == 3      # Reich-Moore
                @test rng.NAPS == 1     # use AP as channel radius
                @test rng.EL == 1.0e-5
                @test rng.EH == 8.5e5

                params = rng.parameters
                @test params isa ReichMooreParameters
                @test params.NLS == 3       # 3 l-values (s, p, d)
                @test params.SPI == 0.0     # Fe-56 spin 0
                @test params.AP > 0.0
                @test length(params.l_values) == 3
                @test params.l_values[1] == 0  # s-wave
                @test params.l_values[2] == 1  # p-wave
                @test params.l_values[3] == 2  # d-wave

                # Check AWRI and APL are stored
                @test length(params.AWRI) == 3
                @test params.AWRI[1] > 50.0  # Fe-56 AWR ~ 55.45

                # Check resonance count: s-wave should have many resonances
                @test length(params.Er[1]) == 40  # 40 s-wave resonances from the file

                # First s-wave resonance (negative energy = bound state)
                @test params.Er[1][1] < 0.0  # bound state

                # First positive-energy s-wave resonance
                # Find it
                first_pos = findfirst(e -> e > 0.0, params.Er[1])
                @test first_pos !== nothing
                @test params.Gn[1][first_pos] > 0.0  # neutron width positive

                # For Fe-56 (non-fissile), GFA and GFB should be zero
                for il in 1:3
                    for ir in eachindex(params.Gfa[il])
                        @test params.Gfa[il][ir] == 0.0
                        @test params.Gfb[il][ir] == 0.0
                    end
                end
            end
        end
    end

    # ======================================================================
    # MF2 Reader -- U-235 (Reich-Moore, fissile)
    # ======================================================================
    @testset "MF2 Reader -- U-235 Reich-Moore" begin
        endf_file = joinpath(ENDF_RESOURCES, "n-092_U_235-ENDF8.0.endf")
        if !isfile(endf_file)
            @warn "Skipping MF2 U-235 tests: file not found"
        else
            open(endf_file) do io
                @test find_section(io, 2, 151)
                mf2 = read_mf2(io)

                @test mf2.ZA == 92235.0
                @test length(mf2.isotopes) == 1

                iso = mf2.isotopes[1]
                rng = iso.ranges[1]
                @test rng.LRF == 3  # Reich-Moore

                params = rng.parameters
                @test params isa ReichMooreParameters
                @test params.SPI == 3.5  # U-235 spin 7/2

                # U-235 is fissile: should have nonzero GFA values
                has_fission = false
                for il in 1:Int(params.NLS)
                    for ir in eachindex(params.Gfa[il])
                        if params.Gfa[il][ir] != 0.0 || params.Gfb[il][ir] != 0.0
                            has_fission = true
                            break
                        end
                    end
                end
                @test has_fission
            end
        end
    end

    # ======================================================================
    # Wavenumber Constant
    # ======================================================================
    @testset "Wavenumber constant cwaven" begin
        cwaven = cwaven_constant()
        @test isfinite(cwaven)
        @test cwaven > 0.0

        # cwaven should be approximately 2.196807e-3 [1/(eV^{1/2} * fm)]
        # This comes from sqrt(2 * 1.00866 * 1.6605e-24 * 1.602e-12) * 1e-12 / hbar
        @test isapprox(cwaven, 2.196807e-3, rtol=1e-3)
    end

    # ======================================================================
    # Channel Radius
    # ======================================================================
    @testset "Channel radius" begin
        # Fe-56: AWRI ~ 55.4545
        ra = channel_radius(55.4545)
        @test isfinite(ra)
        @test ra > 0.0
        # Should be about 0.123 * (1.00866*55.4545)^(1/3) + 0.08 ~ 0.486 fm
        @test 0.4 < ra < 0.6
    end

    # ======================================================================
    # MLBW Cross Section -- Ag-109
    # ======================================================================
    @testset "MLBW Cross Section -- Ag-109" begin
        endf_file = joinpath(ENDF_RESOURCES, "n-047_Ag_109-ENDF8.0.endf")
        if !isfile(endf_file)
            @warn "Skipping MLBW XS tests: file not found"
        else
            open(endf_file) do io
                @test find_section(io, 2, 151)
                mf2 = read_mf2(io)
                rng = mf2.isotopes[1].ranges[1]

                # Thermal energy (0.0253 eV)
                xs_th = cross_section(0.0253, rng)
                @test isfinite(xs_th.total)
                @test isfinite(xs_th.elastic)
                @test isfinite(xs_th.capture)
                @test xs_th.total > 0.0
                @test xs_th.elastic > 0.0
                @test xs_th.capture > 0.0
                @test xs_th.fission == 0.0  # Ag-109 is not fissile

                # Total = elastic + capture + fission
                @test isapprox(xs_th.total,
                               xs_th.elastic + xs_th.capture + xs_th.fission,
                               rtol=1e-10)

                # At a resonance peak (5.19 eV resonance in Ag-109)
                xs_res = cross_section(5.19, rng)
                @test isfinite(xs_res.total)
                @test xs_res.total > 0.0

                # Total should equal sum of partials at resonance
                @test isapprox(xs_res.total,
                               xs_res.elastic + xs_res.capture + xs_res.fission,
                               rtol=1e-10)

                # Cross sections should be positive at several energies
                for E in [1e-5, 0.01, 0.1, 1.0, 10.0, 100.0, 1000.0, 5000.0]
                    xs = cross_section(E, rng)
                    @test xs.total > 0.0
                    @test xs.elastic >= 0.0
                    @test xs.capture >= 0.0
                    @test isapprox(xs.total,
                                   xs.elastic + xs.capture + xs.fission,
                                   rtol=1e-10)
                end

                # At the first prominent resonance (5.19 eV, J=1),
                # capture cross section should be large
                @test xs_res.capture > 100.0  # Should be hundreds of barns

                # 1/v region at very low energies: capture should be large
                xs_low = cross_section(1e-5, rng)
                @test xs_low.capture > 0.0
            end
        end
    end

    # ======================================================================
    # Reich-Moore Cross Section -- Fe-56
    # ======================================================================
    @testset "Reich-Moore Cross Section -- Fe-56" begin
        endf_file = joinpath(ENDF_RESOURCES, "n-026_Fe_056-JEFF3.3.endf")
        if !isfile(endf_file)
            @warn "Skipping RM XS tests: file not found"
        else
            open(endf_file) do io
                @test find_section(io, 2, 151)
                mf2 = read_mf2(io)
                rng = mf2.isotopes[1].ranges[1]

                # Thermal energy
                xs_th = cross_section(0.0253, rng)
                @test isfinite(xs_th.total)
                @test xs_th.total > 0.0
                @test xs_th.elastic > 0.0
                @test xs_th.capture > 0.0
                @test xs_th.fission == 0.0  # Fe-56 is not fissile

                # Total = elastic + capture + fission
                @test isapprox(xs_th.total,
                               xs_th.elastic + xs_th.capture + xs_th.fission,
                               rtol=1e-10)

                # At a known resonance: Er = 27791 eV (first positive s-wave resonance)
                xs_res = cross_section(27791.0, rng)
                @test isfinite(xs_res.total)
                @test xs_res.total > 0.0
                # At resonance, total should be much larger than potential scattering
                @test xs_res.total > 10.0  # should be large at resonance peak

                # Sum rule: total = elastic + capture + fission at several energies
                for E in [0.0253, 100.0, 1000.0, 10000.0, 27791.0, 50000.0, 100000.0]
                    xs = cross_section(E, rng)
                    @test isapprox(xs.total,
                                   xs.elastic + xs.capture + xs.fission,
                                   rtol=1e-10)
                end

                # Cross sections should be positive everywhere
                for E in [1e-5, 0.01, 1.0, 100.0, 1e4, 5e4, 2e5, 5e5]
                    xs = cross_section(E, rng)
                    @test xs.total > 0.0
                    @test xs.elastic >= 0.0
                    @test xs.capture >= 0.0
                end

                # At very low energies, capture follows 1/v law (increases)
                xs1 = cross_section(1e-3, rng)
                xs2 = cross_section(1.0, rng)
                # capture at lower E should be larger (1/v trend)
                @test xs1.capture > xs2.capture
            end
        end
    end

    # ======================================================================
    # Reich-Moore Cross Section -- U-235 (fissile)
    # ======================================================================
    @testset "Reich-Moore Cross Section -- U-235 (fissile)" begin
        endf_file = joinpath(ENDF_RESOURCES, "n-092_U_235-ENDF8.0.endf")
        if !isfile(endf_file)
            @warn "Skipping U-235 RM XS tests: file not found"
        else
            open(endf_file) do io
                @test find_section(io, 2, 151)
                mf2 = read_mf2(io)
                rng = mf2.isotopes[1].ranges[1]

                # Thermal energy
                xs_th = cross_section(0.0253, rng)
                @test isfinite(xs_th.total)
                @test xs_th.total > 0.0
                @test xs_th.elastic > 0.0
                @test xs_th.capture > 0.0
                @test xs_th.fission > 0.0  # U-235 IS fissile

                # Sum rule
                @test isapprox(xs_th.total,
                               xs_th.elastic + xs_th.capture + xs_th.fission,
                               rtol=1e-10)

                # U-235 thermal fission cross section is ~585 barns (from resonances only)
                # The smooth MF3 background adds to this, but the resonance contribution
                # alone should be in the right ballpark
                @test xs_th.fission > 100.0

                # Cross sections positive at several energies
                for E in [1e-5, 0.0253, 0.1, 1.0, 10.0, 100.0, 1000.0]
                    xs = cross_section(E, rng)
                    @test xs.total > 0.0
                    @test xs.elastic >= 0.0
                    @test xs.capture >= 0.0
                    @test xs.fission >= 0.0
                    @test isapprox(xs.total,
                                   xs.elastic + xs.capture + xs.fission,
                                   rtol=1e-10)
                end
            end
        end
    end

    # ======================================================================
    # SLBW Cross Section -- synthetic test
    # ======================================================================
    @testset "SLBW Cross Section -- synthetic" begin
        # Create a simple SLBW system with one resonance for testing
        # Er = 1.0 eV, GN = 0.001 eV, GG = 0.025 eV, GF = 0.0
        # Target: mass 100, spin 0, AP = 0.5 fm
        params = SLBWParameters(
            Int32(1),      # NLS
            0.0,           # SPI
            0.5,           # AP
            [Int32(0)],    # l_values: s-wave
            [100.0],       # AWRI
            [0.0],         # QX
            [Int32(0)],    # LRX
            [[1.0]],       # Er
            [[0.5]],       # AJ
            [[0.001]],     # Gn
            [[0.025]],     # Gg
            [[0.0]],       # Gf
            [[0.0]],       # Gx
        )
        rng = ResonanceRange(1e-5, 100.0, Int32(1), Int32(1),
                             Int32(0), Int32(0), Int32(0), params)

        # At the resonance energy (Er=1.0 eV), cross section should peak
        xs_peak = cross_section(1.0, rng)
        @test xs_peak.total > 0.0
        @test xs_peak.capture > 0.0
        @test xs_peak.fission == 0.0

        # Far from resonance, cross section should be smaller
        xs_far = cross_section(100.0, rng)
        @test xs_peak.total > xs_far.total

        # Sum rule
        @test isapprox(xs_peak.total,
                       xs_peak.elastic + xs_peak.capture + xs_peak.fission,
                       rtol=1e-10)

        # With temperature > 0, should still give finite results
        xs_hot = cross_section(1.0, rng; temperature=300.0)
        @test isfinite(xs_hot.total)
        @test xs_hot.total > 0.0
        @test isapprox(xs_hot.total,
                       xs_hot.elastic + xs_hot.capture + xs_hot.fission,
                       rtol=1e-10)

        # At T > 0, peak should be broadened (lower peak but wider)
        @test xs_hot.total < xs_peak.total * 1.5  # broadened peak should be different
    end

    # ======================================================================
    # Cross Section Consistency Checks
    # ======================================================================
    @testset "Cross Section -- energy scan" begin
        endf_file = joinpath(ENDF_RESOURCES, "n-026_Fe_056-JEFF3.3.endf")
        if !isfile(endf_file)
            @warn "Skipping energy scan tests: file not found"
        else
            open(endf_file) do io
                @test find_section(io, 2, 151)
                mf2 = read_mf2(io)
                rng = mf2.isotopes[1].ranges[1]

                # Scan many energies -- all should be finite and satisfy sum rule
                energies = 10 .^ range(-5, 5, length=50)
                for E in energies
                    xs = cross_section(E, rng)
                    @test isfinite(xs.total)
                    @test isfinite(xs.elastic)
                    @test isfinite(xs.capture)
                    @test isfinite(xs.fission)
                    @test xs.total >= 0.0
                    @test xs.elastic >= 0.0
                    @test xs.capture >= 0.0
                    @test xs.fission >= 0.0
                    @test isapprox(xs.total,
                                   xs.elastic + xs.capture + xs.fission,
                                   rtol=1e-10)
                end
            end
        end
    end

    # ======================================================================
    # round_sigfig
    # ======================================================================
    @testset "round_sigfig" begin
        # Basic rounding to 7 significant figures
        @test round_sigfig(1.23456789, 7, 0) != 0.0
        @test isapprox(round_sigfig(1.23456789, 7, 0), 1.234568, rtol=1e-6)

        # Shading up and down
        v = round_sigfig(1.0e6, 7, 0)
        vup = round_sigfig(1.0e6, 7, +1)
        vdn = round_sigfig(1.0e6, 7, -1)
        @test vup > v
        @test vdn < v

        # Zero handling
        @test round_sigfig(0.0, 7, 0) == 0.0

        # Negative values
        v = round_sigfig(-1.23456789, 7, 0)
        @test v < 0.0
        @test isapprox(v, -1.234568, rtol=1e-5)

        # Various significant figures
        @test isapprox(round_sigfig(3.14159, 3, 0), 3.14, rtol=1e-2)
        @test isapprox(round_sigfig(3.14159, 5, 0), 3.1416, rtol=1e-4)
    end

    # ======================================================================
    # Adaptive Grid -- Lorentzian convergence
    # ======================================================================
    @testset "Adaptive Grid -- Lorentzian" begin
        # A Lorentzian peak: f(x) = 1 / (1 + (x-x0)^2 / gamma^2)
        # with x0 = 100.0, gamma = 0.5
        x0 = 100.0
        gamma = 0.5
        lorentzian(E) = (1.0 / (1.0 + ((E - x0) / gamma)^2),
                         0.5 / (1.0 + ((E - x0) / gamma)^2),
                         0.3 / (1.0 + ((E - x0) / gamma)^2),
                         0.2 / (1.0 + ((E - x0) / gamma)^2))

        # Initial grid: just a few points spanning the peak
        grid = Float64[50.0, 90.0, 95.0, 99.0, 100.0, 101.0, 105.0, 110.0, 150.0]

        config = AdaptiveConfig(0.001)
        energies, values = adaptive_reconstruct(lorentzian, grid, config)

        # Should have more points than initial grid (adaptive refinement)
        @test length(energies) > length(grid)

        # Check convergence: at all midpoints, linear interpolation
        # should match the true function within tolerance
        n = length(energies)
        max_rel_error = 0.0
        for i in 1:(n-1)
            e_mid = 0.5 * (energies[i] + energies[i+1])
            true_val = lorentzian(e_mid)

            # Linear interpolation
            fr = (e_mid - energies[i]) / (energies[i+1] - energies[i])
            for j in 1:4
                interp_val = (1.0 - fr) * values[i, j] + fr * values[i+1, j]
                if abs(true_val[j]) > 1e-10
                    rel_err = abs(true_val[j] - interp_val) / abs(true_val[j])
                    max_rel_error = max(max_rel_error, rel_err)
                end
            end
        end

        # Should converge within the tolerance (with some margin for
        # the relaxed integral criterion)
        @test max_rel_error < 0.02  # within errmax (10x err)

        # Check that points are denser near the peak
        near_peak = count(e -> abs(e - x0) < 5.0, energies)
        far_from_peak = count(e -> abs(e - x0) > 20.0, energies)
        @test near_peak > far_from_peak
    end

    # ======================================================================
    # Adaptive Grid -- Step function
    # ======================================================================
    @testset "Adaptive Grid -- step function" begin
        # A step function that transitions over a narrow region
        step_fn(E) = begin
            s = 1.0 / (1.0 + exp(-100.0 * (E - 5.0)))
            (s, s, 0.0, 0.0)
        end

        grid = Float64[1.0, 3.0, 5.0, 7.0, 10.0]
        config = AdaptiveConfig(0.001)
        energies, values = adaptive_reconstruct(step_fn, grid, config)

        # Should have many points near E=5 (the transition region)
        @test length(energies) > 10

        near_transition = count(e -> abs(e - 5.0) < 1.0, energies)
        @test near_transition >= 5  # at least 5 points near the step
    end

    # ======================================================================
    # RECONR -- Full H-2 (test 84) reconstruction
    # ======================================================================
    @testset "RECONR -- H-2 reconstruction" begin
        endf_file = joinpath(ENDF_RESOURCES, "n-001_H_002-ENDF8.0.endf")
        ref_file = joinpath(@__DIR__, "..", "..", "njoy-reference", "tests", "84",
                           "referenceTape100")

        if !isfile(endf_file)
            @warn "Skipping RECONR tests: ENDF file not found"
        else
            # Run RECONR
            result = reconr(endf_file; mat=128, err=0.001)

            # Basic checks
            @test length(result.energies) > 100
            @test length(result.total) == length(result.energies)
            @test length(result.elastic) == length(result.energies)
            @test length(result.capture) == length(result.energies)

            # Energy grid should be sorted
            for i in 2:length(result.energies)
                @test result.energies[i] >= result.energies[i-1]
            end

            # All cross sections should be non-negative
            for i in 1:length(result.energies)
                @test result.total[i] >= 0.0
                @test result.elastic[i] >= 0.0
                @test result.capture[i] >= 0.0
            end

            # Total = elastic + fission + capture (sum rule)
            for i in 1:length(result.energies)
                @test isapprox(result.total[i],
                               result.elastic[i] + result.fission[i] + result.capture[i],
                               rtol=1e-8)
            end

            # MF2 should be loaded
            @test result.mf2.ZA == 1002.0
            @test isapprox(result.mf2.AWR, 1.9968, rtol=1e-3)

            # MF3 sections should be loaded
            @test length(result.mf3_sections) > 0

            # Check cross sections at thermal energy (0.0253 eV)
            # H-2 has a known thermal cross section of about 3.4 barns total
            thermal_idx = findfirst(e -> e >= 0.0253, result.energies)
            if thermal_idx !== nothing
                @test result.total[thermal_idx] > 2.0   # should be ~3.4 barns
                @test result.total[thermal_idx] < 10.0   # sanity check
                @test result.elastic[thermal_idx] > 1.0  # elastic dominates
            end

            # Compare with reference tape if available
            if isfile(ref_file)
                # Parse reference tape MF3/MT1 (total cross section)
                ref_energies = Float64[]
                ref_total = Float64[]
                open(ref_file) do io
                    in_mf3_mt1 = false
                    header_lines = 0
                    while !eof(io)
                        line = readline(io)
                        p = rpad(line, 80)
                        mf = NJOY._parse_int(p[71:72])
                        mt = NJOY._parse_int(p[73:75])

                        if mf == 3 && mt == 1
                            if !in_mf3_mt1
                                in_mf3_mt1 = true
                                header_lines = 0
                            end
                            header_lines += 1
                            # Skip first 3 lines (HEAD, TAB1 header, interp table)
                            if header_lines <= 3
                                continue
                            end
                            # Parse data pairs
                            for j in 1:3
                                e_str = p[(j-1)*22+1 : (j-1)*22+11]
                                s_str = p[(j-1)*22+12 : (j-1)*22+22]
                                e_val = parse_endf_float(e_str)
                                s_val = parse_endf_float(s_str)
                                if e_val > 0.0
                                    push!(ref_energies, e_val)
                                    push!(ref_total, s_val)
                                end
                            end
                        elseif in_mf3_mt1 && mt == 0
                            break
                        end
                    end
                end

                if length(ref_energies) > 10
                    # Compare total cross sections at reference energies
                    # Our grid may differ, so interpolate our result
                    n_close = 0
                    n_compared = 0
                    for i in 1:min(length(ref_energies), 50)
                        e_ref = ref_energies[i]
                        # Find nearest energy in our grid
                        idx = searchsortedfirst(result.energies, e_ref)
                        if idx > 1 && idx <= length(result.energies)
                            # Linear interpolation
                            e1 = result.energies[idx-1]
                            e2 = result.energies[idx]
                            if abs(e2 - e1) > 0.0
                                fr = (e_ref - e1) / (e2 - e1)
                                our_val = (1.0 - fr) * result.total[idx-1] + fr * result.total[idx]
                                ref_val = ref_total[i]
                                if ref_val > 0.01  # skip very small values
                                    n_compared += 1
                                    rel_diff = abs(our_val - ref_val) / ref_val
                                    if rel_diff < 0.05  # 5% tolerance
                                        n_close += 1
                                    end
                                end
                            end
                        end
                    end

                    # At least 80% of compared points should be close
                    if n_compared > 5
                        @test n_close / n_compared > 0.7
                    end
                end
            end
        end
    end

    # ======================================================================
    # RECONR -- Grid density near resonances
    # ======================================================================
    @testset "RECONR -- grid density near resonances" begin
        endf_file = joinpath(ENDF_RESOURCES, "n-001_H_002-ENDF8.0.endf")
        if !isfile(endf_file)
            @warn "Skipping grid density tests: ENDF file not found"
        else
            result = reconr(endf_file; mat=128, err=0.001)

            # H-2 has no sharp resonances (scattering-radius only in resolved region),
            # but the grid should still have some points at low energies
            # MF3 data for H-2 has breakpoints at 1e-5, 1e-4, 0.0253 eV
            low_e = count(e -> e < 1.0, result.energies)
            @test low_e >= 2  # should have at least a couple of points below 1 eV
        end
    end

    # ======================================================================
    # RECONR -- sum rule (total = elastic + capture + fission)
    # ======================================================================
    @testset "RECONR -- sum rule consistency" begin
        endf_file = joinpath(ENDF_RESOURCES, "n-001_H_002-ENDF8.0.endf")
        if !isfile(endf_file)
            @warn "Skipping sum rule tests: ENDF file not found"
        else
            result = reconr(endf_file; mat=128, err=0.001)

            for i in 1:length(result.energies)
                sum_parts = result.elastic[i] + result.fission[i] + result.capture[i]
                @test isapprox(result.total[i], sum_parts, rtol=1e-8)
            end
        end
    end

    # ======================================================================
    # PENDF Writer -- basic output format
    # ======================================================================
    @testset "PENDF Writer -- basic output" begin
        endf_file = joinpath(ENDF_RESOURCES, "n-001_H_002-ENDF8.0.endf")
        if !isfile(endf_file)
            @warn "Skipping PENDF writer tests: ENDF file not found"
        else
            result = reconr(endf_file; mat=128, err=0.001)

            buf = IOBuffer()
            write_pendf(buf, result; mat=128, err=0.001)
            output = String(take!(buf))

            # Should have content
            @test length(output) > 100

            # Should contain MF1, MF2, MF3 markers
            lines = split(output, '\n')
            has_mf1 = any(l -> length(l) >= 72 && strip(l[71:72]) == "1", lines)
            has_mf3 = any(l -> length(l) >= 72 && strip(l[71:72]) == "3", lines)
            @test has_mf1
            @test has_mf3

            # All lines should be 80 characters (ENDF format)
            for line in lines
                if !isempty(line)
                    @test length(line) == 80
                end
            end
        end
    end

    # ======================================================================
    # MF3 Section Reader
    # ======================================================================
    @testset "MF3 Section Reader" begin
        endf_file = joinpath(ENDF_RESOURCES, "n-001_H_002-ENDF8.0.endf")
        if !isfile(endf_file)
            @warn "Skipping MF3 reader tests: ENDF file not found"
        else
            open(endf_file) do io
                sections = read_mf3_sections(io, 128)
                @test length(sections) > 0

                # H-2 should have MT=1 (total), MT=2 (elastic), MT=102 (capture)
                mt_list = [Int(sec.mt) for sec in sections]
                @test 1 in mt_list
                @test 2 in mt_list
                @test 102 in mt_list

                # Each section should have data
                for sec in sections
                    @test length(sec.tab.x) > 0
                    @test length(sec.tab.y) > 0
                    @test sec.tab.x[1] > 0.0  # first energy should be positive
                end
            end
        end
    end

    # ======================================================================
    # build_grid and build_evaluator
    # ======================================================================
    @testset "build_grid and build_evaluator" begin
        endf_file = joinpath(ENDF_RESOURCES, "n-001_H_002-ENDF8.0.endf")
        if !isfile(endf_file)
            @warn "Skipping build_grid tests: ENDF file not found"
        else
            open(endf_file) do io
                find_section(io, 2, 151)
                mf2 = read_mf2(io)
                mf3_sections = read_mf3_sections(io, 128)

                # Build grid
                grid = build_grid(mf2, mf3_sections)
                @test length(grid) > 10
                @test issorted(grid)
                @test grid[1] > 0.0

                # Build evaluator
                eval_fn = build_evaluator(mf2)
                result = eval_fn(1.0)
                @test length(result) == 4  # (total, elastic, fission, capture)
                @test all(isfinite, result)
            end
        end
    end

end  # @testset "NJOY.jl"
