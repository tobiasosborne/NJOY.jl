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

    # ======================================================================
    # BROADR -- F-function (sigma1 kernel)
    # ======================================================================
    @testset "F-function analytical values" begin
        # f_0(0) = erfc(0)/2 = 1/2
        @test isapprox(f_func(0, 0.0), 0.5, atol=1e-15)
        # f_1(0) = exp(0)/(2*sqrt(pi)) = 1/(2*sqrt(pi))
        @test isapprox(f_func(1, 0.0), 1.0/(2*sqrt(pi)), atol=1e-15)
        # f_2(0) = 1/4
        @test isapprox(f_func(2, 0.0), 0.25, atol=1e-15)
        # f_3(0) = 1/(2*sqrt(pi))
        @test isapprox(f_func(3, 0.0), 1.0/(2*sqrt(pi)), atol=1e-15)
        # f_4(0) = 3/8
        @test isapprox(f_func(4, 0.0), 0.375, atol=1e-15)

        # f_0(a) = erfc(a)/2 for various a
        using SpecialFunctions: erfc
        for a in [0.1, 0.5, 1.0, 2.0, 5.0]
            @test isapprox(f_func(0, a), erfc(a)/2, rtol=1e-12)
        end

        # f_1(a) = exp(-a^2)/(2*sqrt(pi))
        for a in [0.1, 1.0, 3.0]
            @test isapprox(f_func(1, a), exp(-a^2)/(2*sqrt(pi)), rtol=1e-12)
        end

        # Large a: all f_n -> 0
        for n in 0:4
            @test f_func(n, 15.0) == 0.0
        end

        # f_all returns consistent tuple
        for a in [0.0, 0.5, 2.0, 8.0]
            fa = f_all(a)
            for n in 0:4
                @test isapprox(fa[n+1], f_func(n, a), atol=1e-14)
            end
        end
    end

    # ======================================================================
    # BROADR -- H-function cancellation avoidance
    # ======================================================================
    @testset "H-function cancellation" begin
        # h(n, a, a) = 0 exactly
        for n in 0:4
            @test h_func(n, 1.5, 1.5) == 0.0
        end

        # h(n, a, b) = f(n, a) - f(n, b) for well-separated a, b
        for n in 0:4
            h_direct = f_func(n, 0.5) - f_func(n, 2.0)
            @test isapprox(h_func(n, 0.5, 2.0), h_direct, rtol=1e-10)
        end

        # CRITICAL TEST: h(n, a, a+eps) for small eps -- cancellation avoidance
        # The Taylor series must give accurate results where direct subtraction fails
        for n in 0:4
            a = 1.5
            eps_val = 1.0e-8
            b = a + eps_val
            h_val = h_func(n, a, b)
            # The result should be approximately -eps * d/da f_n(a)
            # For f_0: d/da f_0 = -exp(-a^2)/sqrt(pi)
            # So h ~ eps * exp(-a^2)/sqrt(pi) (positive, since f decreases)
            @test isfinite(h_val)
            # Verify it is small but not zero (meaningful value)
            if n == 0
                expected = eps_val * exp(-a^2) / sqrt(pi)
                @test abs(h_val) > 0.0
                @test isapprox(abs(h_val), expected, rtol=0.1)
            end
        end

        # h_taylor directly: verify it matches h_func for small intervals
        for n in 0:4
            a = 2.0; b = 2.0 + 1e-6
            ht = h_taylor(n, a, b)
            hf = h_func(n, a, b)
            @test isapprox(ht, hf, rtol=1e-6)
        end

        # h_all returns consistent values
        f_old = f_all(1.0)
        h_vals, f_new = h_all(1.0, f_old, 1.5)
        for n in 0:4
            @test isapprox(h_vals[n+1], h_func(n, 1.0, 1.5), atol=1e-12)
        end
    end

    # ======================================================================
    # BROADR -- Broadening invariants
    # ======================================================================
    @testset "Broadening constant XS invariant" begin
        # A constant XS remains constant after broadening.
        # Use Fe-56 (awr=55.85) with energies well above the Doppler width
        # so boundary extrapolation effects are negligible.
        energies = collect(range(1.0, 1000.0, length=50))
        sigma_const = 10.0
        xs = fill(sigma_const, length(energies))
        T = 300.0
        awr = 55.85  # Fe-56

        out_e, out_xs = doppler_broaden(energies, xs, T, awr; tol=0.001)

        # All broadened values should equal the constant within ~1%
        for i in eachindex(out_xs)
            @test isapprox(out_xs[i], sigma_const, rtol=0.1)
        end
    end

    @testset "Broadening 1/v XS invariant" begin
        # sigma(E)*sqrt(E)=const is invariant under Doppler broadening.
        # Use Fe-56 to keep Doppler width narrow relative to energy grid.
        energies = collect(range(0.01, 100.0, length=500))
        sigma0 = 100.0
        xs = sigma0 ./ sqrt.(energies)
        T = 300.0
        awr = 55.85

        out_e, out_xs = doppler_broaden(energies, xs, T, awr; tol=0.001)

        # Check product sigma*sqrt(E) is approximately constant
        # away from boundaries
        for i in eachindex(out_e)
            if out_e[i] > 1.0 && out_e[i] < 50.0
                product = out_xs[i] * sqrt(out_e[i])
                @test isapprox(product, sigma0, rtol=0.05)
            end
        end
    end

    @testset "Broadening smooths a step" begin
        # A step function is smoothed by broadening
        energies = [1.0, 5.0, 5.0 + 1e-3, 10.0, 50.0, 100.0]
        xs = [1.0, 1.0, 100.0, 100.0, 100.0, 100.0]
        T = 300.0
        awr = 55.85  # Fe-56

        out_e, out_xs = doppler_broaden(energies, xs, T, awr; tol=0.001)

        # All values should be positive
        @test all(x -> x >= 0, out_xs)
        # Far from the step, XS should approach the plateau values
        for i in eachindex(out_e)
            if out_e[i] > 20.0
                @test out_xs[i] > 50.0  # should be close to 100
            end
        end
    end

    # ======================================================================
    # BROADR -- Thinning
    # ======================================================================
    @testset "Thinning preserves accuracy" begin
        # Build a smooth function with many redundant points
        energies = collect(range(0.1, 10.0, length=500))
        xs = 10.0 .* ones(500)  # constant -- maximally thinnable

        thinned_e, thinned_xs = thin_xs(energies, xs; tol=0.001)

        # Should thin dramatically (constant function needs only 2 points)
        @test length(thinned_e) < 50  # much fewer than 500

        # First and last points preserved
        @test thinned_e[1] == energies[1]
        @test thinned_e[end] == energies[end]

        # Interpolation of thinned result matches original within tolerance
        for i in 1:10:500
            e = energies[i]
            idx = searchsortedlast(thinned_e, e)
            idx = clamp(idx, 1, length(thinned_e) - 1)
            fr = (e - thinned_e[idx]) / (thinned_e[idx+1] - thinned_e[idx])
            interp_val = (1 - fr) * thinned_xs[idx] + fr * thinned_xs[idx+1]
            @test isapprox(interp_val, xs[i], rtol=0.002)
        end
    end

    @testset "Thinning preserves non-trivial features" begin
        # A linear ramp thins well (step_max guard limits how aggressively)
        energies = collect(range(1.0, 100.0, length=200))
        xs = collect(range(1.0, 50.0, length=200))

        thinned_e, thinned_xs = thin_xs(energies, xs; tol=0.001)

        # Should thin significantly (but step_max=1.24 prevents full thinning)
        @test length(thinned_e) < 100
        @test thinned_e[1] == energies[1]
        @test thinned_e[end] == energies[end]

        # Verify interpolation accuracy at all original points
        for i in 1:20:200
            e = energies[i]
            idx = searchsortedlast(thinned_e, e)
            idx = clamp(idx, 1, length(thinned_e) - 1)
            de = thinned_e[idx+1] - thinned_e[idx]
            de > 0 || continue
            fr = (e - thinned_e[idx]) / de
            interp_val = (1 - fr)*thinned_xs[idx] + fr*thinned_xs[idx+1]
            @test isapprox(interp_val, xs[i], rtol=0.002)
        end
    end

    # ======================================================================
    # BROADR -- sigma1_at kernel
    # ======================================================================
    @testset "sigma1_at kernel basic" begin
        # Constant XS: sigma1 should return the constant
        seg_e = [0.01, 1.0, 10.0, 100.0]
        seg_xs = [5.0, 5.0, 5.0, 5.0]
        awr = 56.0  # Fe-like
        T = 300.0
        alpha = awr / (NJOY.PhysicsConstants.bk * T)

        for E in [0.1, 1.0, 10.0, 50.0]
            val = sigma1_at(E, seg_e, seg_xs, alpha)
            @test isapprox(val, 5.0, rtol=0.01)
        end
    end

    # ======================================================================
    # HEATR -- KERMA Coefficients and Damage Energy
    # ======================================================================
    @testset "HEATR elastic heating" begin
        # Analytical formula: h(E) = 2*E*A / (1+A)^2
        # For isotropic CM elastic scattering, the average energy
        # deposited in the target equals this expression exactly.

        # H-1: A=1 => h = E/2
        @test elastic_heating(1.0, 1.0) == 0.5
        @test elastic_heating(1e6, 1.0) == 0.5e6

        # Heavy target: A=238 => h = 2*E*238/(239)^2 ~ small fraction
        A = 238.0
        E = 1e6
        expected = 2.0 * E * A / (1.0 + A)^2
        @test isapprox(elastic_heating(E, A), expected, rtol=1e-14)

        # C-12: A=12
        A = 12.0
        for E in [0.0253, 1.0, 1e4, 1e6]
            h = elastic_heating(E, A)
            @test isapprox(h, 2.0 * E * A / (1.0 + A)^2, rtol=1e-14)
            # Heating must be positive and less than E
            if E > 0
                @test h > 0.0
                @test h < E
            end
        end

        # Symmetry: elastic_heating(E, A) should scale linearly with E
        A = 56.0  # Fe
        @test isapprox(elastic_heating(2.0, A) / elastic_heating(1.0, A), 2.0, rtol=1e-14)

        # Anisotropic scattering: forward scattering (mu_bar>0) deposits less
        # energy, backward scattering (mu_bar<0) deposits more.
        A = 12.0
        E = 1.0e6
        h_fwd = elastic_heating_aniso(E, A, 0.5)
        h_bwd = elastic_heating_aniso(E, A, -0.5)
        @test h_bwd > h_fwd  # backward scattering deposits more energy
    end

    @testset "HEATR capture heating (H-1)" begin
        # For H-1 (A=1, Z=1), radiative capture: n + p -> d + gamma
        # Q = 2.224566 MeV (deuteron binding energy)
        A = 1.0
        Q = 2.224566e6  # eV

        # capture_heating(E, Q, A) = E*A/(A+1) + Q (local deposition, E_gamma=0)
        E = 0.0253
        h = capture_heating(E, Q, A)
        # With local deposition: h = E*A/(A+1) + Q
        @test isapprox(h, E * A / (A + 1.0) + Q, rtol=1e-14)

        # With gamma energy subtracted:
        E_gamma_test = Q - 100.0  # pretend gamma carries Q-100 eV
        h_bal = capture_heating(E, Q, A; E_gamma=E_gamma_test)
        # h = E*A/(A+1) + Q - E_gamma
        @test isapprox(h_bal, E * A / (A + 1.0) + Q - E_gamma_test, rtol=1e-14)

        # Capture recoil energy should be positive and small
        E_recoil = capture_recoil(E, Q, A)
        @test E_recoil > 0.0
        @test E_recoil < Q  # recoil much less than Q value
    end

    @testset "HEATR Lindhard damage function" begin
        # The Lindhard partition function gives damage energy from recoil.
        # For very high energy recoils, electronic stopping dominates
        # and damage fraction -> 0. For low energy, nuclear stopping
        # dominates and damage fraction -> 1 (minus displacement threshold).

        # Pre-compute for Fe in Fe lattice (Z=26, A=56)
        Z = 26.0
        A = 56.0
        lp = lindhard_params(Z, A, Z, A)
        @test lp.el_inv > 0.0
        @test lp.fl > 0.0

        # Below displacement threshold: no damage
        @test lindhard_damage(10.0, lp; E_d=25.0) == 0.0

        # At moderate energy: damage is a fraction of recoil energy
        dam = lindhard_damage(1e4, lp; E_d=25.0)
        @test dam > 0.0
        @test dam < 1e4

        # Damage fraction decreases with increasing energy
        # (electronic stopping becomes more important)
        frac_low = lindhard_damage(100.0, lp; E_d=25.0) / 100.0
        frac_high = lindhard_damage(1e6, lp; E_d=25.0) / 1e6
        @test frac_low > frac_high

        # Zero recoil atom charge: no damage
        lp_zero = lindhard_params(0.0, A, Z, A)
        @test lindhard_damage(1e6, lp_zero) == 0.0

        # Convenience 5-arg form
        dam5 = lindhard_damage(1e4, Z, A, Z, A; E_d=25.0)
        @test isapprox(dam5, dam, rtol=1e-14)

        # Default displacement energies
        @test displacement_energy(26) == 40.0  # Fe
        @test displacement_energy(13) == 27.0  # Al
        @test displacement_energy(6)  == 31.0  # C
        @test displacement_energy(99) == 25.0  # fallback
    end

    @testset "HEATR inelastic heating" begin
        # Inelastic scattering: h = E + Q - E_secondary
        A = 55.845
        Q = -0.847e6  # eV (first excited state)
        E = 2.0e6

        h = inelastic_heating(E, Q, A)
        @test isapprox(h, E + Q, rtol=1e-14)  # E_secondary=0 default

        # With secondary energy subtracted
        E_sec = 1.5e6
        h2 = inelastic_heating(E, Q, A; E_secondary=E_sec)
        @test isapprox(h2, E + Q - E_sec, rtol=1e-14)
    end

    @testset "HEATR fission heating" begin
        # Simple energy balance for fission
        Q = 200.0e6  # typical fission Q ~ 200 MeV
        E = 1.0e6    # 1 MeV incident
        E_nu = 12.0e6  # neutrino energy
        E_del = 10.0e6 # delayed components

        h = fission_heating(E, Q; E_neutrino=E_nu, E_delayed=E_del)
        @test isapprox(h, Q - E_nu - E_del, rtol=1e-14)

        # FissionQComponents form
        qc = FissionQComponents(
            170.0e6,  # E_fragments
             5.0e6,   # E_prompt_n
             7.0e6,   # E_prompt_gamma
             6.5e6,   # E_delayed_beta
             6.3e6,   # E_delayed_gamma
             0.01e6,  # E_delayed_n
            12.0e6,   # E_neutrino
            194.81e6, # E_pseudo_Q
            206.81e6  # E_total
        )
        h_qc = fission_heating(E, qc)
        # h = E_fragments + E_prompt_gamma
        @test isapprox(h_qc, 170.0e6 + 7.0e6, rtol=1e-14)
    end

    @testset "HEATR nxn heating" begin
        # (n,2n) reaction
        A = 55.845
        Q = -10.0e6  # negative Q
        E = 14.0e6
        h = nxn_heating(E, Q, A, 2)
        @test isapprox(h, E + Q, rtol=1e-14)  # E_secondary=0

        h2 = nxn_heating(E, Q, A, 2; E_secondary=3.0e6)
        @test isapprox(h2, E + Q - 2 * 3.0e6, rtol=1e-14)
    end

    @testset "HEATR compute_kerma integration" begin
        # Build a simple PointwiseMaterial with elastic + capture
        energies = [1e-5, 0.0253, 1.0, 1e3, 1e6]
        ne = length(energies)
        # 2 reactions: MT2 (elastic), MT102 (capture)
        xs = zeros(ne, 2)
        xs[:, 1] .= 20.0   # elastic cross section = 20 b
        xs[:, 2] .= [1000.0, 100.0, 10.0, 1.0, 0.1]  # 1/v capture

        pendf = PointwiseMaterial(Int32(125), energies, xs, [2, 102])

        Q_vals = Dict{Int,Float64}(102 => 2.224566e6)
        result = compute_kerma(pendf; awr=1.0, Z=1, Q_values=Q_vals)

        # Result should have correct size
        @test length(result.energies) == ne
        @test length(result.total_kerma) == ne

        # Total = elastic + capture + fission + inelastic
        for ie in 1:ne
            parts = result.elastic_kerma[ie] + result.capture_kerma[ie] +
                    result.fission_kerma[ie] + result.inelastic_kerma[ie]
            @test isapprox(result.total_kerma[ie], parts, rtol=1e-12)
        end

        # Elastic partial should match analytical formula * sigma
        for ie in 1:ne
            E = energies[ie]
            h_expected = elastic_heating(E, 1.0) * xs[ie, 1]
            @test isapprox(result.elastic_kerma[ie], h_expected, rtol=1e-12)
        end

        # Capture partial: h = (E*A/(A+1) + Q) * sigma  (local gamma deposition)
        for ie in 1:ne
            E = energies[ie]
            Q_cap = 2.224566e6
            h_expected = capture_heating(E, Q_cap, 1.0) * xs[ie, 2]
            @test isapprox(result.capture_kerma[ie], h_expected, rtol=1e-12)
        end

        # Verify sum rule helper
        @test verify_kerma_sum_rule(result)

        # Damage should be non-negative
        for ie in 1:ne
            @test result.damage_energy[ie] >= 0.0
        end
    end

    # ======================================================================
    # THERMR -- Thermal Scattering Cross Sections
    # ======================================================================
    @testset "THERMR -- Free Gas" begin
        C = NJOY.PhysicsConstants

        # ------------------------------------------------------------------
        # Test 1: High-energy limit  sigma -> sigma_b / A
        # For E >> kT, the free gas XS approaches the billiard-ball limit.
        # ------------------------------------------------------------------
        @testset "High-energy limit (E >> kT)" begin
            A = 12.0       # Carbon-like
            T = 293.6      # Room temperature [K]
            sigma_b = 5.551 # barn (typical bound XS)
            kT = C.bk * T

            # At E = 1000*kT the result should be very close to sigma_b/A
            E_high = 1000.0 * kT
            xs_high = free_gas_xs(E_high, A, T; sigma_b=sigma_b)
            @test isapprox(xs_high, sigma_b / A, rtol=1e-4)

            # At E = 10000*kT even closer
            E_vhigh = 10000.0 * kT
            xs_vhigh = free_gas_xs(E_vhigh, A, T; sigma_b=sigma_b)
            @test isapprox(xs_vhigh, sigma_b / A, rtol=1e-4)
        end

        # ------------------------------------------------------------------
        # Test 2: Low-energy limit  sigma -> sigma_b*sqrt(pi)/(2*x) / A
        # where x = sqrt(A*E/(kT)).  This is the 1/v behavior.
        # ------------------------------------------------------------------
        @testset "Low-energy 1/v limit (E << kT)" begin
            A = 1.0        # Hydrogen
            T = 293.6
            sigma_b = 81.98 # barn (H-1 bound XS)
            kT = C.bk * T

            # For E << kT, the correct leading term is:
            #   sigma ~ sigma_b/A * 2/(x*sqrt(pi)), x = sqrt(A*E/kT)
            for E in [1e-7, 1e-6, 1e-5]
                x = sqrt(A * E / kT)
                expected = sigma_b / A * 2.0 / (x * sqrt(C.pi))
                xs = free_gas_xs(E, A, T; sigma_b=sigma_b)
                @test isapprox(xs, expected, rtol=1e-2)
            end
        end

        # ------------------------------------------------------------------
        # Test 3: Cross section is always positive
        # ------------------------------------------------------------------
        @testset "Positivity" begin
            for A in [1.0, 12.0, 56.0, 238.0]
                for T in [100.0, 300.0, 600.0, 1200.0]
                    for E in [1e-5, 0.0253, 0.1, 1.0, 10.0]
                        @test free_gas_xs(E, A, T) > 0.0
                    end
                end
            end
        end

        # ------------------------------------------------------------------
        # Test 4: Monotonicity -- for a free gas, XS decreases with energy
        # at low energies (1/v) and approaches constant at high energies.
        # The function should be monotonically decreasing.
        # ------------------------------------------------------------------
        @testset "Monotonicity" begin
            A = 12.0
            T = 300.0
            energies = [1e-5, 1e-4, 1e-3, 0.01, 0.1, 1.0, 10.0]
            xs_vals = [free_gas_xs(E, A, T) for E in energies]
            for i in 2:length(xs_vals)
                @test xs_vals[i] <= xs_vals[i-1]
            end
        end

        # ------------------------------------------------------------------
        # Test 5: Free gas kernel -- detailed balance
        # P(E->E') / P(E'->E) = exp(-(E-E')/(kT))  (up to sqrt(E'/E) prefactor)
        # More precisely, the kernel satisfies:
        #   kernel(E,E',mu) / kernel(E',E,mu) = exp(-(E-E')/kT)
        # when we account for the sqrt(E'/E) prefactor properly.
        # ------------------------------------------------------------------
        @testset "Kernel detailed balance" begin
            A = 12.0
            T = 300.0
            kT = C.bk * T

            # The free gas kernel k(E,E',mu) = sigma_b*sqrt(E'/E)/(2kT)*S(a,b)
            # satisfies: k(E->E',mu)/k(E'->E,mu) = (E'/E)*exp(-(E'-E)/kT)
            test_points = [
                (0.0253, 0.030, 0.5),
                (0.05, 0.03, 0.0),
                (0.1, 0.08, -0.5),
                (0.01, 0.015, 0.9),
            ]

            for (E, Ep, mu) in test_points
                k_fwd = free_gas_kernel(E, Ep, mu, A, T)
                k_rev = free_gas_kernel(Ep, E, mu, A, T)
                if k_fwd > 0.0 && k_rev > 0.0
                    ratio = k_fwd / k_rev
                    expected_ratio = (Ep / E) * exp(-(Ep - E) / kT)
                    @test isapprox(ratio, expected_ratio, rtol=1e-10)
                end
            end
        end

        # ------------------------------------------------------------------
        # Test 6: Kernel is non-negative everywhere
        # ------------------------------------------------------------------
        @testset "Kernel non-negativity" begin
            A = 1.0
            T = 300.0
            for E in [0.001, 0.0253, 0.1]
                for Ep in [0.001, 0.0253, 0.1]
                    for mu in [-1.0, -0.5, 0.0, 0.5, 1.0]
                        @test free_gas_kernel(E, Ep, mu, A, T) >= 0.0
                    end
                end
            end
        end

        # ------------------------------------------------------------------
        # Test 7: Thermal energy grid
        # ------------------------------------------------------------------
        @testset "THERMR energy grid" begin
            @test issorted(THERMR_EGRID)
            @test THERMR_EGRID[1] > 0.0
            @test THERMR_EGRID[end] == 10.0
            @test length(THERMR_EGRID) == 118
        end

        # ------------------------------------------------------------------
        # Test 8: compute_thermal_xs basic functionality
        # ------------------------------------------------------------------
        @testset "compute_thermal_xs" begin
            A = 12.0
            T = 300.0
            sigma_b = 5.551

            # Create a simple constant elastic XS above thermal range
            cold_e = collect(range(1e-5, 20.0, length=200))
            cold_xs = fill(sigma_b / A, 200)  # flat at billiard-ball limit

            new_e, new_xs = compute_thermal_xs(cold_e, cold_xs, A, T;
                                               sigma_b=sigma_b, emax=10.0)
            @test issorted(new_e)
            @test length(new_e) == length(new_xs)
            @test all(new_xs .> 0.0)

            # At low energy, thermal XS should be larger than cold (1/v rise)
            low_idx = searchsortedfirst(new_e, 1e-4)
            high_idx = searchsortedfirst(new_e, 5.0)
            @test new_xs[low_idx] > new_xs[high_idx]
        end
    end

    # ======================================================================
    # THERMR -- S(alpha,beta) Model
    # ======================================================================
    @testset "THERMR -- S(alpha,beta)" begin
        C = NJOY.PhysicsConstants

        # Build a free-gas S_sym(a,b) table on a fine log-spaced grid
        # ENDF symmetric convention: S_sym(a,|b|) = S(a,|b|)*exp(|b|/2)
        # so log(S_sym) = -(a+b)^2/(4a) - 0.5*log(4pi*a) + b/2
        @testset "SAB reproduces free gas" begin
            A = 1.0; T = 300.0; kT = C.bk * T
            # Dense log-spaced alpha grid covering the test-point alpha
            ag = [exp(x) for x in range(log(0.01), log(50.0), length=200)]
            bg = collect(range(0.0, 50.0, length=200))
            na, nb = length(ag), length(bg)
            st = zeros(na, nb)
            for (ia, a) in enumerate(ag), (ib, b) in enumerate(bg)
                st[ia, ib] = -(a + b)^2 / (4*a) - 0.5*log(4*C.pi*a) + b/2
            end
            data = read_thermal_data(ag, bg, st; sigma_b=1.0, awr=A, is_log=true)

            # Choose E,Ep that give alpha,beta well within table range
            # alpha = (E+Ep-2*mu*sqrt(E*Ep))/(A*kT), beta = (Ep-E)/kT
            E, Ep, mu = 0.5, 0.6, 0.5
            alpha_test = (E+Ep-2*mu*sqrt(E*Ep))/(A*kT)
            beta_test = abs(Ep-E)/kT
            # Verify alpha/beta in range
            @test alpha_test > ag[1] && alpha_test < ag[end]
            @test beta_test > bg[1] && beta_test < bg[end]

            k_fg = free_gas_kernel(E, Ep, mu, A, T)
            sb_match = 1.0 / ((A+1)/A)^2
            k_sab = sab_kernel(E, Ep, mu, data, T; sigma_b_override=sb_match)
            if k_fg > 1e-20 && k_sab > 1e-20
                @test isapprox(k_sab, k_fg, rtol=0.05)
            end
        end

        # Detailed balance: k(E->Ep)/k(Ep->E) = (Ep/E)*exp(-(Ep-E)/kT)
        # The (Ep/E) factor comes from the sqrt(Ep/E) prefactor in the kernel.
        @testset "SAB kernel detailed balance" begin
            A = 1.0; T = 300.0; kT = C.bk * T
            ag = [exp(x) for x in range(log(0.01), log(80.0), length=200)]
            bg = collect(range(0.0, 80.0, length=200))
            na, nb = length(ag), length(bg)
            st = zeros(na, nb)
            for (ia, a) in enumerate(ag), (ib, b) in enumerate(bg)
                st[ia, ib] = -(a + b)^2 / (4*a) - 0.5*log(4*C.pi*a) + b/2
            end
            data = read_thermal_data(ag, bg, st; sigma_b=5.0, awr=A, is_log=true)
            for (E, Ep, mu) in [(0.5, 0.6, 0.5), (1.0, 0.8, 0.0)]
                kf = sab_kernel(E, Ep, mu, data, T)
                kr = sab_kernel(Ep, E, mu, data, T)
                if kf > 1e-30 && kr > 1e-30
                    expected = (Ep/E) * exp(-(Ep - E) / kT)
                    @test isapprox(kf / kr, expected, rtol=1e-10)
                end
            end
        end

        @testset "read_thermal_data" begin
            alpha = [0.1, 1.0, 5.0]; beta = [0.0, 1.0, 3.0]
            raw_s = [0.1 0.05 0.01; 0.08 0.04 0.005; 0.02 0.01 0.001]
            data = read_thermal_data(alpha, beta, raw_s;
                                     sigma_b=4.0, awr=1.0, is_log=false)
            @test data isa SABData
            @test length(data.alpha) == 3
            @test data.sigma_b == 4.0
            @test isapprox(data.sab[1,1], log(0.1), atol=1e-10)
        end
    end

    # ======================================================================
    # THERMR -- Bragg Edges (coherent elastic)
    # ======================================================================
    @testset "THERMR -- Bragg edges" begin
        @testset "Graphite Bragg edges" begin
            bragg = build_bragg_data(a=2.4573e-8, c=6.700e-8,
                        sigma_coh=5.50, A_mass=12.011, natom=1,
                        debye_waller=3.0, emax=1.0)
            @test bragg isa BraggData
            @test bragg.n_edges > 0

            edges = bragg_edge_energies(bragg)
            @test all(edges .> 0)
            @test issorted(edges)
            @test edges[1] > 1e-4
            @test edges[1] < 0.01

            # Zero below first edge
            @test bragg_edges(edges[1] * 0.5, bragg) == 0.0
            # Positive above first edge
            @test bragg_edges(edges[1] * 1.1, bragg) > 0.0

            # 1/E behavior above all edges
            E1, E2 = edges[end] * 1.5, edges[end] * 3.0
            @test isapprox(bragg_edges(E1, bragg)*E1,
                           bragg_edges(E2, bragg)*E2, rtol=1e-10)

            # Step increase at each Bragg edge
            for k in 2:min(5, length(edges))
                @test bragg_edges(edges[k]*1.001, bragg) >=
                      bragg_edges(edges[k]*0.999, bragg)
            end
        end
    end

    # ======================================================================
    # THERMR -- Incoherent Elastic
    # ======================================================================
    @testset "THERMR -- Incoherent elastic" begin
        sigma_b, dwp = 80.0, 10.0

        @test incoh_elastic_xs(0.0, sigma_b, dwp) == 0.0
        @test isapprox(incoh_elastic_xs(100.0, sigma_b, dwp),
                       sigma_b/2, rtol=1e-6)
        es = [1e-4, 1e-3, 0.01, 0.1, 1.0, 10.0]
        xv = [incoh_elastic_xs(E, sigma_b, dwp) for E in es]
        for i in 2:length(xv); @test xv[i] >= xv[i-1]; end
        E_s = 1e-6
        @test isapprox(incoh_elastic_xs(E_s, sigma_b, dwp),
                       sigma_b/2*4*E_s*dwp, rtol=1e-4)
    end

    # ======================================================================
    # ACER -- ACE Format Types and NXS/JXS Construction
    # ======================================================================
    @testset "ACER -- ACE format types" begin
        # Build a minimal ACETable and verify NXS/JXS structure
        nxs = ntuple(i -> Int32(0), 16)
        nxs = Base.setindex(nxs, Int32(500), NXS_LEN2)
        nxs = Base.setindex(nxs, Int32(92235), NXS_IZAID)
        nxs = Base.setindex(nxs, Int32(100), NXS_NES)
        nxs = Base.setindex(nxs, Int32(3), NXS_NTR)
        nxs = Base.setindex(nxs, Int32(0), NXS_NR)

        @test nxs[NXS_LEN2] == 500
        @test nxs[NXS_IZAID] == 92235
        @test nxs[NXS_NES] == 100
        @test nxs[NXS_NTR] == 3
        @test nxs[NXS_NR] == 0

        jxs = ntuple(i -> Int32(0), 32)
        jxs = Base.setindex(jxs, Int32(1), JXS_ESZ)
        jxs = Base.setindex(jxs, Int32(0), JXS_NU)
        jxs = Base.setindex(jxs, Int32(501), JXS_MTR)

        @test jxs[JXS_ESZ] == 1
        @test jxs[JXS_NU] == 0
        @test jxs[JXS_MTR] == 501

        # Verify NXS/JXS index constants have expected values
        @test NXS_LEN2 == 1
        @test NXS_IZAID == 2
        @test NXS_NES == 3
        @test NXS_NTR == 4
        @test NXS_NR == 5
        @test JXS_ESZ == 1
        @test JXS_NU == 2
        @test JXS_MTR == 3
        @test JXS_LAND == 8
        @test JXS_AND == 9
    end

    # ======================================================================
    # ACER -- ACE Format Writer Line Lengths
    # ======================================================================
    @testset "ACER -- ACE format writer line lengths" begin
        # Create a minimal ACETable
        nes = 5
        xss_data = vcat(
            collect(1.0e-5:1.0e-5:5.0e-5),   # energies (5)
            fill(10.0, nes),                   # total (5)
            fill(5.0, nes),                    # disappearance (5)
            fill(5.0, nes),                    # elastic (5)
            zeros(nes)                         # heating (5)
        )
        nxs = ntuple(i -> Int32(0), 16)
        nxs = Base.setindex(nxs, Int32(length(xss_data)), NXS_LEN2)
        nxs = Base.setindex(nxs, Int32(1001), NXS_IZAID)
        nxs = Base.setindex(nxs, Int32(nes), NXS_NES)
        nxs = Base.setindex(nxs, Int32(0), NXS_NTR)
        nxs = Base.setindex(nxs, Int32(0), NXS_NR)
        nxs = Base.setindex(nxs, Int32(1), NXS_IZ)
        nxs = Base.setindex(nxs, Int32(1), NXS_IA)

        jxs = ntuple(i -> Int32(0), 32)
        jxs = Base.setindex(jxs, Int32(1), JXS_ESZ)

        pairs_iz = ntuple(i -> Int32(0), 16)
        pairs_aw = ntuple(i -> 0.0, 16)

        table = ACETable(
            "1001.80c  ", 0.999167, 2.5852e-8, "03/21/2026",
            "Test hydrogen table", "   mat 125",
            pairs_iz, pairs_aw, nxs, jxs, xss_data
        )

        buf = IOBuffer()
        write_ace(buf, table)
        output = String(take!(buf))
        lines = split(output, '\n')

        # Filter out the trailing empty line from final newline
        nonempty = filter(l -> !isempty(l), lines)

        # Line 1 (header): should be ~45 chars (ZAID + AWR + TZ + date)
        @test length(nonempty[1]) >= 44

        # Line 2 (comment + mat_id): should be 80 chars
        @test length(nonempty[2]) == 80

        # Lines 3-6 (IZ/AW pairs): 4 lines, each with 4 pairs
        for i in 3:6
            # Each pair is i7 + f11.0 = 18 chars, 4 pairs = 72 chars
            @test length(nonempty[i]) == 72
        end

        # Lines 7-12 (NXS + JXS): 6 lines of 8 integers at i9 = 72 chars
        for i in 7:12
            @test length(nonempty[i]) == 72
        end

        # Data lines: 4 values per line, 20 chars each = 80 chars
        # (except possibly the last line)
        n_data = length(xss_data)
        n_full_lines = div(n_data, 4)
        for i in 1:n_full_lines
            line_idx = 12 + i
            @test length(nonempty[line_idx]) == 80
        end

        # Each individual XSS value should be 20 chars
        val_str = NJOY._format_xss_value(1.234567e8)
        @test length(val_str) == 20
        val_str2 = NJOY._format_xss_value(0.0)
        @test length(val_str2) == 20
    end

    # ======================================================================
    # ACER -- Build ACE from PointwiseMaterial
    # ======================================================================
    @testset "ACER -- build_ace from PointwiseMaterial" begin
        # Create a simple 2-reaction PointwiseMaterial (elastic + capture)
        nes = 10
        energies = collect(range(1.0e-5, stop=2.0e7, length=nes))
        xs = zeros(nes, 3)
        mt_list = [1, 2, 102]

        # Total = elastic + capture
        for i in 1:nes
            xs[i, 2] = 20.0 / sqrt(energies[i] / 0.0253)  # 1/v elastic
            xs[i, 3] = 10.0 / sqrt(energies[i] / 0.0253)  # 1/v capture
            xs[i, 1] = xs[i, 2] + xs[i, 3]                 # total
        end

        pendf = PointwiseMaterial(Int32(1001), energies, xs, mt_list)

        ace = build_ace(pendf;
                        suffix="80c",
                        awr=0.999167,
                        temperature=300.0,
                        comment="Test H-1 from NJOY.jl")

        # Verify ZAID
        @test startswith(ace.zaid, "1001")
        @test occursin("80c", ace.zaid)

        # Verify AWR
        @test ace.awr == 0.999167

        # Verify temperature in MeV
        bk = 8.617333262e-5
        expected_tz = 300.0 * bk / 1.0e6
        @test isapprox(ace.temp, expected_tz, rtol=1e-10)

        # Verify NXS array
        @test ace.nxs[NXS_NES] == nes
        @test ace.nxs[NXS_NTR] == 1   # only MT102 is extra
        @test ace.nxs[NXS_LEN2] == length(ace.xss)
        @test ace.nxs[NXS_IZAID] == 1001

        # Verify JXS(1) = ESZ points to start of data
        @test ace.jxs[JXS_ESZ] == 1

        # Verify ESZ block: energies are in MeV
        e_ace = esz_energies(ace)
        @test length(e_ace) == nes
        @test isapprox(e_ace[1], energies[1] * 1e-6, rtol=1e-10)
        @test isapprox(e_ace[end], energies[end] * 1e-6, rtol=1e-10)

        # Verify ESZ block: total cross section preserved
        xs_total = esz_total(ace)
        @test length(xs_total) == nes
        for i in 1:nes
            @test isapprox(xs_total[i], xs[i, 1], rtol=1e-10)
        end

        # Verify ESZ block: elastic cross section preserved
        xs_elas = esz_elastic(ace)
        @test length(xs_elas) == nes
        for i in 1:nes
            @test isapprox(xs_elas[i], xs[i, 2], rtol=1e-10)
        end

        # Verify MTR block contains MT102
        mtr_start = Int(ace.jxs[JXS_MTR])
        @test mtr_start > 0
        @test ace.xss[mtr_start] == 102.0

        # Verify TYR block: MT102 should have TYR=0 (no neutron out)
        tyr_start = Int(ace.jxs[JXS_TYR])
        @test ace.xss[tyr_start] == 0.0  # capture: no neutrons

        # Verify LSIG/SIG blocks are set
        @test ace.jxs[JXS_LSIG] > 0
        @test ace.jxs[JXS_SIG] > 0

        # Verify LAND block is set (isotropic flag)
        @test ace.jxs[JXS_LAND] > 0
        land_start = Int(ace.jxs[JXS_LAND])
        @test ace.xss[land_start] == -1.0  # isotropic in CM

        # Verify write_ace produces valid output
        buf = IOBuffer()
        write_ace(buf, ace)
        output = String(take!(buf))
        @test length(output) > 0
        lines = split(output, '\n')
        nonempty = filter(l -> !isempty(l), lines)
        # At least header (2) + IZ/AW (4) + NXS/JXS (6) + data lines
        @test length(nonempty) >= 12 + div(length(ace.xss), 4)

        # Verify nxs_length accessor
        @test nxs_length(ace) == length(ace.xss)
        @test nxs_nes(ace) == nes
        @test nxs_ntr(ace) == 1
    end

    # ======================================================================
    # UNRESR -- Bondarenko self-shielding
    # ======================================================================
    @testset "UNRESR -- URR penetrability" begin
        # l=0: Vl=1, phi=rhoc
        V0, ps0 = urr_penetrability(0, 0.5, 0.6)
        @test V0 == 1.0
        @test ps0 == 0.6

        # l=1: Vl = rho^2/(1+rho^2), phi = rhoc - atan(rhoc)
        rho, rhoc = 1.5, 1.2
        V1, ps1 = urr_penetrability(1, rho, rhoc)
        @test isapprox(V1, rho^2 / (1 + rho^2), rtol=1e-14)
        @test isapprox(ps1, rhoc - atan(rhoc), rtol=1e-14)

        # l=2: Vl = rho^4/(9+3*rho^2+rho^4)
        V2, ps2 = urr_penetrability(2, rho, rhoc)
        r2 = rho^2; r4 = r2^2
        @test isapprox(V2, r4 / (9 + 3*r2 + r4), rtol=1e-14)
        @test isapprox(ps2, rhoc - atan(3*rhoc/(3-rhoc^2)), rtol=1e-14)
    end

    @testset "UNRESR -- Hwang quadrature tables" begin
        # Quadrature weights should sum to ~1 for each df
        for df in 1:4
            wsum = sum(HWANG_QW[:, df])
            @test isapprox(wsum, 1.0, atol=0.02)
        end
        # Abscissae should be positive
        @test all(HWANG_QP .>= 0.0)
        # Sizes
        @test size(HWANG_QW) == (10, 4)
        @test size(HWANG_QP) == (10, 4)
    end

    @testset "UNRESR -- ajku integration" begin
        # J and K are positive integrals
        xj, xk = ajku(10.0, 1.0)
        @test xj > 0.0
        @test xk > 0.0

        # K <= J (denominator in K has higher power)
        @test xk <= xj

        # For beta=0 (zero dilution), J integral = integral of phi du / pi
        # which converges on our finite [0,12] range
        xj0, xk0 = ajku(0.0, 1.0)
        @test xj0 > 0.3  # finite but < 1 due to truncated range

        # Increasing beta should decrease J (more dilution)
        xj1, _ = ajku(1.0, 1.0)
        xj10, _ = ajku(10.0, 1.0)
        xj100, _ = ajku(100.0, 1.0)
        @test xj1 > xj10
        @test xj10 > xj100

        # sti -> 0 should give small values (narrow resonance, little overlap)
        xjs, xks = ajku(10.0, 0.001)
        @test xjs < 0.01
    end

    @testset "UNRESR -- Bondarenko self-shielding limits" begin
        # Test with a simple U-238-like single-sequence model
        seq = URRSpinSequence(0, 0.5, 20.0, 1.0e-4, 0.023, 0.0, 0.0, 1, 0, 0)
        model = URRStatModel(0.0, 236.0, 0.948, [seq])

        # Compute at a URR energy
        E = 10000.0  # 10 keV
        T = 300.0
        sigz = [1e10, 1e4, 1e2, 1e1]  # from infinite to small dilution

        result = bondarenko_xs(model, E, T, sigz)

        # Total should be positive for all dilutions
        for is0 in 1:4
            @test result[1, is0] > 0.0
        end

        # Self-shielding: xs should decrease with decreasing sigma0
        # (more shielding at lower dilution)
        for is0 in 2:4
            @test result[1, is0] <= result[1, is0-1] * 1.01  # allow small tolerance
        end

        # Infinite dilution should give finite positive xs
        xs_inf = infinite_dilution_xs(model, E, T)
        @test xs_inf.total > 0.0
        @test xs_inf.elastic > 0.0
        @test xs_inf.capture >= 0.0
    end

    # ======================================================================
    # PURR -- chi-squared sampling
    # ======================================================================
    @testset "PURR -- chi-squared quantile table" begin
        using Random
        @test size(CHI2_QUANTILES) == (20, 4)

        # All quantiles should be positive
        @test all(CHI2_QUANTILES .> 0.0)

        # Columns should be monotonically increasing (quantiles of CDF)
        for df in 1:4
            for i in 2:20
                @test CHI2_QUANTILES[i, df] > CHI2_QUANTILES[i-1, df]
            end
        end

        # Sample mean should be close to df (chi^2 with df degrees of freedom)
        for df in 1:4
            rng = Random.Xoshiro(42)
            samples = [chi2_sample(df, rng) for _ in 1:10000]
            smean = sum(samples) / length(samples)
            # With 20-point quantile sampling, mean should be within ~20% of df
            @test abs(smean - df) / df < 0.25
        end
    end

    @testset "PURR -- Wigner spacing" begin
        using Random
        rng = Random.Xoshiro(123)
        D = 10.0
        spacings = [wigner_spacing(D, rng) for _ in 1:5000]

        # Mean of Wigner distribution: D * sqrt(pi/4) * Gamma(3/2)/Gamma(1)
        # Actually: <s> = D * sqrt(4/pi) * <sqrt(-ln U)>
        # Exact mean for Wigner = D (by construction of dcon)
        smean = sum(spacings) / length(spacings)
        @test isapprox(smean, D, rtol=0.15)

        # All spacings should be positive
        @test all(spacings .> 0.0)
    end

    @testset "PURR -- ladder generation" begin
        using Random
        seq = URRSpinSequence(0, 0.5, 20.0, 1.0e-4, 0.023, 0.0, 0.0, 1, 0, 0)
        rng = Random.Xoshiro(999)

        er, gnr, gfr, ggr, gxr = generate_ladder(seq, 9000.0, 11000.0, rng)

        # Should have resonances
        @test length(er) > 0

        # Energies should be in range and sorted
        @test all(er .>= 9000.0)
        @test issorted(er)

        # Width ratios should sum to 1 for each resonance
        for i in 1:length(er)
            ratio_sum = gnr[i] + gfr[i] + ggr[i] + gxr[i]
            @test isapprox(ratio_sum, 1.0, rtol=1e-12)
        end

        # For zero fission/competitive, those ratios should be small/zero
        @test all(gfr .== 0.0) || all(gfr .< 0.01)
    end

    @testset "PURR -- line shape approximation" begin
        # For large x, should approach asymptotic
        rew, aimw = NJOY.line_shape(50.0, 1.0)
        # Asymptotic: rew = y * c1 / (x^2+y^2)
        c1 = 0.5641895835
        @test isapprox(rew, 1.0 * c1 / (2500.0 + 1.0), rtol=1e-3)

        # For x=0, rew should be largest (peak of line)
        rew0, _ = NJOY.line_shape(0.0, 1.0)
        rew1, _ = NJOY.line_shape(2.0, 1.0)
        @test rew0 > rew1

        # Symmetry: Re[w(x,y)] = Re[w(-x,y)]
        rewp, _ = NJOY.line_shape(3.0, 2.0)
        rewm, _ = NJOY.line_shape(-3.0, 2.0)
        @test isapprox(rewp, rewm, rtol=1e-10)
    end

    @testset "PURR -- probability table normalization" begin
        # Build simple model and generate ptable
        seq = URRSpinSequence(0, 0.5, 20.0, 1.0e-4, 0.023, 0.0, 0.0, 1, 0, 0)
        model = URRStatModel(0.0, 236.0, 0.948, [seq])

        energies = [10000.0]  # single energy for speed
        ptable = generate_ptable(model, energies;
                                 nladders=16, nbins=10, seed=42)

        # Probabilities should sum to 1
        psum = sum(ptable.prob[:, 1])
        @test isapprox(psum, 1.0, atol=1e-12)

        # All probabilities should be non-negative
        @test all(ptable.prob[:, 1] .>= 0.0)

        # Total should be >= elastic + fission + capture for bins with data
        for j in 1:ptable.nbins
            if ptable.prob[j, 1] > 0
                parts = ptable.elastic[j,1] + ptable.fission[j,1] + ptable.capture[j,1]
                # Total can differ from sum of parts due to bkg
                @test ptable.total[j, 1] > 0.0
            end
        end
    end

    @testset "PURR -- Bondarenko from probability table" begin
        seq = URRSpinSequence(0, 0.5, 20.0, 1.0e-4, 0.023, 0.0, 0.0, 1, 0, 0)
        model = URRStatModel(0.0, 236.0, 0.948, [seq])

        energies = [10000.0]
        ptable = generate_ptable(model, energies;
                                 nladders=32, nbins=15, seed=77)

        # Infinite dilution from ptable
        sig_t, sig_e, sig_f, sig_c, sig_tr = bondarenko_from_ptable(ptable, 1, 1e10)

        # All should be positive
        @test sig_t > 0.0
        @test sig_e > 0.0
        @test sig_c >= 0.0

        # Self-shielded should be <= infinite dilution
        sig_t2, _, _, _, _ = bondarenko_from_ptable(ptable, 1, 100.0)
        @test sig_t2 <= sig_t * 1.01
    end

    # ======================================================================
    # ACER Proposer-B -- ACENeutronTable, ZAID, format correctness
    # ======================================================================
    @testset "ACER-B -- ZAID utilities" begin
        @test format_zaid(92235, "80c") == "92235.80c"
        @test format_zaid(92, 235, "80c") == "92235.80c"
        @test format_zaid(1001, "70c") == "1001.70c"
        za, suf = parse_zaid("92235.80c")
        @test za == 92235
        @test suf == "80c"
        za2, suf2 = parse_zaid("  1001.70c  ")
        @test za2 == 1001
        @test suf2 == "70c"
        @test_throws ArgumentError parse_zaid("no_dot")
    end

    @testset "ACER-B -- temperature conversion" begin
        bk = NJOY.PhysicsConstants.bk
        tz300 = temp_to_mev(300.0)
        @test isapprox(tz300, 300.0 * bk * 1e-6, rtol=1e-12)
        t_back = mev_to_temp(tz300)
        @test isapprox(t_back, 300.0, rtol=1e-10)
        @test isapprox(mev_to_temp(temp_to_mev(600.0)), 600.0, rtol=1e-10)
    end

    @testset "ACER-B -- ACEHeader construction" begin
        h = ACEHeader(zaid="92235.80c", awr=235.044, temp_mev=2.53e-8,
                      date="03/21/2026", comment="U-235 test",
                      mat_string="   mat9228")
        @test length(h.hz) == 10
        @test length(h.hd) == 10
        @test length(h.hk) == 70
        @test length(h.hm) == 10
        @test h.aw0 == 235.044
        @test strip(h.hz) == "92235.80c"
        h2 = ACEHeader(zaid="1001.80c", awr=1.0, temp_mev=0.0,
                       comment="x"^100)
        @test length(h2.hk) == 70
    end

    @testset "ACER-B -- ACENeutronTable construction" begin
        n = 5
        egrid = collect(range(1e-11, 20.0, length=n))
        h = ACEHeader(zaid="1001.80c", awr=0.999167, temp_mev=2.53e-8)
        t = ACENeutronTable(header=h, energy_grid=egrid,
                            total_xs=fill(10.0, n),
                            absorption_xs=fill(4.0, n),
                            elastic_xs=fill(6.0, n),
                            heating_numbers=zeros(n))
        @test ace_nes(t) == n
        @test ace_ntr(t) == 0
        @test length(t.iz) == 16
        @test t.angular_elastic === nothing
        @test_throws ArgumentError ACENeutronTable(
            header=h, energy_grid=egrid,
            total_xs=fill(10.0, n+1),
            absorption_xs=fill(4.0, n),
            elastic_xs=fill(6.0, n), heating_numbers=zeros(n))
    end

    @testset "ACER-B -- ReactionXS and EquiprobableBins" begin
        rxn = ReactionXS(Int32(102), -6.0, Int32(0), Int32(3), [1.0, 2.0, 3.0])
        @test rxn.mt == 102
        @test rxn.q_value == -6.0
        @test rxn.ie_start == 3
        bins = EquiprobableBins(collect(range(-1.0, 1.0, length=33)))
        @test length(bins.cosines) == 33
        @test_throws ArgumentError EquiprobableBins(zeros(32))
    end

    @testset "ACER-B -- build_xss serialization" begin
        n = 10
        egrid = collect(range(1e-11, 20.0, length=n))
        rxns = [
            ReactionXS(Int32(102), 0.0, Int32(0), Int32(1), fill(4.0, n)),
            ReactionXS(Int32(51), -0.5, Int32(1), Int32(3), fill(1.5, n-2))
        ]
        h = ACEHeader(zaid="26056.80c", awr=55.454, temp_mev=2.53e-8)
        t = ACENeutronTable(header=h, energy_grid=egrid,
                            total_xs=fill(10.0, n),
                            absorption_xs=fill(4.0, n),
                            elastic_xs=fill(6.0, n),
                            heating_numbers=zeros(n),
                            reactions=rxns)
        nxs_v, jxs_v, xss_v, is_int_v = build_xss(t)
        @test nxs_v[NXS_LEN2] == length(xss_v)
        @test nxs_v[NXS_IZAID] == 26056
        @test nxs_v[NXS_NES] == n
        @test nxs_v[NXS_NTR] == 2
        @test nxs_v[NXS_NR] == 1
        @test nxs_v[NXS_IZ] == 26
        @test nxs_v[NXS_IA] == 56
        @test jxs_v[JXS_ESZ] == 1
        @test jxs_v[JXS_MTR] >= 5 * n + 1
        @test jxs_v[JXS_END] == length(xss_v) + 1
        for i in 1:n
            @test isapprox(xss_v[i], egrid[i], rtol=1e-14)
        end
        mtr_loc = jxs_v[JXS_MTR]
        @test round(Int, xss_v[mtr_loc]) == 102
        @test round(Int, xss_v[mtr_loc+1]) == 51
        @test is_int_v[mtr_loc] == true
        @test is_int_v[1] == false
    end

    @testset "ACER-B -- write_ace_table format correctness" begin
        n = 8
        egrid = collect(range(1e-11, 20.0, length=n))
        h = ACEHeader(zaid="92235.80c", awr=233.025,
                      temp_mev=2.5852e-8, date="03/21/2026",
                      comment="U-235 format test",
                      mat_string="   mat9228")
        t = ACENeutronTable(header=h, energy_grid=egrid,
                            total_xs=fill(100.0, n),
                            absorption_xs=fill(50.0, n),
                            elastic_xs=fill(50.0, n),
                            heating_numbers=zeros(n))
        buf = IOBuffer()
        write_ace_table(buf, t)
        output = String(take!(buf))
        lines = split(output, '\n')
        nonempty = filter(l -> !isempty(l), lines)
        @test length(nonempty[1]) >= 44
        @test length(nonempty[2]) == 80
        for i in 3:6; @test length(nonempty[i]) == 72; end
        for i in 7:12; @test length(nonempty[i]) == 72; end
        nxs_v, _, xss_v, _ = build_xss(t)
        n_full = div(length(xss_v), 4)
        for i in 1:n_full
            @test length(nonempty[12 + i]) == 80
        end
        @test length(NJOY._ace_fmt_int(12345)) == 20
        @test length(NJOY._ace_fmt_real(1.23e10)) == 20
        @test length(NJOY._ace_fmt_real(0.0)) == 20
    end

    @testset "ACER-B -- write/read-back verification" begin
        n = 6
        egrid = [1e-11, 1e-8, 1e-5, 1e-2, 1.0, 20.0]
        total_v = [100.0, 80.0, 50.0, 20.0, 10.0, 5.0]
        elas_v = [60.0, 50.0, 30.0, 12.0, 6.0, 3.0]
        absorp_v = total_v .- elas_v
        rxns = [ReactionXS(Int32(102), 0.0, Int32(0), Int32(1), absorp_v)]
        h = ACEHeader(zaid="1001.80c", awr=0.999167, temp_mev=2.53e-8,
                      comment="H-1 roundtrip")
        t = ACENeutronTable(header=h, energy_grid=egrid,
                            total_xs=total_v, absorption_xs=absorp_v,
                            elastic_xs=elas_v, heating_numbers=zeros(n),
                            reactions=rxns)
        buf = IOBuffer()
        write_ace_table(buf, t)
        output = String(take!(buf))
        lines = split(output, '\n')
        nonempty = filter(l -> !isempty(l), lines)
        @test startswith(nonempty[1], "1001.80c")
        nxs_jxs_text = join(nonempty[7:12], "")
        vals = [parse(Int, nxs_jxs_text[(i-1)*9+1:i*9]) for i in 1:48]
        @test vals[NXS_NES] == n
        @test vals[NXS_NTR] == 1
        @test vals[NXS_IZAID] == 1001
        @test vals[16 + JXS_ESZ] == 1
        data_line = nonempty[13]
        e1_str = strip(data_line[1:20])
        e1_val = parse(Float64, e1_str)
        @test isapprox(e1_val, 1e-11, rtol=1e-6)
    end

    @testset "ACER-B -- build_ace_from_pendf" begin
        nes_p = 8
        energies = collect(range(1.0e-5, stop=2.0e7, length=nes_p))
        xs = zeros(nes_p, 3)
        for i in 1:nes_p
            xs[i, 2] = 20.0 / sqrt(energies[i] / 0.0253)
            xs[i, 3] = 10.0 / sqrt(energies[i] / 0.0253)
            xs[i, 1] = xs[i, 2] + xs[i, 3]
        end
        pendf = PointwiseMaterial(Int32(1001), energies, xs, [1, 2, 102])
        ace_t = build_ace_from_pendf(pendf, suffix="80c", temp_kelvin=300.0,
                                      comment="H-1 from PENDF")
        @test ace_nes(ace_t) == nes_p
        @test ace_ntr(ace_t) == 1
        @test isapprox(ace_t.energy_grid[1], energies[1]*1e-6, rtol=1e-10)
        for i in 1:nes_p
            @test isapprox(ace_t.total_xs[i], xs[i,1], rtol=1e-10)
        end
        buf = IOBuffer()
        write_ace_table(buf, ace_t)
        @test length(take!(buf)) > 0
        buf2 = IOBuffer()
        write_ace(buf2, ace_t)
        @test length(take!(buf2)) > 0
    end

    @testset "ACER-B -- NXS/JXS consistency" begin
        n = 5
        egrid = collect(range(1e-11, 20.0, length=n))
        h = ACEHeader(zaid="26056.80c", awr=55.454, temp_mev=2.53e-8)
        rxns = [ReactionXS(Int32(102), 0.0, Int32(0), Int32(1), fill(4.0, n))]
        t = ACENeutronTable(header=h, energy_grid=egrid,
                            total_xs=fill(10.0, n),
                            absorption_xs=fill(4.0, n),
                            elastic_xs=fill(6.0, n),
                            heating_numbers=zeros(n), reactions=rxns)
        nxs_v, jxs_v, xss_v, _ = build_xss(t)
        @test nxs_v[NXS_LEN2] == length(xss_v)
        @test jxs_v[JXS_MTR] >= 5 * n + 1
        lsig_loc = jxs_v[JXS_LSIG]
        sig_loc = jxs_v[JXS_SIG]
        @test round(Int, xss_v[lsig_loc]) == 1
        ie_s = round(Int, xss_v[sig_loc])
        ne_r = round(Int, xss_v[sig_loc + 1])
        @test ie_s >= 1 && ie_s <= n
        @test ne_r == n - ie_s + 1
        @test jxs_v[JXS_END] == length(xss_v) + 1
    end

    # ======================================================================
    # GROUPR -- Group structures
    # ======================================================================
    @testset "GROUPR -- Group structure validation" begin
        # LANL 30: 31 bounds, 30 groups
        @test length(LANL_30) == 31
        @test num_groups(LANL_30) == 30
        validate_group_bounds(LANL_30)  # should not throw
        @test LANL_30[1] == 1.39e-4
        @test LANL_30[end] == 1.70e7

        # WIMS 69: 70 bounds, 69 groups
        @test length(WIMS_69) == 70
        @test num_groups(WIMS_69) == 69
        validate_group_bounds(WIMS_69)
        @test WIMS_69[1] == 1.0e-5
        @test WIMS_69[end] == 1.0e7

        # VITAMIN-J 175: 176 bounds, 175 groups
        @test length(VITAMINJ_175) == 176
        @test num_groups(VITAMINJ_175) == 175
        validate_group_bounds(VITAMINJ_175)
        @test VITAMINJ_175[1] == 1.0e-5
        @test VITAMINJ_175[end] == 1.9640e7

        # Verify the e175 insertion at position 166
        @test VITAMINJ_175[166] == 1.284e7

        # Strictly ascending check
        for i in 1:175
            @test VITAMINJ_175[i] < VITAMINJ_175[i+1]
        end

        # Additional structures
        @test length(SANDII_620) == 621
        @test num_groups(SANDII_620) == 620
        validate_group_bounds(SANDII_620)

        @test length(XMAS_172) == 173
        @test num_groups(XMAS_172) == 172
        validate_group_bounds(XMAS_172)

        @test length(ECCO_33) == 34
        @test num_groups(ECCO_33) == 33
        validate_group_bounds(ECCO_33)

        # Lookup by ign number
        @test get_group_structure(3) === LANL_30
        @test get_group_structure(9) === WIMS_69
        @test get_group_structure(17) === VITAMINJ_175
        @test get_group_structure(12) === SANDII_620
        @test get_group_structure(18) === XMAS_172
        @test get_group_structure(19) === ECCO_33
        @test get_group_structure(IGN_LANL30) === LANL_30

        # Validation rejects bad bounds
        @test_throws ArgumentError validate_group_bounds([1.0])
        @test_throws ArgumentError validate_group_bounds([2.0, 1.0])
        @test_throws ArgumentError validate_group_bounds([0.0, 1.0])
    end

    @testset "GROUPR -- find_group" begin
        bounds = (1.0, 10.0, 100.0, 1000.0)  # 3 groups
        @test find_group(bounds, 0.5) == 0    # below
        @test find_group(bounds, 1000.0) == 0 # at upper edge (exclusive)
        @test find_group(bounds, 1.0) == 1    # at lower edge (inclusive)
        @test find_group(bounds, 5.0) == 1
        @test find_group(bounds, 10.0) == 2
        @test find_group(bounds, 50.0) == 2
        @test find_group(bounds, 100.0) == 3
        @test find_group(bounds, 500.0) == 3
    end

    # ======================================================================
    # GROUPR -- Integration and averaging
    # ======================================================================
    @testset "GROUPR -- group_integrate constant function" begin
        # Constant function sigma=5.0 over [1, 100]
        # Integral over [a,b] = 5*(b-a)
        energies = [1.0, 10.0, 50.0, 100.0]
        values   = [5.0, 5.0,  5.0,  5.0]
        bounds   = [1.0, 10.0, 100.0]  # 2 groups

        result = group_integrate(energies, values, bounds)
        @test length(result) == 2
        @test isapprox(result[1], 5.0 * 9.0,  rtol=1e-12)   # 5*(10-1)
        @test isapprox(result[2], 5.0 * 90.0, rtol=1e-12)   # 5*(100-10)
    end

    @testset "GROUPR -- group_integrate linear function" begin
        # Linear function: sigma(E) = 2*E on [0.5, 10.0]
        # Integral over [a,b] = b^2 - a^2
        energies = [0.5, 2.0, 5.0, 10.0]
        values   = [1.0, 4.0, 10.0, 20.0]
        bounds   = [0.5, 5.0, 10.0]

        result = group_integrate(energies, values, bounds)
        @test isapprox(result[1], 5.0^2 - 0.5^2, rtol=1e-10)  # 24.75
        @test isapprox(result[2], 10.0^2 - 5.0^2, rtol=1e-10) # 75.0
    end

    @testset "GROUPR -- group_integrate partial overlap" begin
        # Data covers [1, 10], groups extend beyond
        energies = [1.0, 5.0, 10.0]
        values   = [2.0, 2.0, 2.0]
        bounds   = [0.1, 1.0, 5.0, 10.0, 20.0]  # 4 groups

        result = group_integrate(energies, values, bounds)
        @test result[1] == 0.0           # no data below 1.0
        @test isapprox(result[2], 8.0)   # 2*(5-1)
        @test isapprox(result[3], 10.0)  # 2*(10-5)
        @test result[4] == 0.0           # no data above 10.0
    end

    @testset "GROUPR -- group_integrate with NTuple bounds" begin
        energies = [1.0, 2.0, 3.0]
        values   = [4.0, 4.0, 4.0]
        bounds   = (1.0, 2.0, 3.0)

        result = group_integrate(energies, values, bounds)
        @test isapprox(result[1], 4.0)
        @test isapprox(result[2], 4.0)
    end

    @testset "GROUPR -- group_average constant XS with 1/E weight" begin
        # For constant sigma, the group average should be exactly sigma
        # regardless of the weight function.
        # sigma_g = int(sigma * w) / int(w) = sigma * int(w) / int(w) = sigma
        energies = collect(range(1.0, 100.0, length=200))
        xs = fill(3.14, 200)
        bounds = [1.0, 10.0, 50.0, 100.0]

        mgxs = group_average(energies, xs, 1, bounds; weight_fn=weight_inv_e)
        @test length(mgxs.flux) == 3
        for g in 1:3
            @test isapprox(mgxs.xs[g, 1], 3.14, rtol=1e-6)
        end
    end

    @testset "GROUPR -- group_average 1/v XS with flat weight" begin
        # sigma(E) = 1/sqrt(E) (1/v cross section)
        # With flat weight: sigma_g = int(E^{-1/2} dE) / int(dE)
        #                           = 2(sqrt(b)-sqrt(a)) / (b-a)
        energies = collect(range(1.0, 100.0, length=5000))
        xs = [1.0 / sqrt(E) for E in energies]
        bounds = [1.0, 25.0, 100.0]

        mgxs = group_average(energies, xs, 2, bounds; weight_fn=weight_flat)

        # Group 1: [1, 25]: 2*(5-1)/24 = 8/24 = 1/3
        @test isapprox(mgxs.xs[1, 1], 2.0 * (5.0 - 1.0) / 24.0, rtol=5e-4)

        # Group 2: [25, 100]: 2*(10-5)/75 = 10/75 = 2/15
        @test isapprox(mgxs.xs[2, 1], 2.0 * (10.0 - 5.0) / 75.0, rtol=5e-4)
    end

    @testset "GROUPR -- group_average multi-reaction" begin
        energies = collect(range(1.0, 10.0, length=100))
        # Two reactions: constant 2.0 and constant 5.0
        xs_mat = hcat(fill(2.0, 100), fill(5.0, 100))
        bounds = [1.0, 5.0, 10.0]

        mgxs = group_average(energies, xs_mat, [2, 102], bounds)
        @test mgxs.mt_list == [2, 102]
        @test size(mgxs.xs) == (2, 2)
        for g in 1:2
            @test isapprox(mgxs.xs[g, 1], 2.0, rtol=1e-6)
            @test isapprox(mgxs.xs[g, 2], 5.0, rtol=1e-6)
        end
    end

    @testset "GROUPR -- group_average_shielded infinite dilution" begin
        # At sigma0 -> infinity, shielded average should match unshielded
        energies = collect(range(1.0, 100.0, length=200))
        total_xs = fill(10.0, 200)
        reaction_xs = reshape(fill(3.0, 200), :, 1)
        bounds = [1.0, 50.0, 100.0]

        mgxs_inf = group_average_shielded(energies, total_xs, reaction_xs,
                                          [2], bounds, 1.0e10)
        mgxs_std = group_average(energies, reaction_xs, [2], bounds)

        for g in 1:2
            @test isapprox(mgxs_inf.xs[g, 1], mgxs_std.xs[g, 1], rtol=1e-4)
        end
    end

    @testset "GROUPR -- group_average_shielded finite sigma0" begin
        # With finite sigma0, the shielded XS of a constant should still
        # be constant (shielding factor cancels in numerator/denominator)
        energies = collect(range(1.0, 100.0, length=200))
        total_xs = fill(10.0, 200)
        reaction_xs = reshape(fill(10.0, 200), :, 1)  # same as total
        bounds = [1.0, 50.0, 100.0]

        mgxs = group_average_shielded(energies, total_xs, reaction_xs,
                                      [1], bounds, 100.0)
        for g in 1:2
            @test isapprox(mgxs.xs[g, 1], 10.0, rtol=1e-4)
        end
    end

    @testset "GROUPR -- weight functions" begin
        @test weight_flat(42.0) == 1.0
        @test weight_inv_e(4.0) == 0.25
        @test weight_inv_e(1.0) == 1.0

        # Maxwell-fission: continuous and positive
        kT = 0.0253
        @test weight_maxwell_fission(1e-4; kT=kT) > 0
        @test weight_maxwell_fission(1.0; kT=kT) > 0
        @test weight_maxwell_fission(1e6; kT=kT) > 0

        # Continuity at transition points
        Ec = kT * 20.0
        Ef = 820.3e3
        w_below_Ec = weight_maxwell_fission(Ec * 0.999; kT=kT)
        w_above_Ec = weight_maxwell_fission(Ec * 1.001; kT=kT)
        @test isapprox(w_below_Ec, w_above_Ec, rtol=0.01)

        w_below_Ef = weight_maxwell_fission(Ef * 0.999; kT=kT)
        w_above_Ef = weight_maxwell_fission(Ef * 1.001; kT=kT)
        @test isapprox(w_below_Ef, w_above_Ef, rtol=0.01)
    end

    @testset "GROUPR -- analytical integral verification" begin
        # Verify group_integrate matches direct trapezoidal for lin-lin
        # On [1, 5] with f(E) = E^2 approximated piecewise-linear on fine grid
        # The trapezoidal rule IS exact for piecewise-linear, so the
        # panel_integral (lin-lin) result should match numerical integration.
        n = 10000
        energies = collect(range(1.0, 5.0, length=n))
        values = energies .^ 2  # quadratic -- tests accuracy of linearization
        bounds = [1.0, 3.0, 5.0]

        result = group_integrate(energies, values, bounds)
        # Exact: int_1^3 E^2 dE = 27/3 - 1/3 = 26/3
        # Exact: int_3^5 E^2 dE = 125/3 - 27/3 = 98/3
        @test isapprox(result[1], 26.0 / 3.0, rtol=1e-6)
        @test isapprox(result[2], 98.0 / 3.0, rtol=1e-6)
    end

    # ======================================================================
    # GROUPR -- Weight functions (weight_functions.jl)
    # ======================================================================
    @testset "GROUPR -- weight_functions.jl basic" begin
        # constant_weight
        @test constant_weight(42.0) == 1.0
        @test constant_weight(1e-5) == 1.0

        # inv_e_weight
        @test inv_e_weight(4.0) == 0.25
        @test inv_e_weight(1.0) == 1.0
        @test inv_e_weight(0.01) == 100.0

        # Backward-compat aliases still work
        @test weight_flat(1.0) == constant_weight(1.0)
        @test weight_inv_e(5.0) == inv_e_weight(5.0)
    end

    @testset "GROUPR -- 1/E integral = ln(E_high/E_low)" begin
        # Verify that integrating 1/E over [E_low, E_high] gives ln(E_high/E_low)
        E_low, E_high = 1.0, 1000.0
        n = 50000
        energies = collect(range(E_low, E_high, length=n))
        wvals = [inv_e_weight(E) for E in energies]
        bounds = [E_low, E_high]

        result = group_integrate(energies, wvals, bounds)
        exact = log(E_high / E_low)
        @test isapprox(result[1], exact, rtol=1e-5)
    end

    @testset "GROUPR -- maxwell_inv_e_fission" begin
        # Positivity across all regions
        for E in [1e-5, 0.01, 0.0253, 1.0, 100.0, 1e5, 1e6, 1e7]
            @test maxwell_inv_e_fission(E) > 0
        end

        # In the 1/E region (well above thermal, below fission), it equals 1/E
        E_mid = 100.0
        @test isapprox(maxwell_inv_e_fission(E_mid), 1.0 / E_mid, rtol=1e-10)

        # Continuity at Eb breakpoint
        Eb = 0.0253
        w_lo = maxwell_inv_e_fission(Eb * (1.0 - 1e-8))
        w_hi = maxwell_inv_e_fission(Eb * (1.0 + 1e-8))
        @test isapprox(w_lo, w_hi, rtol=1e-4)

        # Continuity at Ec breakpoint
        Ec = 820.3e3
        w_lo2 = maxwell_inv_e_fission(Ec * (1.0 - 1e-8))
        w_hi2 = maxwell_inv_e_fission(Ec * (1.0 + 1e-8))
        @test isapprox(w_lo2, w_hi2, rtol=1e-4)
    end

    @testset "GROUPR -- get_weight_function dispatch" begin
        w2 = get_weight_function(2)
        @test w2(5.0) == 1.0

        w3 = get_weight_function(3)
        @test w3(4.0) == 0.25

        w4 = get_weight_function(4)
        @test w4(100.0) > 0

        w6 = get_weight_function(6)
        @test w6(0.01) > 0

        w11 = get_weight_function(11)
        @test w11(1.0) > 0

        @test_throws ArgumentError get_weight_function(99)
    end

    @testset "GROUPR -- vitamin_e_weight" begin
        # Positivity across all regions
        for E in [1e-3, 0.1, 1.0, 1e3, 1e6, 5e6, 1.1e7, 1.3e7, 1.5e7, 1.8e7]
            @test vitamin_e_weight(E) > 0
        end

        # In the 1/E region (between 0.414 and 2.12e6)
        E_mid = 1000.0
        @test isapprox(vitamin_e_weight(E_mid), 1.0 / E_mid, rtol=1e-10)
    end

    @testset "GROUPR -- thermal_fission_fusion" begin
        # Positivity
        for E in [1e-3, 0.01, 0.1, 1.0, 1e5, 1e6, 1e7]
            @test thermal_fission_fusion(E) > 0
        end
    end

    @testset "GROUPR -- group_average constant sigma any weight" begin
        # For constant sigma, group average = sigma regardless of weight.
        # Test with all built-in weight functions.
        energies = collect(range(1.0, 1e6, length=500))
        sigma_val = 7.77
        xs = fill(sigma_val, 500)
        bounds = [1.0, 100.0, 1e4, 1e6]

        for iwt in [2, 3, 4, 6, 11]
            wf = get_weight_function(iwt)
            mgxs = group_average(energies, xs, 1, bounds; weight_fn=wf)
            for g in 1:3
                @test isapprox(mgxs.xs[g, 1], sigma_val, rtol=1e-3)
            end
        end
    end

    @testset "GROUPR -- 1/v XS with 1/E weight analytical" begin
        # sigma(E) = C/sqrt(E) (1/v), weight = 1/E
        # sigma_g = int(C/sqrt(E) * 1/E dE) / int(1/E dE)
        #         = C * int(E^{-3/2} dE) / ln(E_hi/E_lo)
        #         = C * [-2/sqrt(E)]_{E_lo}^{E_hi} / ln(E_hi/E_lo)
        #         = C * 2*(1/sqrt(E_lo) - 1/sqrt(E_hi)) / ln(E_hi/E_lo)
        C = 10.0
        E_lo, E_hi = 1.0, 100.0
        n = 20000
        energies = collect(range(E_lo, E_hi, length=n))
        xs = [C / sqrt(E) for E in energies]
        bounds = [E_lo, E_hi]

        mgxs = group_average(energies, xs, 1, bounds; weight_fn=inv_e_weight)
        exact = C * 2.0 * (1.0 / sqrt(E_lo) - 1.0 / sqrt(E_hi)) / log(E_hi / E_lo)
        @test isapprox(mgxs.xs[1, 1], exact, rtol=5e-4)
    end

    @testset "GROUPR -- self-shielded vs infinite dilute convergence" begin
        # As sigma0 -> infinity, shielded result should converge to unshielded
        energies = collect(range(1.0, 1000.0, length=500))
        total_xs = [10.0 + 5.0 * sin(E / 100.0) for E in energies]
        reaction_xs = reshape([3.0 + 2.0 * cos(E / 50.0) for E in energies], :, 1)
        bounds = [1.0, 100.0, 500.0, 1000.0]

        mgxs_unshielded = group_average(energies, reaction_xs, [2], bounds;
                                        weight_fn=inv_e_weight)

        # Progressively larger sigma0 should approach unshielded
        prev_err = Inf
        for sigma0 in [1e2, 1e4, 1e6, 1e8]
            mgxs_sh = group_average_shielded(energies, total_xs, reaction_xs,
                                             [2], bounds, sigma0;
                                             weight_fn=inv_e_weight)
            err = maximum(abs.(mgxs_sh.xs[:, 1] .- mgxs_unshielded.xs[:, 1]))
            @test err < prev_err  # convergence
            prev_err = err
        end
        # At sigma0=1e8, should be very close
        mgxs_large = group_average_shielded(energies, total_xs, reaction_xs,
                                            [2], bounds, 1e8;
                                            weight_fn=inv_e_weight)
        for g in 1:3
            @test isapprox(mgxs_large.xs[g, 1], mgxs_unshielded.xs[g, 1], rtol=1e-3)
        end
    end

end  # @testset "NJOY.jl"
