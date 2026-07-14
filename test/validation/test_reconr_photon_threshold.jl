using Test
using NJOY

# Synthetic oracle for the negative-Q branch in RECONR lunion/emerge.
# Ref: njoy-reference/src/reconr.f90:1911-1943,4726-4753,4791-4794.
# lunion replaces the first TAB1 energy with the upward-rounded physical
# threshold and cascades any following non-increasing breakpoints before
# emerge evaluates the scratch TAB1.
@testset "RECONR photon negative-Q threshold follows lunion scratch TAB1" begin
    interp = InterpolationTable(Int32[3], [LinLin])
    tab = TabulatedFunction(interp, [1.0, 1.05, 2.0], [0.0, 0.5, 2.0])
    sec = MF3Section(Int32(51), 0.0, -1.0, tab, Int32(13), 10.0,
                     Int32(0), Int32(0))

    threshold = NJOY.round_sigfig(1.1, 7, +1)
    next_energy = NJOY.round_sigfig(threshold, 7, +1)
    grid = [threshold, next_energy, 1.5, 2.0]
    result = only(NJOY._reconstruct_photon_sections([sec], grid))

    expected_mid = NJOY.round_sigfig(
        terp1(next_energy, 0.5, 2.0, 2.0, 1.5, LinLin), 7, 0)
    @test result.energies == grid
    @test result.xs[1] == 0.0
    @test result.xs[3] == expected_mid
end

# Ref: reconr.f90:1875-1877. These skips must be shared by lunion input and
# emerge output; otherwise a section can be omitted from the grid but emitted.
@testset "RECONR photon redundant-section skips precede lunion and emerge" begin
    interp = InterpolationTable(Int32[2], [LinLin])
    tab = TabulatedFunction(interp, [1.0, 2.0], [0.0, 1.0])
    mt460 = MF3Section(460, 0.0, 0.0, tab, 12)
    mt501 = MF3Section(501, 0.0, 0.0, tab, 13)
    mt522 = MF3Section(522, 0.0, 0.0, tab, 13)
    subshell = MF3Section(534, 0.0, 0.0, tab, 23)

    @test isempty(NJOY._eligible_photon_sections([mt460, mt501], MF3Section[]))
    @test isempty(NJOY._eligible_photon_sections([mt522], [subshell]))
    @test NJOY._eligible_photon_sections([mt522], MF3Section[]) == [mt522]
end
