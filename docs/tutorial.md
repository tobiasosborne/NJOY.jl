# Tutorial: Processing U-238 from ENDF to ACE

This tutorial walks through the complete NJOY.jl processing chain for U-238,
starting from an ENDF-6 evaluated data file and producing an ACE library
suitable for MCNP continuous-energy transport.

The pipeline mirrors the classic NJOY2016 sequence
RECONR -> BROADR -> HEATR -> ACER, but expressed as composable Julia
function calls instead of card-image input decks.

## Prerequisites

```julia
using NJOY
```

You will need an ENDF-6 file for U-238. The standard distribution is
`n-092_U_238.endf` from the ENDF/B-VIII.0 library.

## Step 1: Reconstruct pointwise cross sections (RECONR)

The `reconstruct` function reads the ENDF file, extracts MF2 resonance
parameters and MF3 background cross sections, and produces a linearized
pointwise representation (PENDF) using adaptive grid refinement.

```julia
endf_file = "n-092_U_238.endf"

# Reconstruct at 0.1% tolerance (default)
pendf = reconstruct(endf_file; err=0.001)
```

The return value is a `PointwiseMaterial` containing:
- `pendf.energies` -- energy grid in eV
- `pendf.cross_sections` -- matrix of shape (n_energies, n_reactions)
- `pendf.mt_list` -- reaction identifiers: `[1, 2, 18, 102]`
  (total, elastic, fission, capture)

You can inspect the result:

```julia
println("Grid points: ", length(pendf.energies))
println("Energy range: ", pendf.energies[1], " -- ", pendf.energies[end], " eV")
println("Reactions: ", pendf.mt_list)
```

### How it works internally

The RECONR pipeline is a chain of composable functions:

1. `read_mf2(io)` -- parse File 2 resonance parameters into typed structs
2. `build_evaluator(mf2)` -- return a closure `f(E) -> (total, elastic, fission, capture)`
3. `build_grid(mf2, mf3_sections)` -- union of resonance nodes and MF3 breakpoints
4. `adaptive_reconstruct(f, grid, config)` -- generic adaptive linearization
5. `merge_background!(energies, values, mf3, mf2)` -- add MF3 backgrounds

Each step is independently testable. The evaluator closure captures the MF2
data and dispatches to `cross_section_slbw`, `cross_section_mlbw`, or
`cross_section_rm` depending on the resonance formalism.

## Step 2: Doppler-broaden to reactor temperature (BROADR)

Real applications require cross sections at a finite temperature. The
`doppler_broaden` function applies the exact SIGMA1 kernel to broaden from
0 K to the target temperature.

```julia
# U-238 atomic weight ratio
awr = 236.0058  # AWR for U-238 (ratio to neutron mass)

# Broaden to 300 K (room temperature)
pendf_300K = doppler_broaden(pendf, 300.0; awr=awr, tol=0.001)

println("Broadened grid points: ", length(pendf_300K.energies))
```

The broadened result is a new `PointwiseMaterial` on a refined grid. All four
reaction channels (total, elastic, fission, capture) are broadened
simultaneously on a shared adaptive grid.

For a single cross section vector, you can also use the lower-level interface:

```julia
# Broaden just elastic (column 2)
new_e, new_xs = doppler_broaden(
    pendf.energies,
    pendf.cross_sections[:, 2],  # elastic
    300.0, awr;
    tol=0.001
)
```

## Step 3: Compute heating numbers (HEATR)

KERMA (Kinetic Energy Release in MAterial) coefficients give the energy
deposited per reaction. NJOY.jl provides pure-function heating kernels for
each reaction type.

```julia
# Q-values for U-238 reactions (eV)
Q_values = Dict{Int,Float64}(
    102 => 4.806e6,   # radiative capture Q
    18  => 197.0e6,   # fission Q (approximate total)
)

kerma = compute_kerma(pendf_300K;
    awr=awr,
    Z=92,               # atomic number for damage calculation
    Q_values=Q_values,
    E_d=40.0             # displacement energy for uranium [eV]
)
```

The `KERMAResult` contains per-energy arrays:

```julia
println("Total KERMA at 1 eV: ", kerma.total_kerma[10], " eV-barn")
println("Damage energy at 1 eV: ", kerma.damage_energy[10], " eV-barn")

# Verify the sum rule: total = elastic + capture + fission + inelastic
@assert verify_kerma_sum_rule(kerma)
```

Individual heating functions are also available for direct use:

```julia
# Elastic heating at 1 MeV for U-238
h_el = elastic_heating(1.0e6, awr)  # eV per collision

# Capture heating
h_cap = capture_heating(1.0e6, 4.806e6, awr)

# Lindhard damage partition
lp = lindhard_params(92.0, 238.0, 92.0, 238.0)
E_dam = lindhard_damage(1000.0, lp; E_d=40.0)
```

## Step 4: Produce an ACE file for MCNP (ACER)

The final step converts the broadened PENDF data to ACE format. NJOY.jl
supports both a structured builder and a flat-table builder.

### Using the flat-table builder

```julia
# Build ACE table (eV -> MeV conversion is automatic)
ace = build_ace(pendf_300K;
    suffix="80c",
    awr=awr,
    temperature=300.0,
    comment="U-238 from NJOY.jl, 300K",
    date="03/21/2026"
)

# Write Type 1 (ASCII) ACE file
open("92238.80c", "w") do io
    write_ace(io, ace)
end
```

### Using the structured builder

```julia
# Build structured ACE neutron table
ace_table = build_ace_from_pendf(pendf_300K;
    suffix="80c",
    temp_kelvin=300.0,
    comment="U-238 from NJOY.jl"
)

# Write to file
open("92238.80c", "w") do io
    write_ace_table(io, ace_table)
end
```

### Verifying the ACE file

Check the contents of the generated ACE table:

```julia
# Number of energy points
println("NES = ", length(ace_table.energy_grid))

# Energy range (in MeV, ACE convention)
println("E_min = ", ace_table.energy_grid[1], " MeV")
println("E_max = ", ace_table.energy_grid[end], " MeV")

# Number of non-elastic reactions
println("NTR = ", length(ace_table.reactions))
```

## Step 5 (optional): Write PENDF output

If you want a PENDF tape for downstream processing:

```julia
open("u238.pendf", "w") do io
    write_pendf(io, pendf_300K; temperature=300.0, err=0.001)
end
```

Or with the original ENDF file as a source for MF2 copying:

```julia
open("u238.pendf", "w") do io
    open(endf_file, "r") do endf_io
        write_pendf(io, pendf_300K;
            endf_source=endf_io,
            temperature=300.0)
    end
end
```

## Complete script

Putting it all together:

```julia
using NJOY

endf_file = "n-092_U_238.endf"
awr = 236.0058
T = 300.0

# RECONR: reconstruct pointwise cross sections
pendf = reconstruct(endf_file; err=0.001)

# BROADR: Doppler-broaden to target temperature
pendf_T = doppler_broaden(pendf, T; awr=awr, tol=0.001)

# HEATR: compute KERMA coefficients
kerma = compute_kerma(pendf_T;
    awr=awr, Z=92,
    Q_values=Dict(102 => 4.806e6, 18 => 197.0e6),
    E_d=40.0)
@assert verify_kerma_sum_rule(kerma)

# ACER: produce ACE file for MCNP
ace = build_ace(pendf_T; suffix="80c", awr=awr, temperature=T,
                comment="U-238 processed by NJOY.jl")
open("92238.80c", "w") do io
    write_ace(io, ace)
end

println("Done. ACE file written to 92238.80c")
println("  Grid points: ", length(pendf_T.energies))
println("  Reactions: ", pendf_T.mt_list)
```

## Additional processing modules

Beyond the basic RECONR -> BROADR -> HEATR -> ACER chain shown above,
NJOY.jl provides several additional processing capabilities:

### Thermal scattering (THERMR)

```julia
# Replace elastic channel with free-gas thermal scattering below 10 eV
pendf_thermal = compute_thermal(pendf_T, T, awr;
    model=:free_gas, emax=10.0)
```

### Multigroup averaging (GROUPR)

```julia
# Collapse to LANL 30-group structure with 1/E weighting
mg = group_average(
    pendf_T.energies,
    pendf_T.cross_sections,
    pendf_T.mt_list,
    LANL_30;
    weight_fn=inv_e_weight
)
```

### Unresolved resonance self-shielding (UNRESR)

```julia
# Bondarenko self-shielded cross sections at background dilutions
sigma0_values = [1e10, 1e4, 1e3, 1e2, 1e1]
result = bondarenko_xs(model, E, T, sigma0_values)
```

### Covariance processing (ERRORR)

```julia
# Process MF33 covariance data into multigroup matrices
cov_matrices = process_covariance(endf_file, LANL_30)
```

### Tape management (MODER)

```julia
# Scan a tape for its contents
open(endf_file, "r") do io
    dir = read_tape_directory(io)
    println("Materials: ", materials(dir))
end

# Extract one material
open(endf_file, "r") do io
    text = extract_material(io, 9237)
end
```
