module EspressoDFTVerify

using LinearAlgebra
using PseudoPotentialData
using QuantumEspresso_jll
using SHA
using TOML

export SPEC_ID, CM1_PER_ATOMIC_FREQUENCY, BOHR_TO_ANGSTROM, QEFixture,
       si_fixture, nacl_fixture, aln_fixture, oracle_fixtures,
       pseudopotential_paths, pseudopotential_hashes, write_pw_input,
       run_oracle, check_reference, load_reference

const SPEC_ID = "espressodft-v0.2-qe7.5-2026-07-21"
const QE_VERSION = "7.5.0+0"
const AMU_TO_ELECTRON_MASS = 1822.888486209
const CM1_PER_ATOMIC_FREQUENCY = 219474.63136320
const BOHR_TO_ANGSTROM = 0.529177210903
const EV_TO_HARTREE = 1 / 27.211386245988
const DENSITY_G_SAMPLES = [
    (0, 0, 0), (1, 0, 0), (0, 1, 0), (0, 0, 1),
    (1, 1, 0), (1, 0, 1), (0, 1, 1), (1, 1, 1),
]

struct QEFixture
    name::String
    family::String
    lattice_bohr::Matrix{Float64}      # vectors are columns
    species::Vector{Symbol}
    masses_amu::Vector{Float64}
    positions_fractional::Matrix{Float64}
    ecut_ry::Float64
    kgrid::NTuple{3,Int}
    input_dft::String
    q_reduced::Vector{NTuple{3,Float64}}
end

function primitive_fcc(a::Real)
    a = Float64(a)
    [0.0 a / 2 a / 2; a / 2 0.0 a / 2; a / 2 a / 2 0.0]
end

si_fixture() = QEFixture(
    "si-private",
    "dojo.nc.sr.lda.v0_4_1.standard.upf",
    primitive_fcc(10.20),
    [:Si, :Si],
    [28.0855, 28.0855],
    [0.0 0.25; 0.0 0.25; 0.0 0.25],
    20.0,
    (4, 4, 4),
    "LDA",
    [(0.0, 0.0, 0.0), (0.25, 0.25, 0.25)],
)

nacl_fixture() = QEFixture(
    "nacl-private",
    "dojo.nc.sr.pbe.v0_4_1.standard.upf",
    primitive_fcc(10.60),
    [:Na, :Cl],
    [22.98976928, 35.45],
    [0.0 0.5; 0.0 0.5; 0.0 0.5],
    60.0,
    (4, 4, 4),
    "PBE",
    [(0.0, 0.0, 0.0)],
)

function aln_fixture()
    a, c, u = 5.90, 9.65, 0.382
    lattice = [a -a / 2 0.0; 0.0 sqrt(3) * a / 2 0.0; 0.0 0.0 c]
    QEFixture(
        "aln-heldout",
        "dojo.nc.sr.pbe.v0_4_1.standard.upf",
        lattice,
        [:Al, :Al, :N, :N],
        [26.9815385, 26.9815385, 14.0067, 14.0067],
        [0.0 2 / 3 0.0 2 / 3;
         0.0 1 / 3 0.0 1 / 3;
         0.0 1 / 2 u 1 / 2 + u],
        60.0,
        (4, 4, 2),
        "PBE",
        [(0.0, 0.0, 0.0), (0.25, 0.0, 0.0)],
    )
end

oracle_fixtures() = (si_fixture(), nacl_fixture(), aln_fixture())

function pseudopotential_paths(fixture::QEFixture)
    family = PseudoPotentialData.PseudoFamily(fixture.family)
    unique_species = unique(fixture.species)
    Dict(element => PseudoPotentialData.pseudofile(family, element)
         for element in unique_species)
end

pseudopotential_hashes(fixture::QEFixture) =
    Dict(String(element) => bytes2hex(open(sha256, path))
         for (element, path) in pseudopotential_paths(fixture))

fmt(x::Real) = string(Float64(x))

function write_pw_input(path::AbstractString, fixture::QEFixture,
                        pseudo_paths::AbstractDict, outdir::AbstractString)
    mkpath(outdir)
    unique_species = unique(fixture.species)
    open(path, "w") do io
        println(io, "&CONTROL")
        println(io, "  calculation = 'scf'")
        println(io, "  prefix = '", fixture.name, "'")
        println(io, "  pseudo_dir = '", dirname(first(values(pseudo_paths))), "'")
        println(io, "  outdir = '", outdir, "'")
        println(io, "  tprnfor = .true.")
        println(io, "  tstress = .true.")
        println(io, "/")
        println(io, "&SYSTEM")
        println(io, "  ibrav = 0")
        println(io, "  nat = ", length(fixture.species))
        println(io, "  ntyp = ", length(unique_species))
        println(io, "  ecutwfc = ", fmt(fixture.ecut_ry))
        println(io, "  occupations = 'fixed'")
        println(io, "  nspin = 1")
        println(io, "  input_dft = '", fixture.input_dft, "'")
        println(io, "/")
        println(io, "&ELECTRONS")
        println(io, "  conv_thr = 1.0d-12")
        println(io, "  electron_maxstep = 120")
        println(io, "  mixing_beta = 0.5")
        println(io, "/")
        println(io, "ATOMIC_SPECIES")
        for element in unique_species
            atom = findfirst(==(element), fixture.species)
            println(io, element, " ", fmt(fixture.masses_amu[atom]), " ",
                    basename(pseudo_paths[element]))
        end
        println(io, "CELL_PARAMETERS bohr")
        for vector in eachcol(fixture.lattice_bohr)
            println(io, join(fmt.(vector), " "))
        end
        println(io, "ATOMIC_POSITIONS crystal")
        for atom in eachindex(fixture.species)
            println(io, fixture.species[atom], " ",
                    join(fmt.(fixture.positions_fractional[:, atom]), " "))
        end
        println(io, "K_POINTS automatic")
        println(io, join(fixture.kgrid, " "), " 0 0 0")
    end
    path
end

function q_to_qe_cartesian(fixture::QEFixture, q::NTuple{3,<:Real})
    alat = norm(fixture.lattice_bohr[:, 1])
    reduced_to_qe = alat * inv(fixture.lattice_bohr)'
    Tuple(reduced_to_qe * collect(Float64, q))
end

function write_ph_input(path::AbstractString, fixture::QEFixture,
                        outdir::AbstractString, dynpath::AbstractString,
                        q::NTuple{3,<:Real})
    qe_q = q_to_qe_cartesian(fixture, q)
    gamma = all(iszero, q)
    open(path, "w") do io
        println(io, "&INPUTPH")
        println(io, "  tr2_ph = 1.0d-14")
        println(io, "  prefix = '", fixture.name, "'")
        println(io, "  outdir = '", outdir, "'")
        println(io, "  fildyn = '", dynpath, "'")
        println(io, "  epsil = ", gamma ? ".true." : ".false.")
        println(io, "  ldisp = .false.")
        println(io, "/")
        println(io, join(fmt.(qe_q), " "))
    end
    path
end

function write_pp_input(path::AbstractString, fixture::QEFixture,
                        outdir::AbstractString, cube_path::AbstractString)
    plot_path = joinpath(dirname(cube_path), "density.plot")
    open(path, "w") do io
        println(io, "&INPUTPP")
        println(io, "  prefix = '", fixture.name, "'")
        println(io, "  outdir = '", outdir, "'")
        println(io, "  filplot = '", plot_path, "'")
        println(io, "  plot_num = 0")
        println(io, "/")
        println(io, "&PLOT")
        println(io, "  nfile = 1")
        println(io, "  filepp(1) = '", plot_path, "'")
        println(io, "  weight(1) = 1.0")
        println(io, "  iflag = 3")
        println(io, "  output_format = 6")
        println(io, "  fileout = '", cube_path, "'")
        println(io, "/")
    end
    path
end

parse_fortran_float(value::AbstractString) = parse(Float64, replace(value, r"[dD]" => "E"))
numbers(line::AbstractString) = parse_fortran_float.(getproperty.(collect(eachmatch(
    r"[-+]?(?:\d+\.\d*|\.\d+|\d+)(?:[EeDd][-+]?\d+)?", line)), :match))

function parse_pw_output(text::AbstractString)
    energy_matches = collect(eachmatch(r"!\s+total energy\s+=\s+([-+0-9.EeDd]+)\s+Ry", text))
    isempty(energy_matches) && error("QE pw.x output has no final total energy")
    energy_ha = parse_fortran_float(energy_matches[end].captures[1]) / 2

    force_rows = [numbers(m.match)[end-2:end] ./ 2 for m in eachmatch(
        r"(?m)^\s*atom\s+\d+\s+type\s+\d+\s+force\s*=.*$", text)]
    forces = isempty(force_rows) ? Matrix{Float64}(undef, 3, 0) :
             reduce(hcat, force_rows)

    stress_marker = findlast("total   stress  (Ry/bohr**3)", text)
    stress_marker === nothing && error("QE pw.x output has no stress block")
    tail = text[first(stress_marker):end]
    stress_lines = split(tail, '\n')
    data_start = findfirst(line -> length(numbers(line)) >= 6, stress_lines)
    data_start === nothing && error("QE stress block is malformed")
    stress_rows = [numbers(stress_lines[data_start + i])[1:3] ./ 2 for i in 0:2]
    stress = reduce(hcat, stress_rows)'

    gamma_bands = Float64[]
    lines = split(text, '\n')
    for (index, line) in pairs(lines)
        marker = match(r"k\s*=\s*([-+0-9.]+)\s+([-+0-9.]+)\s+([-+0-9.]+).*bands \(ev\):", line)
        marker === nothing && continue
        k = parse.(Float64, marker.captures)
        norm(k) <= 1e-10 || continue
        started = false
        for band_line in lines[index + 1:end]
            values = numbers(band_line)
            if isempty(values)
                started && break
                continue
            end
            occursin(r"[A-Za-z]", replace(band_line, r"[Ee][+-]?\d+" => "")) && break
            append!(gamma_bands, values .* EV_TO_HARTREE)
            started = true
        end
        break
    end
    isempty(gamma_bands) && error("QE pw.x output has no Gamma-point bands")

    Dict(
        "energy_ha" => energy_ha,
        "forces_ha_per_bohr" => [collect(forces[:, atom]) for atom in axes(forces, 2)],
        "stress_ha_per_bohr3" => [collect(stress[row, :]) for row in axes(stress, 1)],
        "gamma_eigenvalues_ha" => gamma_bands,
    )
end

function parse_cube_density(text::AbstractString)
    lines = split(text, '\n')
    length(lines) >= 7 || error("QE pp.x cube output is truncated")
    origin = numbers(lines[3])
    length(origin) >= 4 || error("QE cube atom/origin row is malformed")
    nat = abs(round(Int, origin[1]))
    dims = ntuple(axis -> abs(round(Int, numbers(lines[3 + axis])[1])), 3)
    first_value_line = 7 + nat
    values = Float64[]
    for line in lines[first_value_line:end]
        append!(values, numbers(line))
    end
    prod(dims) == length(values) ||
        error("QE cube density size mismatch: $(length(values)) != $(prod(dims))")
    permutedims(reshape(values, dims[3], dims[2], dims[1]), (3, 2, 1))
end

function density_coefficients(values::AbstractArray{<:Real,3})
    dims = size(values)
    coefficients = Dict{String,Any}()
    for G in DENSITY_G_SAMPLES
        coefficient = zero(ComplexF64)
        for i in axes(values, 1), j in axes(values, 2), k in axes(values, 3)
            reduced = ((i - 1) / dims[1], (j - 1) / dims[2], (k - 1) / dims[3])
            phase = -2pi * (G[1] * reduced[1] + G[2] * reduced[2] + G[3] * reduced[3])
            coefficient += values[i, j, k] * cis(phase)
        end
        coefficient /= length(values)
        coefficients[join(G, ",")] = [real(coefficient), imag(coefficient)]
    end
    coefficients
end

function parse_frequencies_cm1(text::AbstractString)
    [parse_fortran_float(m.captures[1]) for m in eachmatch(
        r"freq\s*\(\s*\d+\)\s*=.*?=\s*([-+0-9.EeDd]+)\s*\[cm-1\]", text)]
end

function parse_dielectric(text::AbstractString)
    marker = findlast("Dielectric constant in cartesian axis", text)
    marker === nothing && return nothing
    tail = split(text[first(marker):end], '\n')
    rows = Vector{Vector{Float64}}()
    for line in tail[2:end]
        vals = numbers(line)
        if length(vals) == 3
            push!(rows, vals)
            length(rows) == 3 && break
        end
    end
    length(rows) == 3 || error("QE dielectric block is malformed")
    rows
end

function parse_asr_born_charges(text::AbstractString, nat::Integer)
    marker = findlast("with asr applied:", text)
    marker === nothing && return nothing
    rows = Vector{Vector{Float64}}()
    for line in split(text[first(marker):end], '\n')
        occursin(r"E\*[xyz]", line) || continue
        vals = numbers(line)
        length(vals) >= 3 && push!(rows, vals[end-2:end])
        length(rows) == 3nat && break
    end
    length(rows) == 3nat || error("QE ASR Born-charge block is malformed")
    [[rows[3(atom - 1) + direction] for direction in 1:3] for atom in 1:nat]
end

function parse_first_force_constant(text::AbstractString, fixture::QEFixture)
    lines = split(text, '\n')
    start = findfirst(line -> occursin("Dynamical  Matrix in cartesian axes", line), lines)
    start === nothing && error("QE dynamical-matrix file has no matrix")
    nat = length(fixture.species)
    phi = zeros(ComplexF64, 3nat, 3nat)
    line = start + 1
    blocks = 0
    while line <= length(lines)
        occursin("Dielectric Tensor", lines[line]) && break
        block = match(r"^\s*(\d+)\s+(\d+)\s*$", lines[line])
        if block === nothing
            line += 1
            continue
        end
        atom_i, atom_j = parse.(Int, block.captures)
        (1 <= atom_i <= nat && 1 <= atom_j <= nat) || break
        for direction_i in 1:3
            line += 1
            vals = numbers(lines[line])
            length(vals) == 6 || error("malformed QE dynamical-matrix row")
            for direction_j in 1:3
                row = 3(atom_i - 1) + direction_i
                col = 3(atom_j - 1) + direction_j
                phi[row, col] = complex(vals[2direction_j - 1], vals[2direction_j])
            end
        end
        blocks += 1
        line += 1
        blocks == nat^2 && break
    end
    blocks == nat^2 || error("incomplete first QE dynamical matrix")

    masses = fixture.masses_amu .* AMU_TO_ELECTRON_MASS
    dynamical = similar(phi)
    for atom_i in 1:nat, atom_j in 1:nat, a in 1:3, b in 1:3
        row, col = 3(atom_i - 1) + a, 3(atom_j - 1) + b
        dynamical[row, col] = 0.5phi[row, col] / sqrt(masses[atom_i] * masses[atom_j])
    end
    dynamical
end

matrix_rows(matrix::AbstractMatrix) = [collect(matrix[row, :]) for row in axes(matrix, 1)]

function run_program(product::Function, arguments::Cmd, stdout_path::AbstractString,
                     stderr_path::AbstractString)
    product() do executable
        run(pipeline(`$executable $arguments`, stdout=stdout_path, stderr=stderr_path))
    end
end

function run_oracle(fixture::QEFixture)
    mktempdir(prefix="espressodft-qe75-") do workdir
        outdir = joinpath(workdir, "scratch")
        mkpath(outdir)
        pseudos = pseudopotential_paths(fixture)
        pw_input = write_pw_input(joinpath(workdir, "pw.in"), fixture, pseudos, outdir)
        pw_output = joinpath(workdir, "pw.out")
        run_program(QuantumEspresso_jll.pwscf, `-in $pw_input`, pw_output,
                    joinpath(workdir, "pw.err"))
        result = parse_pw_output(read(pw_output, String))
        result["pseudopotential_sha256"] = pseudopotential_hashes(fixture)
        result["family"] = fixture.family
        result["ecut_ry"] = fixture.ecut_ry
        result["kgrid"] = collect(fixture.kgrid)

        pp_input = joinpath(workdir, "pp.in")
        pp_output = joinpath(workdir, "pp.out")
        cube_path = joinpath(workdir, "density.cube")
        write_pp_input(pp_input, fixture, outdir, cube_path)
        run_program(QuantumEspresso_jll.pp, `-in $pp_input`, pp_output,
                    joinpath(workdir, "pp.err"))
        cube_density = parse_cube_density(read(cube_path, String))
        result["density_grid"] = collect(size(cube_density))
        result["density_fourier"] = density_coefficients(cube_density)
        result["density_electron_count"] =
            result["density_fourier"]["0,0,0"][1] * abs(det(fixture.lattice_bohr))

        phonons = Dict{String,Any}()
        for (index, q) in enumerate(fixture.q_reduced)
            ph_input = joinpath(workdir, "ph-$index.in")
            ph_output = joinpath(workdir, "ph-$index.out")
            dynpath = joinpath(workdir, "dyn-$index")
            write_ph_input(ph_input, fixture, outdir, dynpath, q)
            run_program(QuantumEspresso_jll.phonon, `-in $ph_input`, ph_output,
                        joinpath(workdir, "ph-$index.err"))
            ph_text = read(ph_output, String)
            dyn = parse_first_force_constant(read(dynpath, String), fixture)
            key = join(fmt.(q), ",")
            entry = Dict{String,Any}(
                "q_reduced" => collect(q),
                "frequencies_cm1_raw" => parse_frequencies_cm1(ph_text),
                "dynamical_real" => matrix_rows(real(dyn)),
                "dynamical_imag" => matrix_rows(imag(dyn)),
            )
            if all(iszero, q)
                entry["dielectric"] = parse_dielectric(ph_text)
                entry["born_charges_asr"] = parse_asr_born_charges(
                    ph_text, length(fixture.species))
            end
            phonons[key] = entry
        end
        result["phonons"] = phonons
        result
    end
end

load_reference(path::AbstractString=joinpath(@__DIR__, "..", "oracle", "reference.toml")) =
    TOML.parsefile(path)

function oracle_atol(path::AbstractString, expected=nothing)
    if occursin("frequencies_cm1", path)
        # Square-rooting tiny acoustic eigenvalues amplifies otherwise accepted
        # 1e-9-scale matrix noise. Protect the matrix separately and use a
        # still-sub-candidate 1 cm^-1 reproducibility gate only for soft modes.
        return expected isa Number && abs(expected) < 20 ? 1.0 : 2e-1
    end
    occursin("born_charges", path) && return 5e-5
    # QE 7.5 differs by about 1.92e-6 for the largest Si dielectric entry
    # between the pinned x86_64 Linux and aarch64 macOS artifacts. Keep the
    # reproducibility gate three orders tighter than the candidate tolerance.
    occursin("dielectric", path) && return 5e-6
    occursin("dynamical_", path) && return 3e-9
    occursin("density_fourier", path) && return 2e-5
    occursin("density_electron_count", path) && return 2e-6
    occursin("gamma_eigenvalues", path) && return 5e-6
    occursin("energy_ha", path) && return 5e-8
    occursin("forces_", path) && return 5e-8
    occursin("stress_", path) && return 5e-8
    0.0
end

function collect_reference_errors!(errors::Vector{String}, got, expected,
                                   path::AbstractString)
    if expected isa AbstractDict
        for key in keys(expected)
            if !haskey(got, key)
                push!(errors, "oracle output is missing $path.$key")
                continue
            end
            collect_reference_errors!(errors, got[key], expected[key], "$path.$key")
        end
    elseif expected isa AbstractVector
        if length(got) != length(expected)
            push!(errors,
                  "oracle shape drift at $path: $(length(got)) != $(length(expected))")
            return errors
        end
        for index in eachindex(expected)
            collect_reference_errors!(errors, got[index], expected[index],
                                      "$path[$index]")
        end
    elseif expected isa Number
        difference = abs(got - expected)
        tolerance = oracle_atol(path, expected)
        difference <= tolerance || push!(
            errors,
            "oracle drift at $path: |$got - $expected| = $difference > $tolerance",
        )
    else
        got == expected || push!(errors,
                                 "oracle metadata drift at $path: $got != $expected")
    end
    errors
end

function compare_reference(got, expected, path::AbstractString)
    errors = String[]
    collect_reference_errors!(errors, got, expected, path)
    isempty(errors) || error(join(errors, '\n'))
    true
end

function check_reference(reference::AbstractDict, generated::AbstractDict)
    reference["spec_id"] == SPEC_ID || error("reference specification mismatch")
    reference["qe_version"] == QE_VERSION || error("reference QE version mismatch")
    errors = String[]
    for fixture in oracle_fixtures()
        expected = reference["fixtures"][fixture.name]
        got = generated["fixtures"][fixture.name]
        collect_reference_errors!(errors, got, expected, "fixtures.$(fixture.name)")
    end
    isempty(errors) || error("oracle reference mismatch:\n" * join(errors, '\n'))
    true
end

end
