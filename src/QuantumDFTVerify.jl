module QuantumDFTVerify

using LinearAlgebra
using PseudoPotentialData
using QuantumEspresso_jll
using SHA
using TOML

export SPEC_ID, CM1_PER_ATOMIC_FREQUENCY, QEFixture, si_fixture, nacl_fixture,
       pseudopotential_paths, pseudopotential_hashes, write_pw_input,
       run_oracle, check_reference, load_reference

const SPEC_ID = "quantumdft-v0.1-qe7.5-2026-07-21"
const QE_VERSION = "7.5.0+0"
const AMU_TO_ELECTRON_MASS = 1822.888486209
const CM1_PER_ATOMIC_FREQUENCY = 219474.63136320

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

    Dict(
        "energy_ha" => energy_ha,
        "forces_ha_per_bohr" => [collect(forces[:, atom]) for atom in axes(forces, 2)],
        "stress_ha_per_bohr3" => [collect(stress[row, :]) for row in axes(stress, 1)],
    )
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
    mktempdir(prefix="quantumdft-qe75-") do workdir
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

function oracle_atol(path::AbstractString)
    occursin("frequencies_cm1", path) && return 5e-2
    occursin("born_charges", path) && return 5e-5
    occursin("dielectric", path) && return 1e-6
    occursin("dynamical_", path) && return 5e-10
    occursin("energy_ha", path) && return 5e-8
    occursin("forces_", path) && return 5e-8
    occursin("stress_", path) && return 5e-8
    0.0
end

function compare_reference(got, expected, path::AbstractString)
    if expected isa AbstractDict
        for key in keys(expected)
            haskey(got, key) || error("oracle output is missing $path.$key")
            compare_reference(got[key], expected[key], "$path.$key")
        end
    elseif expected isa AbstractVector
        length(got) == length(expected) ||
            error("oracle shape drift at $path: $(length(got)) != $(length(expected))")
        for index in eachindex(expected)
            compare_reference(got[index], expected[index], "$path[$index]")
        end
    elseif expected isa Number
        difference = abs(got - expected)
        tolerance = oracle_atol(path)
        difference <= tolerance ||
            error("oracle drift at $path: |$got - $expected| = $difference > $tolerance")
    else
        got == expected || error("oracle metadata drift at $path: $got != $expected")
    end
    true
end

function check_reference(reference::AbstractDict, generated::AbstractDict)
    reference["spec_id"] == SPEC_ID || error("reference specification mismatch")
    reference["qe_version"] == QE_VERSION || error("reference QE version mismatch")
    for fixture in (si_fixture(), nacl_fixture())
        expected = reference["fixtures"][fixture.name]
        got = generated["fixtures"][fixture.name]
        compare_reference(got, expected, "fixtures.$(fixture.name)")
    end
    true
end

end
