module EspressoDFT

using LinearAlgebra

export Crystal, KSModel, PlaneWaveBasis, SCFOptions, QEInput, AtomicDisplacement,
       read_qe_input, run_qe_input, ground_state, energy, forces, stress,
       density, eigenvalues, occupations, response, dynamical_matrix,
       phonon_modes, born_effective_charges, dielectric_tensor

struct Crystal
    lattice
    species
    positions
    masses
end

function Crystal(lattice, species, positions; masses, positions_are_fractional=true)
    fractional = positions_are_fractional ? positions : lattice \ positions
    Crystal(copy(lattice), copy(species), mod.(fractional, 1), copy(masses))
end

struct KSModel
    crystal
    electron_count
    xc
end
KSModel(crystal; pseudopotentials, xc=:pbe, charge=0, spin=:unpolarized) =
    KSModel(crystal, 2length(crystal.species), xc)

struct PlaneWaveBasis
    model
    Ecut
    kpoints
    kweights
    G_vectors
    fft_size
end
function PlaneWaveBasis(model; Ecut, kgrid, fft_size=nothing)
    nk = prod(kgrid)
    PlaneWaveBasis(model, Ecut, fill((0.0, 0.0, 0.0), nk), fill(1 / nk, nk),
                   [NTuple{3,Int}[] for _ in 1:nk], something(fft_size, (4, 4, 4)))
end

Base.@kwdef struct SCFOptions
    energy_tolerance::Float64 = 1e-10
    density_tolerance::Float64 = 1e-8
    maxiter::Int = 100
    extra_bands::Int = 4
end

struct QEInput
    model
    basis
    options
end

struct AtomicDisplacement
    atom::Int
    direction::Int
    q::NTuple{3,Float64}
end
AtomicDisplacement(atom, direction, q) =
    AtomicDisplacement(Int(atom), Int(direction), Tuple(Float64.(q)))

struct DummyGroundState
    basis
end

ground_state(basis; options=SCFOptions()) = DummyGroundState(basis)
energy(gs) = 0.0
forces(gs) = zeros(3, length(gs.basis.model.crystal.species))
stress(gs) = zeros(3, 3)
density(gs) = (values=zeros(gs.basis.fft_size), cell_volume=abs(det(gs.basis.model.crystal.lattice)))
eigenvalues(gs) = [zeros(4) for _ in gs.basis.kpoints]
occupations(gs) = [zeros(4) for _ in gs.basis.kpoints]
response(gs, perturbation; tolerance=1e-8, maxiter=200) =
    (delta_density=zeros(ComplexF64, gs.basis.fft_size), residual_norm=0.0, converged=true)
dynamical_matrix(gs, q; tolerance=1e-8, maxiter=200) =
    zeros(ComplexF64, 3length(gs.basis.model.crystal.species),
          3length(gs.basis.model.crystal.species))
function phonon_modes(gs, q; tolerance=1e-8, maxiter=200)
    n = 3length(gs.basis.model.crystal.species)
    (frequencies=zeros(n), eigenvectors=Matrix{Float64}(I, n, n))
end
born_effective_charges(gs; tolerance=1e-8, maxiter=200) =
    zeros(length(gs.basis.model.crystal.species), 3, 3)
dielectric_tensor(gs; tolerance=1e-8, maxiter=200) = Matrix{Float64}(I, 3, 3)

read_qe_input(path_or_io) = error("sentinel parser is intentionally incomplete")
run_qe_input(input) = input isa QEInput ? ground_state(input.basis) :
                      run_qe_input(read_qe_input(input))

end
