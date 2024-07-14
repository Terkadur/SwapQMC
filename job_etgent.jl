include("./qmc_etgent.jl")

const Lx, Ly = 2, 2
const T = hopping_matrix_Hubbard_2d(Lx, Ly, 1.0)

const U = 8.0
@show U

const filling = 0.75
const N_up = convert(Int64, Lx * Ly * filling)
const N_dn = N_up
@show filling

const Nₖ = 32
const λₖ = 16 / Nₖ
@show λₖ

seed = 1234
Random.seed!(seed)
@show seed

const system = GenericHubbard(
    # (Nx, Ny), (N_up, N_dn)
    (Lx, Ly, 1), (N_up, N_dn),
    # t, U
    T, U,
    # μ
    0.0,
    # β, L
    50.0, 500,
    # data type of the system
    sys_type=ComplexF64,
    # if use charge decomposition
    useChargeHST=true,
    # if use first-order Trotteriaztion
    useFirstOrderTrotter=false
)

# 5s/warmup
# 51s/measure
const qmc = QMC(
    system,
    # number of warm-ups, samples and measurement interval
    1, 1, 10,
    # stablization and update interval
    10, 10,
    # if force spin symmetry
    forceSymmetry=true,
    # debugging flag
    saveRatio=false
)

const φ₀_up = trial_wf_free(system, 1, T)
const φ₀ = [φ₀_up, copy(φ₀_up)]

const partitionFraction = 0.5
const partitionSize = convert(Int64, Lx * Ly * partitionFraction)
const Aidx = collect(1:partitionSize)

const extsys = ExtendedSystem(system, Aidx, subsysOrdering=false)

path = "./tarek_data/wtf"
filename = "EtgEnt_N$(sum(system.N))_U$(system.U)_lambda$(λₖ)_beta$(system.β)_seed$(seed)_smaller_dont.jld"
@time run_incremental_sampling_gs(extsys, qmc, φ₀, λₖ, Nₖ, path, filename)
