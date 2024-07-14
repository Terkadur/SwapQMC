include("./qmc_pq2.jl")

const Lx, Ly = 4, 3
const T = hopping_matrix_Hubbard_2d(Lx, Ly, 1.0)

const U = 4.0
@show U

const filling = 0.75
const N_up = convert(Int64, Lx * Ly * filling)
const N_dn = N_up
@show filling

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
    18.0, 180,
    # data type of the system
    sys_type=ComplexF64,
    # if use charge decomposition
    useChargeHST=true,
    # if use first-order Trotteriaztion
    useFirstOrderTrotter=false
)

# 1.5s/warmup
# 15s/measure
const qmc = QMC(
    system,
    # number of warm-ups, samples and measurement interval
    16, 32, 10,
    # stablization and update interval
    10, 10,
    # if force spin symmetry
    forceSymmetry=true,
    # debugging flag
    saveRatio=false
)

const φ₀_up = trial_wf_free(system, 1, T)
const φ₀ = [φ₀_up, copy(φ₀_up)]

const partitionFraction = 0.25
const partitionSize = convert(Int64, Lx * Ly * partitionFraction)
const Aidx = collect(1:partitionSize)

const extsys = ExtendedSystem(system, Aidx, subsysOrdering=false)

path = "./tarek_data/wtf"
filename = "Pq2_N$(sum(system.N))_U$(system.U)_beta$(system.β)_seed$(seed)_smaller_use.jld"

@time run_replica_sampling_gs(extsys, qmc, φ₀, path, filename)
