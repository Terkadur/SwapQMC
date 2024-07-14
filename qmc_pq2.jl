using JLD, SwapQMC

function run_regular_sampling_gs(
    extsys::ExtendedSystem{Sys}, qmc::QMC, φ₀::Vector{Wf},
    path::String, filename::String
) where {Sys<:Hubbard,Wf<:AbstractMatrix}

    system = extsys.system
    bins = qmc.measure_interval

    walker = HubbardWalker(system, qmc, φ₀)

    # measurements
    sampler = EtgSampler(extsys, qmc)

    # thermalization step
    sweep!(system, qmc, walker, loop_number=qmc.nwarmups)

    for i in 1:qmc.nsamples
        sweep!(system, qmc, walker, loop_number=bins)
        measure_Pn!(sampler, walker, forwardMeasurement=true)
    end

    # store the measurement
    jldopen("$(path)/$(filename)", "w") do file
        write(file, "Pn_up", sampler.Pn₊)
        write(file, "Pn_dn", sampler.Pn₋)
    end

    return nothing
end

function run_replica_sampling_gs(
    extsys::ExtendedSystem{Sys}, qmc::QMC, φ₀::Vector{Wf},
    path::String, filename::String
) where {Sys<:Hubbard,Wf<:AbstractMatrix}

    system = extsys.system

    walker1 = HubbardWalker(system, qmc, φ₀)
    walker2 = HubbardWalker(system, qmc, φ₀)

    replica = Replica(extsys, walker1, walker2)

    bins = qmc.measure_interval

    # measurements
    sampler = EtgSampler(extsys, qmc)

    for i in 1:qmc.nwarmups
        (i % 16 == 0) && println(i)
        if (i - 1) % 256 < 127
            sweep!(system, qmc, replica, walker1, 1, loop_number=1, jumpReplica=false)
        elseif (i - 1) % 256 == 127
            sweep!(system, qmc, replica, walker1, 1, loop_number=1, jumpReplica=true)
        elseif 127 < (i - 1) % 256 < 255
            sweep!(system, qmc, replica, walker2, 2, loop_number=1, jumpReplica=false)
        else
            sweep!(system, qmc, replica, walker2, 2, loop_number=1, jumpReplica=true)
        end
    end

    for i in 1:qmc.nsamples
        (i % 16 == 0) && println(i)
        if (i - 1) % 256 < 127
            sweep!(system, qmc, replica, walker1, 1, loop_number=bins, jumpReplica=false)
        elseif (i - 1) % 256 == 127
            sweep!(system, qmc, replica, walker1, 1, loop_number=bins, jumpReplica=true)
        elseif 127 < (i - 1) % 256 < 255
            sweep!(system, qmc, replica, walker2, 2, loop_number=bins, jumpReplica=false)
        else
            sweep!(system, qmc, replica, walker2, 2, loop_number=bins, jumpReplica=true)
        end

        measure_Pn2!(sampler, replica, forwardMeasurement=true)
    end
    
    # store the measurement
    jldopen("$(path)/$(filename)", "w") do file
        write(file, "Pn2_up", sampler.Pn₊)
        write(file, "Pn2_dn", sampler.Pn₋)
    end

    return nothing
end
