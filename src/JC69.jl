using Distributions
using LinearAlgebra

struct JC69
    α::Float64
    transition_matrix::Matrix{Float64}
    nucleotides::Vector{String}
end

function create_JC69(; α=1e-8)
    if 3α >= 1
        ArgumentError("α must be less than 1/3")
    end

    trans_mat = fill(α, 4, 4)
    trans_mat[diagind(trans_mat)] .= 1-3α

    JC69(α, trans_mat, ["A", "C", "G", "T"])
end

function mut_simulation(model::JC69, seq::String; verbose=false)

    total_mutations = 0
    pos = Vector{Int}()
    refs = Vector{String}()
    vars = Vector{String}()

    for (position, nucleo) in enumerate(seq)
        nucleo = string(nucleo)
        nucleo_idx = findfirst(x -> x == nucleo, model.nucleotides)
        mut_prob = model.transition_matrix[:, nucleo_idx]
        mut_dist = Categorical(mut_prob)
        new_nucleo = model.nucleotides[rand(mut_dist)]

        if nucleo != new_nucleo
            push!(pos, position)
            push!(refs, nucleo)
            push!(vars, new_nucleo)

            total_mutations += 1
            if verbose
                println("Position: $position,\tRef: $nucleo,\tVar: $new_nucleo")
            end
        end
    end

    if verbose
        println("Total mutations: $total_mutations\n")
    end

    seq_mut_results(pos, refs, vars, total_mutations)
end