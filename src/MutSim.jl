module MutSim

struct seq_mut_results
    position::Vector{Integer}
    Ref::Vector{String}
    Var::Vector{String}
    n_mut::Integer
end

include("sequences.jl")
include("JC69.jl")
include("cell_mutation.jl")
include("cell_proliferation.jl")
include("cancer_cell_sim.jl")
include("point_mutate_sim.jl")


export
    read_fa,
    get_seq,
    get_description,
    create_JC69,
    mut_simulation,
    cell_mutation,
    print_num_mut_per_gene,
    sim_cell_mut_branching,
    simulate,
    create_animation,
    simulate_dna_mutations
end