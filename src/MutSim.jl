module MutSim

include("sequences.jl")
include("JC69.jl")
include("cell_mutation.jl")
include("cell_proliferation.jl")

struct seq_mut_results
    position::Vector{Int64}
    Ref::Vector{String}
    Var::Vector{String}
end

export
    read_fa,
    read_fq,
    get_seq,
    get_description,
    create_JC69,
    mut_simulation
end