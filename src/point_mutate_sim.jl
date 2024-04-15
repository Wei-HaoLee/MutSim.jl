using Random
using Distributions

# Function to simulate DNA mutations
function simulate_dna_mutations(L::Float64, mutation_rate::Float64, mutated_positions::Vector{Int64})
    num_mutations = rand(Poisson(mutation_rate * L))  # Sample the number of mutations
    mutation_positions = rand(1:L, num_mutations)

    for mp in mutation_positions
        push!(mutated_positions, floor(Int, mp))
    end

    unique(mutated_positions)
end