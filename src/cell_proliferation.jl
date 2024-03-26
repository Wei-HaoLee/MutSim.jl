using UUIDs
using Random

struct cell_node
    ancestor_cell::UUID
    cell::UUID
    nth_generation::Integer
    variant::Dict{String, seq_mut_results}
end

function sim_cell_mut_branching(model::Union{JC69}, 
    sequences::Vector{FASTX.FASTA.Record},
    n_generation::Integer)

    rng = MersenneTwister(612)
    
    # store each cell's variant and information
    cells_variant = Dict{UUID, cell_node}()
    
    # create an ancestor cell (root)
    ancestor_UUID = uuid1(rng)
    ancestor_variant = cell_mutation(model, sequences)
    cells_variant[ancestor_UUID] = cell_node(ancestor_UUID, ancestor_UUID, 1, ancestor_variant)
    ancestor_cells = [ancestor_UUID]

    print(ancestor_UUID)
    # start each generation
    for i in 2:n_generation
        for ac_UUID in ancestor_cells
            ac_variant = cells_variant[ac_UUID].variant

            # create two new cells from ancestor cell
            for _ in 1:2
                cell_UUID = uuid1(rng)
                cell_variant = cell_mutation(model, sequences, ac_variant)
                cells_variant[cell_UUID] = cell_node(ac_UUID, cell_UUID, i, cell_variant)

                # add new cell to ancestor cells
                push!(ancestor_cells, cell_UUID)
            end

            # remove their parents from ancestor cells
            filter!(e -> e ≠ ac_UUID, ancestor_cells)
        end
    end

    cells_variant
end