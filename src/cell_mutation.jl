function cell_mutation(model::Union{JC69}, 
    sequences::Vector{FASTX.FASTA.Record})::Dict{String, seq_mut_results}

    sequences_mutations = Dict{String, seq_mut_results}()
    for s in sequences
       seq_header =  description(s)

       sequences_mutations[seq_header] = seq_mut_results(Vector{Integer}(), Vector{String}(), Vector{String}(), 0)
    end

    return sequences_mutations
end

function seq_modify(seq::String, variant::seq_mut_results)
    # string to vector of characters
    seq_mut = collect(seq)

    for i in 1:variant.n_mut
        pos = variant.position[i]
        var = variant.Var[i]

        # var needs to be converted to chracter
        seq_mut[pos] = collect(var)[1]
    end

    return join(seq_mut)
end

function cell_mutation(model::Union{JC69}, 
    sequences::Vector{FASTX.FASTA.Record},
    variant::Dict{String, seq_mut_results})::Dict{String, seq_mut_results}
    
    sequences_mutations = Dict{String, seq_mut_results}()
    for s in sequences
        seq_header =  description(s)
        seq = sequence(String, s)
        seq_mut = seq_modify(seq, variant[seq_header])

        new_variant = mut_simulation(model, seq_mut)
        total_variant = seq_mut_results(
            [variant[seq_header].position; new_variant.position],
            [variant[seq_header].Ref; new_variant.Ref],
            [variant[seq_header].Var; new_variant.Var],
            variant[seq_header].n_mut + new_variant.n_mut,
        )
        sequences_mutations[seq_header] = total_variant
    end

    return sequences_mutations
end

function print_num_mut_per_gene(variant::Dict{String, seq_mut_results})
    for (k, v) in variant
        println("$k: $(v.n_mut) mutations")
    end
end