using MutSim

path = "data/sequence.fa"
records = read_fa(path)


# single cell mutations
model = create_JC69(α=1e-5)
variant = cell_mutation(model, records)
print_num_mut_per_gene(variant)

# second n_generation cell mutations
second_variant = cell_mutation(model, records, variant)
print_num_mut_per_gene(second_variant)


# cell generation mutation
sim_cell_mut_branching(model, records, 2)