using MutSim

path = "data/sequence.fa"
records = read_fa(path)
sequences = get_seq(records)
seq = sequences[1]

model = create_JC69(α=1e-5)
foreach(x -> mut_simulation(model, x), sequences)
res = mut_simulation(model, seq)