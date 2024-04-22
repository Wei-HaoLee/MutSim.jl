using Plots
using DataStructures


# Initialize the grid with all normal cells
function initialize_grid(size)
    ones(Int, size, size)  # All cells start as normal
end

# Combined function to find neighbors
function find_neighbors(grid, x, y, condition)
    neighbors = []
    for dx in -1:1, dy in -1:1
        nx, ny = x + dx, y + dy
        if nx >= 1 && ny >= 1 && nx <= size(grid, 1) && ny <= size(grid, 2) && condition(grid[nx, ny])
            push!(neighbors, (nx, ny))
        end
    end
    neighbors
end

# Apply rules for tumor growth with parameters
using Random

# Precompute random choices to reduce the number of calls to the RNG
function random_choice(choices)
    return choices[rand(1:end)]
end


function apply_rules(grid, mutation_profiles, normal_mutation_rate, cancer_mutation_rate, genome_size, 
                     takeover_probability, death_rate, mutation_per_bp)
    grid_size = size(grid, 1)
    new_grid = copy(grid)  # Create a new grid to avoid modifying the grid in place

    # Process the grid in a single pass
    for x in 1:grid_size, y in 1:grid_size
        
        cell_type = grid[x, y]
        cell_index = (x-1) * grid_size + y
        
        if cell_type == 1 # normal cell
            if rand() < death_rate
                new_grid[x, y] = 0
            else
                empty_neighbors = find_neighbors(grid, x, y, v -> v == 0)

                if !isempty(empty_neighbors)
                    nx, ny = random_choice(empty_neighbors)
                    mutation_profiles[cell_index] = simulate_dna_mutations(genome_size, normal_mutation_rate, mutation_profiles[cell_index])

                    if length(mutation_profiles[cell_index]) >= mutation_per_bp * genome_size
                        new_grid[x, y] = 2
                    end

                    new_grid[nx, ny] = new_grid[x, y]
                    mutation_profiles[(nx-1) * grid_size + ny] = copy(mutation_profiles[cell_index])
                end
            end
        elseif cell_type == 2 # cancer cell
            empty_neighbors = find_neighbors(grid, x, y, v -> v == 0)
            if !isempty(empty_neighbors)
                nx, ny = random_choice(empty_neighbors)
                new_grid[nx, ny] = 2

                mutation_profiles[cell_index] = simulate_dna_mutations(genome_size, cancer_mutation_rate, mutation_profiles[cell_index])
                mutation_profiles[(nx-1) * grid_size + ny] = copy(mutation_profiles[cell_index])

            elseif rand() < takeover_probability
                normal_neighbors = find_neighbors(grid, x, y, v -> v == 1)
                if !isempty(normal_neighbors)
                    nx, ny = random_choice(normal_neighbors)
                    new_grid[nx, ny] = 2

                    mutation_profiles[cell_index] = simulate_dna_mutations(genome_size, cancer_mutation_rate, mutation_profiles[cell_index])
                    mutation_profiles[(nx-1) * grid_size + ny] = copy(mutation_profiles[cell_index])
                end
            end
        end
    end

    return new_grid
end



# Visualization of the grid
function plot_grid(grid)
    heatmap(grid, c=cgrad(["#d8e9ef", "#79a8a9", "#fc9d9a"], categorical=true), clim=(0, 2), size=(600, 600), axis=false, legend=false)
end

# Main simulation function
function simulate(grid_size, steps, normal_mutation_rate, cancer_mutation_rate, genome_size,
                  takeover_probability, death_rate, mutation_per_bp)
    grids = [initialize_grid(grid_size)]
    mutation_profiles = [Vector{Int64}() for _ in 1:grid_size^2]

    for i in 1:steps
        grid = apply_rules(grids[end], mutation_profiles,
                                              normal_mutation_rate, cancer_mutation_rate, genome_size,
                                              takeover_probability, death_rate,
                                              mutation_per_bp)
        
        push!(grids, grid)
    end

    grids, mutation_profiles
end

# Main simulation function
function simulate_cancer(grid_size, steps, normal_mutation_rate, cancer_mutation_rate, genome_size,
    takeover_probability, death_rate, mutation_per_bp, timepoint)
    grids = [initialize_grid(grid_size)]
    mutation_profiles = [Vector{Int64}() for _ in 1:grid_size^2]

    for i in 1:steps
        grid = apply_rules(grids[end], mutation_profiles,
                                        normal_mutation_rate, cancer_mutation_rate, genome_size,
                                        takeover_probability, death_rate,
                                        mutation_per_bp)

        if i == timepoint
            grid[rand(1:grid_size), rand(1:grid_size)] = 2
        end
        push!(grids, grid)
    end

    grids, mutation_profiles
end

function simulate_continue(grid, mutation_profiles, steps, normal_mutation_rate, cancer_mutation_rate, 
                           genome_size, takeover_probability, death_rate, mutation_per_bp)
    grids = [grid]

    for i in 1:steps
        grid = apply_rules(grids[end], mutation_profiles,
                                              normal_mutation_rate, cancer_mutation_rate, genome_size,
                                              takeover_probability, death_rate,
                                              mutation_per_bp)
        push!(grids, grid)
    end

    grids, mutation_profiles
end

function get_mutation_frequency(mutation_profiles)
    pos_count = OrderedDict()
    for cell_mutation in mutation_profiles
        for pos in cell_mutation
            if haskey(pos_count, pos)
                pos_count[pos] += 1
            else
                pos_count[pos] = 1
            end
        end
    end

    n_mut_freq = OrderedDict()
    for n_mut in values(pos_count)
        if haskey(n_mut_freq, n_mut)
            n_mut_freq[n_mut] += 1
        else
            n_mut_freq[n_mut] = 1
        end
    end
    n_mut_freq
end

function create_animation(grids, filename="cancer-cell-anim.gif", fps=10)
    anim = @animate for grid in grids
        plot_grid(grid)
    end
    gif(anim, filename, fps=fps)
end