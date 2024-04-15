using Plots

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
function apply_rules(grid, mutation_profiles, normal_mutation_rate, cancer_mutation_rate, genome_size,
    takeover_probability, death_rate, proliferation_rate_normal, proliferation_rate_tumor)
    for x in 1:size(grid, 1), y in 1:size(grid, 2)
        cell_type = grid[x, y]
        cell_index = (x-1) * size(grid, 1) + y

            if cell_type == 1 || cell_type == 2  # Apply mutations only if the cell is normal or tumor
            mutation_profiles[cell_index] = simulate_dna_mutations(genome_size, cell_type == 1 ? normal_mutation_rate : cancer_mutation_rate, mutation_profiles[cell_index])
            end

                if cell_type == 1 && length(mutation_profiles[cell_index]) >= (0.1/1e+6) * genome_size  # Check if normal cell transitions to a tumor cell
                    grid[x, y] = 2
                elseif cell_type == 1 && rand() < death_rate  # Normal cell dies with a certain probability
                    grid[x, y] = 0
                else
                    empty_neighbors = find_neighbors(grid, x, y, v -> v == 0)
                if !isempty(empty_neighbors) && rand() < (cell_type == 1 ? proliferation_rate_normal : proliferation_rate_tumor)
                    nx, ny = empty_neighbors[rand(1:end)]
                    grid[nx, ny] = cell_type
                    mutation_profiles[(nx-1) * size(grid, 1) + ny] = copy(mutation_profiles[cell_index])
                end
                if cell_type == 2 && rand() < takeover_probability
                    normal_neighbors = find_neighbors(grid, x, y, v -> v == 1)
                if !isempty(normal_neighbors)
                    nx, ny = normal_neighbors[rand(1:end)]
                    grid[nx, ny] = 2
                    mutation_profiles[(nx-1) * size(grid, 1) + ny] = copy(mutation_profiles[cell_index])
                end
            end
        end
    end
end



# Visualization of the grid
function plot_grid(grid)
    heatmap(grid, c=cgrad(["#d8e9ef", "#79a8a9", "#fc9d9a"], categorical=true), clim=(0, 2), size=(600, 600), axis=false, legend=false)
end

# Main simulation function
function simulate(grid_size, steps, normal_mutation_rate, cancer_mutation_rate, genome_size,
                  takeover_probability, death_rate, proliferation_rate_normal, proliferation_rate_tumor, 
                  show_step=50)
    grids = [initialize_grid(grid_size)]
    mutation_profiles = [Vector{Int64}() for _ in 1:grid_size^2]

    for i in 1:steps
        if i % show_step == 0
            println("Current step: $i")
        end
        apply_rules(grids[end], mutation_profiles,
                    normal_mutation_rate, cancer_mutation_rate, genome_size,
                    takeover_probability, death_rate,
                    proliferation_rate_normal, proliferation_rate_tumor)
        push!(grids, copy(grids[end]))
    end

    grids, mutation_profiles
end

function simulate_continue(grid, mutation_profiles, steps, normal_mutation_rate, cancer_mutation_rate, genome_size,
    takeover_probability, death_rate, proliferation_rate_normal, proliferation_rate_tumor, show_step=50)
    grids = [grid]

    for i in 1:steps
        if i % show_step == 0
            println("Current step: $i")
        end
        
        apply_rules(grids[end], mutation_profiles,
            normal_mutation_rate, cancer_mutation_rate, genome_size,
            takeover_probability, death_rate,
            proliferation_rate_normal, proliferation_rate_tumor)
        push!(grids, copy(grids[end]))
    end

    grids, mutation_profiles
end

function create_animation(grids, filename="cancer-cell-anim.gif", fps=10)
    anim = @animate for grid in grids
        plot_grid(grid)
    end
    gif(anim, filename, fps=fps)
end