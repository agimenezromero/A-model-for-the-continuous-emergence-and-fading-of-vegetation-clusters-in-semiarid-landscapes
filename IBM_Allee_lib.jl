using Distances
using Distributions
using StatsBase

using DelimitedFiles

mutable struct Parameters{T<:Real}

    constants::Array{T} # γ_c, σ_c, γ_p, σ_p, σ_d, c, p, d, f

    total_rates::Array{T} # D, P (total Death, total Proliferation)

    L::T #System length

end

mutable struct Variables{T<:Real}

    P::Array{T} #2xN matrix for positions of particles

    dist::Array{T} #NxN matrix for distances between particles

    rates::Array{T} # D_squared, P, D (Death_squared, Proliferation, Movement, Death)    

    N_particles::Int #Number of particles

    switching_times::Array{T} #Array of switching times

end

# Implement periodic boundary conditions
function periodic(x, L)
    if x < -L / 2.0
        return x + L
    elseif x >= L / 2.0
        return x - L
    else
        return x
    end
end

# Competition kernel (Gaussian)
@views function competition_kernel(ξ_squared, params)

    return @. params.constants[1] * exp(-ξ_squared / (2.0 * params.constants[2]^2))

end

# Proliferation kernel (Gaussian)
@views function proliferation_kernel(ξ_squared, params)

    return @. params.constants[3] * exp(-ξ_squared / (2.0 * params.constants[4]^2))

end

# Initialization of variables (distances and rates from given initial positions)
@views function initialize_vars(vars, params)

    @inbounds for i in 1:vars.N_particles

        @inbounds for j in 1:vars.N_particles

            #Compute distance to other particles
            xi_x = periodic(vars.P[1, j] - vars.P[1, i], params.L)
            xi_y = periodic(vars.P[2, j] - vars.P[2, i], params.L)

            ξ_squared = xi_x^2 + xi_y^2

            #Fill distance matrix
            vars.dist[i, j] = sqrt(ξ_squared)

            if j ≠ i

                vars.rates[3, i] += competition_kernel(ξ_squared, params)
                vars.rates[2, i] += proliferation_kernel(ξ_squared, params)

            end

        end

        vars.rates[1, i] = params.constants[6] .+ vars.rates[3, i] .^ 2 #D_squared_n
        vars.rates[2, i] += params.constants[7]  #P_n
    end

end

# Recompute rates (when changing the kernels)
@views function recompute_rates(vars, params)

    @inbounds for i in 1:vars.N_particles

        @inbounds for j in 1:vars.N_particles

            if j ≠ i

                vars.rates[3, i] += competition_kernel(vars.dist[i, j]^2, params)
                vars.rates[2, i] += proliferation_kernel(vars.dist[i, j]^2, params)

            end

        end

        vars.rates[1, i] = params.constants[6] .+ vars.rates[3, i] .^ 2 #D_squared_n
        vars.rates[2, i] += params.constants[7]  #P_n
    end

end

# Implement offspring event (adds a new particle at a given position from the parent particle)
@views function add_particle(idx, vars, params)

    vars.N_particles += 1

    #Length distribution for proliferation events
    dist = rand(MultivariateNormal([0.0, 0.0], params.constants[5]))

    #Implement proliferation
    new_x = periodic(vars.P[1, idx] + dist[1], params.L)
    new_y = periodic(vars.P[2, idx] + dist[2], params.L)

    new_pos = [new_x, new_y]

    #Update position matrix
    vars.P = hcat(vars.P, new_pos)

    #Increase size of rates matrix
    vars.rates = [vars.rates zeros(3, 1)]

    #Update distance matrix
    vars.dist = [vars.dist zeros(vars.N_particles - 1, 1)]
    vars.dist = [vars.dist; zeros(1, vars.N_particles)]

    #Compute distance to other particles
    @inbounds for i in 1:vars.N_particles

        #Compute distance to other particles
        xi_x = periodic(vars.P[1, i] - new_pos[1], params.L)
        xi_y = periodic(vars.P[2, i] - new_pos[2], params.L)

        ξ_squared = xi_x^2 + xi_y^2

        # Fill distance matrix
        vars.dist[i, vars.N_particles] = sqrt(ξ_squared)
        vars.dist[vars.N_particles, i] = sqrt(ξ_squared)

        # Compute rates for new particle
        if i ≠ vars.N_particles

            # Update rates
            ρ_c = competition_kernel(ξ_squared, params) #X_n

            ρ_p = proliferation_kernel(ξ_squared, params) #Y_n

            vars.rates[2, vars.N_particles] += ρ_p  #P_n
            vars.rates[3, vars.N_particles] += ρ_c #D_n

            # Compute change in rates for other particles
            vars.rates[3, i] += ρ_c #D_n
            vars.rates[1, i] = params.constants[6] .+ vars.rates[3, i] .^ 2 #D_squared_n
            vars.rates[2, i] += ρ_p  #P_n

        end

    end

    vars.rates[1, vars.N_particles] = params.constants[6] + vars.rates[3, vars.N_particles]^2 #D_squared_n
    vars.rates[2, vars.N_particles] += params.constants[7]  #P_n

end

# Implement death event (removes a particle at a given position)
@views function remove_particle(idx, vars, params)

    old_pos = vars.P[:, idx]

    vars.N_particles -= 1

    #Update position matrix
    vars.P = hcat(vars.P[:, 1:idx-1], vars.P[:, idx+1:vars.N_particles+1])

    #Update distance matrix
    vars.dist = [vars.dist[:, 1:idx-1] vars.dist[:, idx+1:vars.N_particles+1]]
    vars.dist = [vars.dist[1:idx-1, :]; vars.dist[idx+1:vars.N_particles+1, :]]

    # Decrease size of rates matrix
    vars.rates = [vars.rates[:, 1:idx-1] vars.rates[:, idx+1:vars.N_particles+1]]

    # Compute change in rates for other particles
    @inbounds for i in 1:vars.N_particles

        #Compute distance to other particles
        xi_x = periodic(vars.P[1, i] - old_pos[1], params.L)
        xi_y = periodic(vars.P[2, i] - old_pos[2], params.L)

        ξ_squared = xi_x^2 + xi_y^2

        # Update rates
        ρ_c = competition_kernel(ξ_squared, params) #X_n

        ρ_p = proliferation_kernel(ξ_squared, params) #Y_n

        # Compute change in rates for other particles
        vars.rates[3, i] -= ρ_c #D_n
        vars.rates[1, i] = params.constants[6] .+ vars.rates[3, i] .^ 2 #D_squared_n
        vars.rates[2, i] -= ρ_p  #P_n

    end

end

# Compute waiting time and total rate for next event (Death, Proliferation). See Gillespie algorithm.
@views function compute_waiting_time(vars, params)

    params.total_rates[1] = sum(vars.rates[1, :]) #Total death rate
    params.total_rates[2] = sum(vars.rates[2, :]) #Total proliferation rate

    W = sum(params.total_rates)

    tau = -1 / W * log(rand()) #Waiting time

    return tau, W

end

# Implement proliferation event (using the previous add_particle function)
@views function proliferation_event(vars, params)

    #Select random particle
    cum_sum = cumsum(vars.rates[2, :])
    R = rand() * params.total_rates[2]
    idx = findfirst(cum_sum .> R)

    add_particle(idx, vars, params)
end

# Implement death event (using the previous remove_particle function)
@views function death_event(vars, params)

    #Select random particle
    cum_sum = cumsum(vars.rates[1, :])
    R = rand() * params.total_rates[1]
    idx = findfirst(cum_sum .> R)

    remove_particle(idx, vars, params)

end

# Choose event to happen (Death, Proliferation)
@views function choose_apply_event(reactions, orders, W, params, vars)

    U = rand() * W

    @inbounds for i in 1:2 #Three event types: Death, Proliferation, Movement

        if U < sum(params.total_rates[orders[1:i]])

            #Execute corresponding reaction
            reactions[i](vars, params)

            if i != 1

                aux_r = copy(reactions)
                aux_o = copy(orders)

                reactions[i-1] = aux_r[i]
                reactions[i] = aux_r[i-1]

                orders[i-1] = aux_o[i]
                orders[i] = aux_o[i-1]

            end

            break

        end

    end

    return reactions, orders

end

# Gillespie algorithm (returns the number of particles at given times)
@views function IBM(t_max, vars, params)

    #Initialization
    N_t = []
    t_res = []

    reactions = [death_event, proliferation_event]

    orders = [1, 2]

    t = 0

    initialize_vars(vars, params)

    count = 1

    # N_t[count] = vars.N_particles
    push!(N_t, vars.N_particles)
    push!(t_res, t)

    count2 = 1

    kernel_type = "LP"

    t_hp = 0.0

    while t < t_max

        if (kernel_type == "LP") & (t > vars.switching_times[count2])

            println("Switching to HP at t = ", t)

            # Change the facilitation kernel strength
            params.constants[3] = params.constants[3] * params.constants[9]

            # Recompute rates
            recompute_rates(vars, params)

            kernel_type = "HP"

            count2 += 1

        elseif (kernel_type == "HP") & (t_hp > params.constants[8])

            println("Switching to LP at t = ", t)

            # Change the facilitation kernel strength
            params.constants[3] = params.constants[3] / params.constants[9]

            # Recompute rates
            recompute_rates(vars, params)

            kernel_type = "LP"

            t_hp = 0.0

        else

            # Should not happen I think, but just in case
            if t > vars.switching_times[count2]

                println("As t = ", t, " > switching time ", vars.switching_times[count2], " adding one count")

                count2 += 1

            end
            
        end

        if (vars.N_particles == 0)

            println("Extinction at t = ", t, " with ", vars.N_particles, " particles")

            break

        end

        #Compute total rate
        tau, W = compute_waiting_time(vars, params)

        t += tau

        if kernel_type == "HP"

            t_hp += tau

        end

        # Choose event to happen
        reactions, orders = choose_apply_event(reactions, orders, W, params, vars)

        # Record values of interest
        push!(N_t, vars.N_particles)
        push!(t_res, t)

    end

    return N_t, t_res

end

# Gillespie algorithm (saves a txt file with the number of particles at given times and several txt files with the positions of the particles at given times)
@views function IBM_save_all(t_max, vars, params; folder="Results", save_at=1.0)

    #Initialization
    N_t = []
    t_res = []

    reactions = [death_event, proliferation_event]

    orders = [1, 2]

    t = 0.0
    t_save_at = 0

    initialize_vars(vars, params)

    count = 1

    count2 = 1

    kernel_type = "LP"

    t_hp = 0.0

    push!(N_t, vars.N_particles)
    push!(t_res, t)

    # Save positions of particles
    writedlm(folder * "/P_" * string(t) * ".txt", vars.P)

    while t < t_max

        if (kernel_type == "LP") & (t > vars.switching_times[count2])

            println("Switching to HP at t = ", t)

            # Change the facilitation kernel strength
            params.constants[3] = params.constants[3] * params.constants[9]

            # Recompute rates
            recompute_rates(vars, params)

            kernel_type = "HP"

            count2 += 1

        elseif (kernel_type == "HP") & (t_hp > params.constants[8])

            println("Switching to LP at t = ", t)

            # Change the facilitation kernel strength
            params.constants[3] = params.constants[3] / params.constants[9]

            # Recompute rates
            recompute_rates(vars, params)

            kernel_type = "LP"

            t_hp = 0.0

        else

            # Should not happen I think, but just in case
            if t > vars.switching_times[count2]

                println("As t = ", t, " > switching time ", vars.switching_times[count2], " adding one count")

                count2 += 1

            end
            
        end

        if (vars.N_particles == 0)

            println("Extinction at t = ", t, " with ", vars.N_particles, " particles")

            break

        end

        #Compute total rate
        tau, W = compute_waiting_time(vars, params)

        t += tau

        if kernel_type == "HP"

            t_hp += tau

        end

        #Choose event to happen
        reactions, orders = choose_apply_event(reactions, orders, W, params, vars)

        # Record values of interest
        if t - t_save_at > save_at

            # Save positions of particles
            writedlm(folder * "/P_" * string(t) * ".txt", vars.P)

            push!(N_t, vars.N_particles)
            push!(t_res, t)

            t_save_at += save_at

        end


    end

    # Save number of particles
    writedlm(folder * "/N_t.txt", N_t)
    writedlm(folder * "/t_res.txt", t_res)

end