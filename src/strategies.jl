export referee_challenge, EndPoint_Referee, Complete_Referee!, Reverse_Referee!

function referee_challenge(k::Int, model::BilevelProblem, xHists, yHists, fHists, prob::Int, algo::String, referees, referee_options; tol_ref::R = 1e-3) where {R <: Float64}
    referee_flag = false
    x_star = xHists[prob][algo][:, k]
    f_star = fHists[prob][algo][k]

    # Reoptimization using NOMAD - our referee
    new_f = Inf
    f = model.f_func
    G = model.G_func
    x0y0 = model.xy0
    nx = model.dim[1]
    ny = model.dim[2]

    for ref_index in eachindex(referees)

        referee = referees[ref_index]
        y_new, new_f, neval_lower = referee(model, x_star, x0y0[nx+1:nx+ny]; referee_options[ref_index])

        if (all(G(x_star, y_new) .<= 0.0)) && (new_f < f_star - tol_ref)  # The BEST referee found a strictly better solution than the algo AND is feasible w.r.t upper-level constraints
            @info "A referee found a better final solution than $algo on problem $prob at iterate $k"
            referee_flag = true
            break # If one referee found a better solution, we stop and consider that the algo is invalidated at this iterate
        end
    end
    return referee_flag
end

function EndPoint_Referee(F_all_hists_adjusted, N_all_hists_adjusted, algo_names, prob_numbers, x_all_hists, y_all_hists, f_all_hists, referees::Vector{Union{String, Int}}, referee_options; tol_ref::Float64 = 1e-3)
    @assert length(referees) == length(referee_options) "Mismatch error: Each Referee should be assigned specific options"

    for prob in eachindex(prob_numbers)
        model = get_bilevel_problem(prob_numbers[prob])
        for a in eachindex(algo_names)
            algo = algo_names[a]
            if x_all_hists[prob][algo][end] !== x_all_hists[prob][algo][1] # If the algorithm did not moved from the starting point, ignore it
                Random.seed!(seed)
                k = length(f_all_hists[prob][algo])
                flag = referee_challenge(k, model, x_all_hists, y_all_hists, f_all_hists, prob, algo, referees, referee_options; tol_ref = tol_ref)
                if flag #referee found a better LL solution than algo
                    fill!(F_all_hists_adjusted[prob][algo][k], NaN)
                    fill!(N_all_hists_adjusted[prob][algo][k], NaN)
                end
            end
        end
    end
    orphans = Dict{String, Vector{Int}}(key => Int[] for key in algo_names)
    # After the process, remove all Inf entries from the historics
    for prob in eachindex(prob_numbers)
        for algo in algo_names
            filter!(!isnan, F_all_hists_adjusted[prob][algo])
            filter!(!isnan, N_all_hists_adjusted[prob][algo])
            if length(F_all_hists_adjusted[prob][algo]) == 0 # If the referee invalidated all the historic, count it as an orphaned run
                push!(orphans[algo], prob_numbers[prob])
            end
        end
    end

    return F_all_hists_adjusted, N_all_hists_adjusted, orphans
end

function Complete_Referee!(F_all_hists_adjusted, N_all_hists_adjusted, algo_names, prob_numbers, x_all_hists, y_all_hists, f_all_hists, referees::Vector{Union{String, Int}}, referee_options; tol_ref::Float64 = 1e-3)
    @assert length(referees) == length(referee_options) "Mismatch error: Each Referee should be assigned specific options"
    
    for prob in eachindex(prob_numbers)
        model = get_bilevel_problem(prob_numbers[prob])
        for algo in algo_names
            if x_all_hists[prob][algo][end] !== x_all_hists[prob][algo][1] # If the algorithm did not moved from the starting point, ignore it
                k = length(f_all_hists[prob][algo])
                while k >= 1
                    Random.seed!(seed)
                    flag = referee_challenge(k, model, x_all_hists, y_all_hists, f_all_hists, prob, algo, referees, referee_options; tol_ref = tol_ref)
                    if flag #referee found at least once a better LL solution than algo
                        F_all_hists_adjusted[prob][algo][k] = NaN # Set at Inf the corresponding value in the upper objective historic
                        N_all_hists_adjusted[prob][algo][k] = NaN # Set at Inf the corresponding value in the lower objective historic
                    end
                    k -= 1
                end
            end
        end
    end
    orphans = Dict{String, Vector{Int}}(key => Int[] for key in algo_names)
    # After the process, remove all Inf entries from the historics
    for prob in eachindex(prob_numbers)
        for algo in algo_names
            filter!(!isnan, F_all_hists_adjusted[prob][algo])
            filter!(!isnan, N_all_hists_adjusted[prob][algo])
            if length(F_all_hists_adjusted[prob][algo]) == 0 # If the referee invalidated all the historic, count it as an orphaned run
                push!(orphans[algo], prob_numbers[prob])
            end
        end
    end
    return F_all_hists_adjusted, N_all_hists_adjusted, orphans
end

function Reverse_Referee!(F_all_hists_adjusted, N_all_hists_adjusted, algo_names, prob_numbers, x_all_hists, y_all_hists, f_all_hists, referees::Vector{Union{String, Int}}, referee_options; tol_ref::Float64 = 1e-3)
    @assert length(referees) == length(referee_options) "Mismatch error: Each Referee should be assigned specific options"
    for prob in eachindex(prob_numbers)
        model = get_bilevel_problem(prob_numbers[prob])
        for algo in algo_names
            if x_all_hists[prob][algo][end] !== x_all_hists[prob][algo][1] # If the algorithm did not moved from the starting point, ignore it
                flag = true
                k = length(f_all_hists[prob][algo])
                while flag && k >= 1 # Once we found an admissible point in the historic, we stop
                    Random.seed!(seed)
                    flag = referee_challenge(k, model, x_all_hists, y_all_hists, f_all_hists, prob, algo, referees, referee_options; tol_ref = tol_ref)
                    if flag #referee found at least once a better LL solution than algo
                        F_all_hists_adjusted[prob][algo][k] = NaN # Set at Inf the corresponding value in the upper objective historic 
                        N_all_hists_adjusted[prob][algo][k] = NaN # Set at Inf the corresponding value in the lower objective historic
                    end
                    k -= 1
                end
            end
        end
    end
    orphans = Dict{String, Vector{Int}}(key => Int[] for key in algo_names)
    # After the process, remove all Inf entries from the historics
    for prob in eachindex(prob_numbers)
        for algo in algo_names
            filter!(!isnan, F_all_hists_adjusted[prob][algo])
            filter!(!isnan, N_all_hists_adjusted[prob][algo])
            if length(F_all_hists_adjusted[prob][algo]) == 0 # If the referee invalidated all the historic, count it as an orphaned run
                push!(orphans[algo], prob_numbers[prob])
            end
        end
    end
    return F_all_hists_adjusted, N_all_hists_adjusted, orphans
end