export compute_lambda, generate_lambda_list

function compute_lambda(model::BilevelProblem; λ_choice::String, nb_runs::Int = 100)
    @assert λ_choice ∈ ["LL", "UL"] "λ_choice error: Select a choice for computing λ among UL or LL"

    # Get the model and initial points

    nx = model.dim[1]
    ny = model.dim[2]
    x0y0 = model.xy0
    
    xk = x0y0[1:nx]
    yk = x0y0[nx+1:nx+ny]
    F = model.F_func
    f = model.f_func
    G = model.G_func
    g = model.g_func

    t_UL = 0
    t_LL = 0

    # computes the mean of the times on 100 runs
    for i in 1:nb_runs
        # Begin procedure to compute time of UL
        t_UL_start = time()
        F(xk, yk)
        G(xk, yk)
        t_UL_i = time() - t_UL_start
        t_UL += t_UL_i

        # Begin procedure to compute time of LL
        t_LL_start = time()
        f(xk, yk)
        g(xk, yk)
        t_LL_i = time() - t_LL_start
        t_LL += t_LL_i
    end

    # Compute the ratio for λ
    λ = (λ_choice == UL) ? t_UL / t_LL : t_LL / t_UL
    return λ, t_UL, t_LL
end

function generate_lambda_list(prob_list; λ_choice::String, nb_runs::Int = 100)
    t_UL_list = zeros(length(prob_list))
    t_LL_list = similar(t_UL_list)
    λ_list = similar(t_UL_list)

    @inbounds for i in eachindex(prob_list)
        # Get the value of λ and other CPU times
        prob = prob_list[i]
        model = get_bilevel_problem(prob)
        λ, t_UL, t_LL = compute_lambda(model; λ_choice = λ_choice, nb_runs = nb_runs)

        t_UL_list[prob] = t_UL/nb_runs
        t_LL_list[prob] = t_LL/nb_runs
        λ_list[prob] =  t_UL / t_LL
    end

    return λ_list, t_UL_list, t_LL_list
end