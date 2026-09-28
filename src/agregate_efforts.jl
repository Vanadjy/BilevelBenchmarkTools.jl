export agregate_efforts!

function agregate_efforts!(N_all_hists, N_UL_all_hists, N_LL_all_hists, prob_list; λ_choice::String, nb_runs::Int = 100)
    λ_list, _, _ = generate_lambda_list(prob_list; λ_choice = λ_choice, nb_runs = nb_runs)

    @assert length(λ_list) == length(N_UL_all_hists) "Dimension mismatch error: length of λ_list inapropriate."
    @assert length(N_all_hists) == length(N_UL_all_hists) "Dimension mismatch error: length of N_all_hists inapropriate."
    @assert length(N_LL_all_hists) == length(N_UL_all_hists) "Dimension mismatch error: length of N_LL_all_hists inapropriate."
    @assert λ_choice ∈ ["UL", "LL"] "λ_choice error: Select a choice for computing λ among UL or LL"

    for prob in eachindex(prob_numbers)
        for algo in keys(N_UL_all_hists[prob])
            if λ_choice == "LL"
                @assert effort_choice ∈ ["LL", "Agregate"] "Invalid effort choice. Must be one of: LL or Agregate when λ_choice is LL."
                N_all_hists[prob][algo] .= (λ_list[prob] .* N_UL_all_hists[prob][algo]) .+ N_LL_all_hists[prob][algo]
            else
                @assert effort_choice ∈ ["UL", "Agregate"] "Invalid effort choice. Must be one of: UL or Agregate when λ_choice is UL."
                N_all_hists[prob][algo] .= N_UL_all_hists[prob][algo] .+ (N_LL_all_hists[prob][algo])./λ_list[prob]
            end
        end
    end
end

function agregate_efforts(N_UL_all_hists, N_LL_all_hists, prob_list; λ_choice::String, nb_runs::Int = 100)
    N_all_hists = copy(N_UL_all_hists)
    agregate_efforts!(N_all_hists, N_UL_all_hists, N_LL_all_hists, prob_list; λ_choice = λ_choice, nb_runs = nb_runs)
    return N_all_hists
end