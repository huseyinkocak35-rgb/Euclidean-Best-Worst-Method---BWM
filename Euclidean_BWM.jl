using JuMP
using Ipopt

function BWM(BO_input::Vector{T}, OW_input::Vector{T}) where T <: Real
    BO = Float64.(BO_input)
    OW = Float64.(OW_input)
    
    n = length(BO)
    if n != length(OW)
        error("Sizes of BO and OW vectors must be equal and match the number of criteria.")
    end
    println("Number of criteria: ", n)
    
    abest = minimum(BO)
    aworst = maximum(BO)
    abest2 = maximum(OW)
    aworst2 = minimum(OW)
    
    if abest != 1.0 error("Best-to-best comparison should be 1.") end
    if aworst2 != 1.0 error("Worst-to-worst comparison should be 1.") end
    if aworst != abest2 error("Value of the best-worst comparison should be same in BO and OW.") end
    
    # Identify indices of best and worst criteria
    bob = findall(x -> x == abest, BO)
    bow = findall(x -> x == aworst, BO)
    owb = findall(x -> x == abest2, OW)
    oww = findall(x -> x == aworst2, OW)
    
    if length(bob) != length(owb) error("Number of best criteria should be same in BO and OW.") end
    if length(bow) != length(oww) error("Number of worst criteria should be same in BO and OW.") end
    if bob != owb error("Best criteria are not same in BO and OW.") end
    if bow != oww error("Worst criteria are not same in BO and OW.") end
    
    # Display single or multiple best/worst criteria
    if length(bob) == 1
        println("The best criterion: C", bob[1])
    else
        println("The best criteria: ", join(["C$j" for j in bob], ", "))
    end

    if length(bow) == 1
        println("The worst criterion: C", bow[1])
    else
        println("The worst criteria: ", join(["C$j" for j in bow], ", "))
    end
    println("Value of the best-worst comparison: ", aworst)
    println("-"^50)

    # Reference the last identified index (consistent with original formulation)
    best_idx = bob[end]
    worst_idx = bow[end]

    # ==========================================
    # 1. EUCLIDEAN BWM COMPUTATION
    # ==========================================
    model_euc = Model(Ipopt.Optimizer)
    set_silent(model_euc)
    @variable(model_euc, w[1:n] >= 1e-6)
    @constraint(model_euc, sum(w) == 1)
    
    @NLobjective(model_euc, Min, 
        sum((w[best_idx]/w[i] - BO[i])^2 for i in 1:n if i != best_idx) + 
        sum((w[i]/w[worst_idx] - OW[i])^2 for i in 1:n if i != best_idx && i != worst_idx)
    )
    optimize!(model_euc)
    w_euc = value.(w)
    xi_euc_val = sqrt(objective_value(model_euc))
    
    # Euclidean Consistency Index (CI) Model
    model_euc_ci = Model(Ipopt.Optimizer)
    set_silent(model_euc_ci)
    @variable(model_euc_ci, w1[1:n] >= 1e-6)
    @constraint(model_euc_ci, sum(w1) == 1)
    
    @NLobjective(model_euc_ci, Min, 
        sum((w1[best_idx]/w1[i] - aworst)^2 for i in 1:n if i != best_idx) + 
        sum((w1[i]/w1[worst_idx] - aworst)^2 for i in 1:n if i != best_idx && i != worst_idx)
    )
    optimize!(model_euc_ci)
    ci_euc_val = sqrt(objective_value(model_euc_ci))
    
    # Total Deviation (TD) for Euclidean BWM
    sumTD11 = sum((BO[i] - (w_euc[best_idx]/w_euc[i]))^2 for i in 1:n)
    sumTD12 = sum((OW[i] - (w_euc[i]/w_euc[worst_idx]))^2 for i in 1:n)
    td_euc = (sumTD11 + sumTD12) / (2 * n)
    
    println("The Euclidean BWM results are as follows:")
    for i in 1:n
        println("\t w", i, " = ", round(w_euc[i], digits=5))
    end
    println("\t xi = ", round(xi_euc_val, digits=5))
    println("\t consistency index = ", round(ci_euc_val, digits=5))
    println("\t the consistency ratio = xi/consistency index = ", round(xi_euc_val / ci_euc_val, digits=5))
    println("\t Total Deviation = ", round(td_euc, digits=5))
    println("-"^50)

    # ==========================================
    # 2. LINEAR CHEBYSHEV BWM COMPUTATION
    # ==========================================
    model_cheb = Model(Ipopt.Optimizer)
    set_silent(model_cheb)
    @variable(model_cheb, w3[1:n] >= 0)
    @variable(model_cheb, x >= 0)
    @constraint(model_cheb, sum(w3) == 1)
    
    for i in 1:n
        if i != best_idx
            @constraint(model_cheb, w3[best_idx] - BO[i] * w3[i] <= x)
            @constraint(model_cheb, w3[best_idx] - BO[i] * w3[i] >= -x)
        end
        if i != best_idx && i != worst_idx
            @constraint(model_cheb, w3[i] - OW[i] * w3[worst_idx] <= x)
            @constraint(model_cheb, w3[i] - OW[i] * w3[worst_idx] >= -x)
        end
    end
    @objective(model_cheb, Min, x)
    optimize!(model_cheb)
    w_cheb = value.(w3)
    xi_cheb_val = value(x)
    
    # Linear Chebyshev Consistency Index (CI) Model
    model_cheb_ci = Model(Ipopt.Optimizer)
    set_silent(model_cheb_ci)
    @variable(model_cheb_ci, w4[1:n] >= 0)
    @variable(model_cheb_ci, y >= 0)
    @constraint(model_cheb_ci, sum(w4) == 1)
    
    for i in 1:n
        if i != best_idx
            @constraint(model_cheb_ci, w4[best_idx] - aworst * w4[i] <= y)
            @constraint(model_cheb_ci, w4[best_idx] - aworst * w4[i] >= -y)
        end
        if i != best_idx && i != worst_idx
            @constraint(model_cheb_ci, w4[i] - aworst * w4[worst_idx] <= y)
            @constraint(model_cheb_ci, w4[i] - aworst * w4[worst_idx] >= -y)
        end
    end
    @objective(model_cheb_ci, Min, y)
    optimize!(model_cheb_ci)
    max_xi_cheb_val = value(y)
    
    # Total Deviation (TD) for Chebyshev BWM
    sumTD21 = sum((BO[i] - (w_cheb[best_idx]/w_cheb[i]))^2 for i in 1:n)
    sumTD22 = sum((OW[i] - (w_cheb[i]/w_cheb[worst_idx]))^2 for i in 1:n)
    td_cheb = (sumTD21 + sumTD22) / (2 * n)
    
    println("The linear Chebyshev BWM results are as follows:")
    for i in 1:n
        println("\t w", i, " = ", round(w_cheb[i], digits=5))
    end
    println("\t xi = ", round(xi_cheb_val, digits=5))
    println("\t max xi = ", round(max_xi_cheb_val, digits=5))
    println("\t the actual consistency ratio = xi/max xi = ", round(xi_cheb_val / max_xi_cheb_val, digits=5))
    println("\t Total Deviation = ", round(td_cheb, digits=5))
end
