# plans coupants classiques
using CPLEX
using JuMP
include("robust.jl")


function plans_coupants_classique(n, m; Gamma = 3, c = 1000)
    A = loadData(n, m)
    println("A = ", A)
    model = Model(CPLEX.Optimizer)

    @variable(model, x[1:n] >= 0)

    @constraint(model, [i in 1:m], sum(A[i,j] * x[j] for j in 1:n) <= c)

    @objective(model, Max, sum(x[i] for i in 1:n))

    optimize!(model)
    println("non robuste : ",objective_value(model))
    println("x non robuste : ", value.(x))

    condition = true
    while condition
        a_bar = copy(A)
        x_val = value.(x)
        println("a_bar = ", a_bar)
        for i in 1:m
            ax = [A[i,j] * x_val[j] for j in 1:n]
            ind_ax_dec = sortperm(ax, rev = true) # plus grand en premier
            for k in 1:Gamma
                a_bar[i, ind_ax_dec[k]] = 1.01 * a_bar[i, ind_ax_dec[k]] 
            end
        end
        println("a_bar = ", a_bar)
        a_barx = [sum(a_bar[i,j] * x_val[j] for j in 1:n) for i in 1:m]
        if maximum(a_barx) <= c + 1e-5
            condition = false
        else
            i_max = argmax(a_barx)
            @constraint(model, sum(a_bar[i_max, j] * x[j] for j in 1:n) <= c) 
            optimize!(model)
        end
    end

    return value.(x), objective_value(model)
end

function plans_coupants_projection(n, m; Gamma = 3, c = 1000)
    function get_t(x_real, d, i)
        return((c - sum(A[i,j] * x_real[j] for j in 1:n)) / sum(A[i,j] * d[j] for j in 1:n))
    end

    A = loadData(n, m)
    model = Model(CPLEX.Optimizer)

    @variable(model, x[1:n] >= 0)

    @constraint(model, [i in 1:m], sum(A[i,j] * x[j] for j in 1:n) <= c)

    d = [1 for i in 1:n]
    @objective(model, Max, sum(x[i] * d[i] for i in 1:n))

    
    x_real = [0 for i in 1:n]

    condition = true
    while condition
        best_t = Inf
        best_ind_a = 0
        a_bar = copy(A)
        for i in 1:m
            if  d' * A[i,:] > 1e-5
                t = get_t(x_real, d, i)
                ax_td = [A[i,j] * (x_real[j] + t * d[j]) for j in 1:n]
                ind_ax_dec = sortperm(ax_td, rev = true) # plus grand en premier
                for k in 1:Gamma
                    a_bar[i, ind_ax_dec[k]] = 1.01 * a_bar[i, ind_ax_dec[k]] 
                end
                nv_t = (c - sum(a_bar[i,j] * x_real[j] for j in 1:n)) / sum(a_bar[i,j] * d[j] for j in 1:n)
                if (sum(a_bar[i,j] * (x_real[j] + nv_t * d[j]) for j in 1:n) > c - 1e-5) && (nv_t < best_t)
                    best_t = nv_t
                    best_ind_a = i
                end
            end
        end
        if best_t == Inf
            println("on sort de la boucle")
            condition = false
        else
            #@objective(model, Max, sum(x[i] * d[i] for i in 1:n))
            @constraint(model, sum(a_bar[best_ind_a, j] * x[j] for j in 1:n) <= c) 
            optimize!(model)

            x_bord = value.(x)
            println("x_real = ", x_real)
            println("d = ", d)  
            x_real = x_real + best_t/10 * d
            println("LB = ", sum(x_real[i] for i in 1:n))
            println("UB = ", sum(x_bord[i] for i in 1:n))
            println("x_real = ", x_real)
            println("x_bord = ", x_bord)
            d = x_bord - x_real
            println("d = ", d)  
        end
    end

    return x_real, sum(x_real[i] for i in 1:n)
end
