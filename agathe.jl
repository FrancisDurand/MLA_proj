# plans coupants classiques
using CPLEX
using JuMP
include("robust.jl")


function plans_coupants_Q1(n, m; Gamma = 3; c = 1000)
    A = loadData(n, m)
    model = Model(CPLEX.Optimizer)

    @variable(model, x[1:n] >= 0)

    @objective(model, Min, sum(x[i] for i in 1:n))

    optimize!(model)

    condition = true
    while condition
        x_val = value.(x)

        ss_pb = Model(CPLEX.Optimizer)

        ax = [A[i,:]'*x_val for i in 1:m]
        a_etoile = A[argmax(ax),:]
        println("a_etoile = ", a_etoile)
        ind_ax_dec = sortperm(ax, rev = true) # plus grand en premier
        for i in 1:Gamma
            a_etoile[ind_ax_dec[i]] = 1.01 * a_etoile[ind_ax_dec[i]] 
        end
        println("a_etoile = ", a_etoile)
        @constraint(model, sum(a_etoile[i] * x[i] for i in 1:n) <= c)
    end
end

function plans_coupants_Q2(n, m; Gamma = 3; c = 1000)
    A = loadData(n, m)
    model = Model(CPLEX.Optimizer)

    @variable(model, x[1:n] >= 0)

    @objective(model, Min, sum(x[i] for i in 1:n))

    optimize!(model)

    condition = true
    while condition
        x_val = value.(x)
        ax = [A[i,:]'*x_val for i in 1:m]
        a_etoile = A[argmax(ax),:]
        println("a_etoile = ", a_etoile)
        ind_ax_dec = sortperm(ax, rev = true) # plus grand en premier
        for i in 1:Gamma
            a_etoile[ind_ax_dec[i]] = 1.01 * a_etoile[ind_ax_dec[i]] 
        end
        println("a_etoile = ", a_etoile)
        @constraint(model, sum(a_etoile[i] * x[i] for i in 1:n) <= c)
    end
end