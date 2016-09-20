using JuMP

x0 = rand(9)

function foo(x...)
    sum(collect(x).^2)
end

function jump_energy(x0)
    model = Model()
    @defVar(model, x[i=1:length(x0)], start=x0[i])
    registerNLFunction(:foo, length(x0), foo, autodiff=true)
    @setNLObjective(model, Min, foo(x[1], x[2], x[3], x[4], x[5], x[6], x[7], x[8], x[9]))
    solve(model)
    @show getValue(x)
end

jump_energy(x0)