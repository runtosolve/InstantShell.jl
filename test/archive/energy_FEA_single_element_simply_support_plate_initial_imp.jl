using Symbolics, StaticArrays, SymbolicUtils, HCubature, Optimization, OptimizationOptimJL, CairoMakie, NonlinearSolve




function calculate_strain_energy(ac, bc, μc, δc, tc, Ec)
    @variables x y δ a b μ
    w = δ * sin(π * x/(a)) * sin(π * y/(b))

    Dxx = Differential(x) * Differential(x)
    Dyy = Differential(y) * Differential(y)
    Dxy = Differential(x) * Differential(y)

    #calculate strain energy
    dU = Dxx(w)^2 + Dyy(w)^2 + 2 * μ * Dxx(w) * Dyy(w) + 2 * (1 - μ) * Dxy(w)^2

    dU = expand_derivatives(dU)

    dU = substitute(dU, (Dict(μ => μc, a => ac, b => bc, δ => δc)))

    f_expr = Base.remove_linenums!(build_function(dU, [x, y], expression=Val{false}))

    # f_expr = build_function(dU, [x, y])
    f = eval(f_expr)
    # f(SA[0.0, 0.15])

    D = Ec * tc^3 / (12 * (1 - μc^2))

    U = hcubature(f, SA[0.0, 0.0], SA[ac, bc]) 

    return U[1] * D/2

end

function calculate_external_work(ac, bc, μc, δc, tc, δoc, σxc, σyc, τxyc)
    #calculate external work
    @variables x y δ a b μ σx σy τxy δo
    w = δ * sin(π * x/(a)) * sin(π * y/(b))
    wo = δo * sin(π * x/(a)) * sin(π * y/(b))

    Dx = Differential(x)
    Dy = Differential(y)

    dV = σx * (Dx(w) + Dx(wo))^2 + σy * (Dy(w) + Dy(wo))^2 + 2 * τxy * (Dx(w)*Dy(w) + Dx(wo)*Dy(wo))

    dV = expand_derivatives(dV)

    dV = substitute(dV, (Dict(μ => μc, a => ac, b => bc, δ => δc, δo => δoc, σx => σxc, σy => σyc, τxy => τxyc)))
   
    f_expr = Base.remove_linenums!(build_function(dV, [x, y], expression=Val{false}))
    f = eval(f_expr)

    V = hcubature(f, SA[0.0, 0.0], SA[ac, bc])

    return V[1] * tc/2

end

function potential_energy(u, p)
    
    (ac, bc, μc, tc, δoc, σxc, σyc, τxyc) = p

    U = calculate_strain_energy(ac, bc, μc, u[1], tc, Ec)

    V = calculate_external_work(ac, bc, μc, u[1], tc, δoc, σxc, σyc, τxyc)

    Π = U + V

    return Π

end


ac = 100.0
bc = 100.0
μc = 0.3
δoc = 2.0
Ec = 203000  #N/mm^2
tc = 1.2
σxc = -100.0 #N/mm^2
σyc = 0.0 #N/mm^2
τxyc = 0.0


σxc = range(0.0, -60.0, 10) #N/mm^2
u = Vector{Float64}(undef, length(σxc))

for i in eachindex(σxc)
    p = (ac, bc, μc, tc, δoc, σxc[i], σyc, τxyc)
    uspan = (0.000001, 20.0)
    prob_int = IntervalNonlinearProblem(potential_energy, uspan, p)
    sol = solve(prob_int)
    u[i] = sol.u
end




f = Figure()
ax = Axis(f[1, 1])
scatterlines!(ax, u .+ δoc, -σxc)
f



