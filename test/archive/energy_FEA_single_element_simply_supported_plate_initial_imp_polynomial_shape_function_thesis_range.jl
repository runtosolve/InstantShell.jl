using Symbolics, StaticArrays, SymbolicUtils, HCubature, Optimization, OptimizationOptimJL, CairoMakie, NonlinearSolve, RuntimeGeneratedFunctions


include("test/A.jl")

function calculate_strain_energy(ac, bc, μc, δc, tc, Ec)
    @variables x y a b μ a1 a2 a3 a4 a5 a6 a7 a8 a9 a10 a11 a12
    w = a1 + a2 * x + a3 * y + a4 * x^2 + a5 * x * y + a6 * y^2 + a7 * x^3 + a8 * x^2 * y + a9 * x * y^2 + a10 * y^3 + a11 * x^3 * y + a12 * x * y^3

    Dxx = Differential(x) * Differential(x)
    Dyy = Differential(y) * Differential(y)
    Dxy = Differential(x) * Differential(y)

    #calculate strain energy
    dU = Dxx(w)^2 + Dyy(w)^2 + 2 * μ * Dxx(w) * Dyy(w) + 2 * (1 - μ) * Dxy(w)^2

    dU = expand_derivatives(dU)

    Ac = A(ac, bc)

    ac = Ac \ δ

    dU = substitute(dU, (Dict(μ => μc, a => ac, b => bc, a1 => ac[1], a2 => ac[2], a3 => ac[3], a4 => ac[4], a5 => ac[5], a6 => ac[6], a7 => ac[7], a8 => ac[8], a9 => ac[9], a10 => ac[10], a11 => ac[11], a12 => ac[12])))

    f_expr = Base.remove_linenums!(build_function(dU, [x, y], expression=Val{false}))

    # f_expr = build_function(dU, [x, y])
    f = eval(f_expr)
    # f(SA[0.0, 0.15])

    D = Ec * tc^3 / (12 * (1 - μc^2))

    U = hcubature(f, SA[-ac, bc], SA[ac, bc]) 

    return U[1] * D/2

end

function calculate_external_work(ac, bc, μc, δc, tc, δoc, σxc, σyc, τxyc)
    #calculate external work
    @variables x y a b μ σx σy τxy δo a1 a2 a3 a4 a5 a6 a7 a8 a9 a10 a11 a12
    w = a1 + a2 * x + a3 * y + a4 * x^2 + a5 * x * y + a6 * y^2 + a7 * x^3 + a8 * x^2 * y + a9 * x * y^2 + a10 * y^3 + a11 * x^3 * y + a12 * x * y^3
    wo = δo * sin(π * x/(2a)) * sin(π * y/(2b))

    Dx = Differential(x)
    Dy = Differential(y)

    dV = σx * (Dx(w) + Dx(wo))^2 + σy * (Dy(w) + Dy(wo))^2 + 2 * τxy * (Dx(w)*Dy(w) + Dx(wo)*Dy(wo))

    dV = expand_derivatives(dV)



    dV = substitute(dV, (Dict(μ => μc, a => ac, b => bc, a1 => a1c, a2 => a2c, a3 => a3c, a4 => a4c, a5 => a5c, a6 => ac6, a7 => a7c, a8 => a8c, a9 => a9c, a10 => a10c, a11 => a11c, a12 => a12c, δo => δoc, σx => σxc, σy => σyc, τxy => τxyc)))
   
    f_expr = Base.remove_linenums!(build_function(dV, [x, y], expression=Val{false}))
    f = eval(f_expr)

    V = hcubature(f, SA[-ac, -bc], SA[ac, bc])

    return V[1] * tc/2

end

function potential_energy(u, p)
    
    (ac, bc, μc, tc, δoc, σxc, σyc, τxyc) = p

    U = calculate_strain_energy(ac, bc, μc, u[1], u[2], u[3], u[4], u[5], u[6], u[7], u[8], u[9], u[10], u[11], u[12], tc, Ec)

    V = calculate_external_work(ac, bc, μc, u[1], u[2], u[3], u[4], u[5], u[6], u[7], u[8], u[9], u[10], u[11], u[12], tc, δoc, σxc, σyc, τxyc)

    Π = U + V

    return Π

end


ac = 50.0
bc = 50.0
μc = 0.3
δoc = 2.0
Ec = 203000  #N/mm^2
tc = 1.2
σxc = -100.0 #N/mm^2
σyc = 0.0 #N/mm^2
τxyc = 0.0


p = (ac, bc, μc, tc, δoc, σxc, σyc, τxyc)
u0 = zeros(Float64, 12)
# prob_int = IntervalNonlinearProblem(potential_energy, uspan, p)

prob = OptimizationProblem(potential_energy, u0, p)

# Import a solver package and solve the optimization problem
using OptimizationOptimJL
sol = solve(prob, NelderMead())

# OptimizationProblem(rosenbrock, u0, p)
# sol = solve(prob_int)
# u[i] = sol.u



# σxc = range(0.0, -100.0, 10) #N/mm^2
# u = Vector{Float64}(undef, length(σxc))

# for i in eachindex(σxc)
#     p = (ac, bc, μc, tc, δoc, σxc[i], σyc, τxyc)
#     uspan = (0.000001, 20.0)
#     prob_int = IntervalNonlinearProblem(potential_energy, uspan, p)
#     sol = solve(prob_int)
#     u[i] = sol.u
# end




f = Figure()
ax = Axis(f[1, 1])
scatterlines!(ax, u .+ δoc, -σxc)
f



