using Symbolics, StaticArrays, SymbolicUtils, HCubature, Optimization, OptimizationOptimJL, CairoMakie


function calculate_strain_energy(ac, bc, μc, δc, tc, Ec)
    @variables x y δ a b μ
    # w = δ * sin(2π * x / a) * sin(2π * y / b)
    w = δ * (1 - cos(2π * x / a)) * (1 - cos(2π * y / b))

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

U_Chajes(D, A, a) = (16 * D * π^4 * A^2) / a^2

ac = 100.0
bc = 100.0
μc = 0.3
δoc = 1.0
Ec = 203000  #N/mm^2
tc = 1.2
σyc = 0.0 #N/mm^2
τxyc = 0.0

D = Ec * tc^3 / (12 * (1 - μc^2))

δc = 1.0
U = calculate_strain_energy(ac, bc, μc, δc, tc, Ec)

U_theory = U_Chajes(D, δc, ac)

U / U_theory

