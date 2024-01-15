using Symbolics, StaticArrays, SymbolicUtils, HCubature, Optimization, OptimizationOptimJL, CairoMakie



function calculate_external_work(ac, bc, μc, δc, tc, δoc, σxc, σyc, τxyc)
    #calculate external work
    @variables x y δ a b μ σx σy τxy δo
    w = δ * (1 - cos(2π * x / a)) * (1 - cos(2π * y / b))

    Dx = Differential(x)
    Dy = Differential(y)

    dV = σx * (Dx(w))^2 + σy * (Dy(w))^2 + 2 * τxy * (Dx(w)*Dy(w))

    dV = expand_derivatives(dV)

    dV = substitute(dV, (Dict(μ => μc, a => ac, b => bc, δ => δc, δo => δoc, σx => σxc, σy => σyc, τxy => τxyc)))
   
    f_expr = Base.remove_linenums!(build_function(dV, [x, y], expression=Val{false}))
    f = eval(f_expr)

    V = hcubature(f, SA[0.0, 0.0], SA[ac, bc])

    return V[1] * tc/2

end


  



V_Chajes(Nx, A) = -(3 * Nx * π^2 * A^2) / 2

ac = 100.0
bc = 100.0
μc = 0.3
δoc = 1.0
Ec = 203000  #N/mm^2
tc = 1.2
σxc = -100.0 #N/mm^2
σyc = 0.0 #N/mm^2
τxyc = 0.0


δc = 1.0
V = calculate_external_work(ac, bc, μc, δc, tc, δoc, σxc, σyc, τxyc)

V_theory = V_Chajes(-σxc * tc, δc)



