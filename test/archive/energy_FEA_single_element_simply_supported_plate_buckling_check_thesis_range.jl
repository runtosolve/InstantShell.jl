using Symbolics, StaticArrays, SymbolicUtils, HCubature, Optimization, OptimizationOptimJL, CairoMakie




function calculate_strain_energy(ac, bc, μc, δc, tc, Ec)
    @variables x y δ a b μ
    # w = δ * sin(2π * x / a) * sin(2π * y / b)
    # w = δ * (1 - cos(2π * x / a)) * (1 - cos(2π * y / b))
    w = δ * sin(π * x/(2a)) * sin(π * y/(2b))

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

    U = hcubature(f, SA[-ac, -bc], SA[ac, bc]) 

    return U[1] * D/2

end

function calculate_external_work(ac, bc, μc, δc, tc, δoc, σxc, σyc, τxyc)
    #calculate external work
    @variables x y δ a b μ σx σy τxy δo
    w = δ * sin(π * x/(2a)) * sin(π * y/(2b))

    Dx = Differential(x)
    Dy = Differential(y)

    dV = σx * (Dx(w))^2 + σy * (Dy(w))^2 + 2 * τxy * (Dx(w)*Dy(w))

    dV = expand_derivatives(dV)

    dV = substitute(dV, (Dict(μ => μc, a => ac, b => bc, δ => δc, δo => δoc, σx => σxc, σy => σyc, τxy => τxyc)))
   
    f_expr = Base.remove_linenums!(build_function(dV, [x, y], expression=Val{false}))
    f = eval(f_expr)

    V = hcubature(f, SA[-ac, -bc], SA[ac, bc])

    return V[1] * tc/2

end


ac = 50.0
bc = 50.0
μc = 0.3
δoc = 1.0
Ec = 203000  #N/mm^2
tc = 1.2
σyc = 0.0 #N/mm^2
τxyc = 0.0


fcr = 4.0 * π^2 * Ec / (12 * (1-μc^2))*(tc/bc)^2

# fcr = (10.07 * D * π^2 / ac^2) / tc

σxc = range(0, -150.0, 50) #N/mm^2

δc = [1.0, 1.000000000001]

δΠ = Vector{Float64}(undef, length(σxc))

for i in eachindex(σxc)

    U = calculate_strain_energy.(ac, bc, μc, δc, tc, Ec)
    V = calculate_external_work.(ac, bc, μc, δc, tc, δoc, σxc[i], σyc, τxyc)

    δΠ[i] = (U[2]-U[1]) + (V[2]-V[1])

end

f = Figure()
ax = Axis(f[1, 1])
scatterlines!(ax, σxc, δΠ)
f



