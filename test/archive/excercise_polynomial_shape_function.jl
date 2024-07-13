using Symbolics, StaticArrays, SymbolicUtils, HCubature, Optimization, OptimizationOptimJL, CairoMakie



@variables x y a1 a2 a3 a4 a5 a6 a7 a8 a9 a10 a11 a12 a b

Dx = Differential(x)
Dy = Differential(y)

a_terms = [a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11, a12]

# define shape function
w = a1 + a2 * x + a3 * y + a4 * x^2 + a5 * x * y + a6 * y^2 + a7 * x^3 + a8 * x^2 * y + a9 * x * y^2 + a10 * y^3 + a11 * x^3 * y + a12 * x * y^3

A = Matrix{Any}(undef, (12, 12))

w1 = substitute(w, (Dict(x => -a, y => -b)))
A[1, :] = [Symbolics.coeff(w1, a_terms[i]) for i in eachindex(a_terms)]

w1y = substitute(expand_derivatives(Dx(w)), (Dict(x => -a, y => -b)))
A[2, :] = [Symbolics.coeff(w1y, a_terms[i]) for i in eachindex(a_terms)]

w1x = substitute(expand_derivatives(Dy(w)), (Dict(x => -a, y => -b)))
A[3, :] = [Symbolics.coeff(w1x, a_terms[i]) for i in eachindex(a_terms)]

w2 = substitute(w, (Dict(x => -a, y => b)))
A[4, :] = [Symbolics.coeff(w2, a_terms[i]) for i in eachindex(a_terms)]

w2y = substitute(expand_derivatives(Dx(w)), (Dict(x => -a, y => b)))
A[5, :] = [Symbolics.coeff(w2y, a_terms[i]) for i in eachindex(a_terms)]

w2x = substitute(expand_derivatives(Dy(w)), (Dict(x => -a, y => b)))
A[6, :] = [Symbolics.coeff(w2x, a_terms[i]) for i in eachindex(a_terms)]

w3 = substitute(w, (Dict(x => a, y => b)))
A[7, :] = [Symbolics.coeff(w3, a_terms[i]) for i in eachindex(a_terms)]

w3y = substitute(expand_derivatives(Dx(w)), (Dict(x => a, y => b)))
A[8, :] = [Symbolics.coeff(w3y, a_terms[i]) for i in eachindex(a_terms)]

w3x = substitute(expand_derivatives(Dy(w)), (Dict(x => a, y => b)))
A[9, :] = [Symbolics.coeff(w3x, a_terms[i]) for i in eachindex(a_terms)]

w4 = substitute(w, (Dict(x => a, y => -b)))
A[10, :] = [Symbolics.coeff(w4, a_terms[i]) for i in eachindex(a_terms)]

w4y = substitute(expand_derivatives(Dx(w)), (Dict(x => a, y => -b)))
A[11, :] = [Symbolics.coeff(w4y, a_terms[i]) for i in eachindex(a_terms)]

w4x = substitute(expand_derivatives(Dy(w)), (Dict(x => a, y => -b)))
A[12, :] = [Symbolics.coeff(w4x, a_terms[i]) for i in eachindex(a_terms)]

#w = Aa

#ri = [wi, θxi, θyi]

Δ = [0.0, -0.01,-0.01, 0.0, -0.01, 0.01, 0.0, 0.01, 0.01, 0.0, 0.01, -0.01]

a_element = 1.0
b_element = 1.0

A_num = reshape([substitute(A[i,j], (Dict([a => a_element, b => b_element]))) for i=1:12 for j=1:12], (12, 12))'

a_sol = A_num \ Δ

w_part = substitute(w, (Dict([a1 => a_sol[1], a2 => a_sol[2], a3 => a_sol[3], a4 => a_sol[4], a5 => a_sol[5], a6 => a_sol[6], a7 => a_sol[7], a8 => a_sol[8], a9 => a_sol[9], a10 => a_sol[10], a11 => a_sol[11], a12 => a_sol[12]])))

x_range = LinRange(-a_element, a_element, 10)
y_range = LinRange(-b_element, b_element, 10)


w_shape = [substitute(w_part, (Dict([x=> x_range[i], y=> y_range[j]]))).val for i in eachindex(x_range), j in eachindex(y_range)]


using CairoMakie

f = surface(x_range, y_range, w_shape, axis=(type=Axis3,))
