using Symbolics



@variables x y a1 a2 a3 a4 a5 a6 a7 a8 a9 a10 a11 a12 a b

Dx = Differential(x)
Dy = Differential(y)

a_terms = [a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11, a12]

# define shape function
w = a1 + a2 * x + a3 * y + a4 * x^2 + a5 * x * y + a6 * y^2 + a7 * x^3 + a8 * x^2 * y + a9 * x * y^2 + a10 * y^3 + a11 * x^3 * y + a12 * x * y^3

C = Matrix{Any}(undef, (12, 12))

w1 = substitute(w, (Dict(x => 0.0, y => 0.0)))
C[1, :] = [Symbolics.coeff(w1, a_terms[i]) for i in eachindex(a_terms)]

w1x = substitute(expand_derivatives(Dy(w)), (Dict(x => 0.0, y => 0.0)))
C[2, :] = [Symbolics.coeff(w1x, a_terms[i]) for i in eachindex(a_terms)]

w1y = substitute(expand_derivatives(Dx(w)), (Dict(x => 0.0, y => 0.0)))
C[3, :] = [Symbolics.coeff(w1y, a_terms[i]) for i in eachindex(a_terms)]

w2 = substitute(w, (Dict(x => a, y => 0.0)))
C[4, :] = [Symbolics.coeff(w2, a_terms[i]) for i in eachindex(a_terms)]

w2x = substitute(expand_derivatives(Dy(w)), (Dict(x => a, y => 0.0)))
C[5, :] = [Symbolics.coeff(w2x, a_terms[i]) for i in eachindex(a_terms)]

w2y = substitute(expand_derivatives(Dx(w)), (Dict(x => a, y => 0.0)))
C[6, :] = [Symbolics.coeff(w2y, a_terms[i]) for i in eachindex(a_terms)]

w3 = substitute(w, (Dict(x => a, y => b)))
C[7, :] = [Symbolics.coeff(w3, a_terms[i]) for i in eachindex(a_terms)]

w3x = substitute(expand_derivatives(Dy(w)), (Dict(x => a, y => b)))
C[8, :] = [Symbolics.coeff(w3x, a_terms[i]) for i in eachindex(a_terms)]

w3y = substitute(expand_derivatives(Dx(w)), (Dict(x => a, y => b)))
C[9, :] = [Symbolics.coeff(w3y, a_terms[i]) for i in eachindex(a_terms)]

w4 = substitute(w, (Dict(x => a, y => 0.0)))
C[10, :] = [Symbolics.coeff(w4, a_terms[i]) for i in eachindex(a_terms)]

w4x = substitute(expand_derivatives(Dy(w)), (Dict(x => a, y => 0.0)))
C[11, :] = [Symbolics.coeff(w4x, a_terms[i]) for i in eachindex(a_terms)]

w4y = substitute(expand_derivatives(Dx(w)), (Dict(x => a, y => 0.0)))
C[12, :] = [Symbolics.coeff(w4y, a_terms[i]) for i in eachindex(a_terms)]



Symbolics.get_variables(w)

expand_derivatives(Dy(w))

#w = Aa


