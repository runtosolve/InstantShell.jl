using Symbolics, StaticArrays, HCubature


@variables x y a1 a2 a3 a4 a5 a6 a7 a8 a9 a10 a11 a12 a b E ν G t 

Dx = Differential(x)
Dy = Differential(y)

a_terms = Symbolics.variables(:a, 1:12)

w_terms = [1.0, x, y, x^2, x*y, y^2, x^3, x^2*y, x*y^2, y^3, x^3 * y, x * y^3]

# a_membrane_terms = []
# uv_membrane_terms = [1.0, x, y, xy]

# #triangular
# uv_membrane_terms = [1.0, x, y]

w = w_terms' * a_terms


C = Matrix{Num}(undef, (12, 12))

w1 = substitute(w, (Dict(x => 0.0, y => 0.0)))
C[1, :] = convert(Vector{Num}, [Symbolics.coeff(w1, a_terms[i]) for i in eachindex(a_terms)])

w1x = substitute(expand_derivatives(Dy(w)), (Dict(x => 0.0, y => 0.0)))
C[2, :] = convert(Vector{Num}, [Symbolics.coeff(w1x, a_terms[i]) for i in eachindex(a_terms)])

w1y = substitute(expand_derivatives(Dx(w)), (Dict(x => 0.0, y => 0.0)))
C[3, :] = -convert(Vector{Num}, [Symbolics.coeff(w1y, a_terms[i]) for i in eachindex(a_terms)])

w2 = substitute(w, (Dict(x => a, y => 0.0)))
C[4, :] = convert(Vector{Num}, [Symbolics.coeff(w2, a_terms[i]) for i in eachindex(a_terms)])

w2x = substitute(expand_derivatives(Dy(w)), (Dict(x => a, y => 0.0)))
C[5, :] = convert(Vector{Num}, [Symbolics.coeff(w2x, a_terms[i]) for i in eachindex(a_terms)])

w2y = substitute(expand_derivatives(Dx(w)), (Dict(x => a, y => 0.0)))
C[6, :] = -convert(Vector{Num},[Symbolics.coeff(w2y, a_terms[i]) for i in eachindex(a_terms)])

w3 = substitute(w, (Dict(x => a, y => b)))
C[7, :] = convert(Vector{Num},[Symbolics.coeff(w3, a_terms[i]) for i in eachindex(a_terms)])

w3x = substitute(expand_derivatives(Dy(w)), (Dict(x => a, y => b)))
C[8, :] = convert(Vector{Num},[Symbolics.coeff(w3x, a_terms[i]) for i in eachindex(a_terms)])

w3y = substitute(expand_derivatives(Dx(w)), (Dict(x => a, y => b)))
C[9, :] = -convert(Vector{Num},[Symbolics.coeff(w3y, a_terms[i]) for i in eachindex(a_terms)])

w4 = substitute(w, (Dict(x => 0.0, y => b)))
C[10, :] = convert(Vector{Num},[Symbolics.coeff(w4, a_terms[i]) for i in eachindex(a_terms)])

w4x = substitute(expand_derivatives(Dy(w)), (Dict(x => 0.0, y => b)))
C[11, :] = convert(Vector{Num}, [Symbolics.coeff(w4x, a_terms[i]) for i in eachindex(a_terms)])

w4y = substitute(expand_derivatives(Dx(w)), (Dict(x => 0.0, y => b)))
C[12, :] = -convert(Vector{Num}, [Symbolics.coeff(w4y, a_terms[i]) for i in eachindex(a_terms)])


P = [ w_terms'

      expand_derivatives(Dy.(w_terms))'

      -expand_derivatives(Dx.(w_terms))'

]


N = P * inv(C, laplace=false)



Q = Matrix{Num}(undef, (3, 12))

κx = -expand_derivatives(Dx(Dx(w)))
κy = -expand_derivatives(Dy(Dy(w)))
κxy = expand_derivatives(Dx(Dy(w)))
Q[1, :] = convert(Vector{Num}, [Symbolics.coeff(κx, a_terms[i]) for i in eachindex(a_terms)])
Q[2, :] = convert(Vector{Num}, [Symbolics.coeff(κy, a_terms[i]) for i in eachindex(a_terms)])
Q[3, :] = -2 * convert(Vector{Num}, [Symbolics.coeff(κxy, a_terms[i]) for i in eachindex(a_terms)])


D = t^3/12  * [E/(1-ν^2)      (ν*E)/(1-ν^2)        0
           
               (ν*E)/(1-ν^2)     E/(1-ν^2)         0

                  0           0                 G]


dK = transpose(Q * inv(C, laplace=false)) * D * (Q * inv(C, laplace=false))                  


dK = substitute(dK, (Dict(E => 1.0, ν => 0.3, t => 1.0, G => 0.5, a => 1.0, b => 1.0)))



f_expr = Base.remove_linenums!(build_function(dK[1, 1], [x, y], expression=Val{false}))

f = build_function(dK[1, 1], [x, y], expression=Val{true})

# f_expr = build_function(dU, [x, y])
f = eval(f_expr)
# f(SA[0.0, 0.15])


K = hcubature(f, SA[0.0, 0.0], SA[a, b]) 
