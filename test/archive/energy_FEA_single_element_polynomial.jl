using Symbolics, StaticArrays, SymbolicUtils, HCubature




# @variables x y μ a1 a2 a3 a4 a5 a6 a7 a8 a9 a10 a11 a12 δo a b
# @variables x y μ δ δo a b
# Dx = Differential(x)
# Dy = Differential(y)
# Dxx = Differential(x) * Differential(x)
# Dyy = Differential(y) * Differential(y)
# Dxy = Differential(x) * Differential(y)

#define shape function
# w = a1 + a2 * x + a3 * y + a4 * x^2 + a5 * x * y + a6 * y^2 + a7 * x^3 + a8 * x^2 * y + a9 * x * y^2 + a10 * y^3 + a11 * x^3 * y + a12 * x * y^3
