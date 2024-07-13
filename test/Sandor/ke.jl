using Symbolics

@variables E ν G t x y a b

Duv = t * [E / (1 - ν^2)        (ν * E) / (1 - ν^2)     0
           
          (ν * E) / (1 - ν^2)   E / (1 - ν^2)           0

          0                     0                       G]


Dwθ = (t^3 / 12) * Duv


D = [   Duv             zeros(Num, (3, 3))


        zeros(Num, (3, 3))   Dwθ]




Luv = [Differential(x)          0
      
       0                        Differential(y)  


       Differential(y)          Differential(x)]


Lwθ = [-1 * Differential(x) * Differential(x)

        -1 * Differential(y) * Differential(y)

        -2 * Differential(x) * Differential(y)]



L = [   Luv                     zeros(Num, (3, 1))


        zeros(Num, (3, 2))      Lwθ     ]



# N_o2_x_1 = 1 - 3x / a + (2x^2)/a^2
# N_o2_x_2 = -x / a + 2x^2 / a^2
# N_o1_y_1 = 1 - y / b
# N_o1_y_2 = y / b
# N_o3_x_1 = 1 - 3x^2 / a^2 + 2x^3/a^3
# N_o3_x_2 = x - 2x^2 / a + x^3 / a^2
# N_o3_y_1 = 1 - 3y^2 / b^2 + 2y^3 / b^3
# N_o3_y_2 = x - 2y^2 / b + y^3 / b^2


Nx1 = 1 - x/a 
Nx2 = x/a
Ny1 = 1 - y/b
Ny2 = y/b 

Nu = [Nx1 * Ny1       Nx1 * Ny2     Nx2 * Ny1     Nx2 * Ny2]
Nv = [Nx1 * Ny1       Nx1 * Ny2     Nx2 * Ny1     Nx2 * Ny2]


Nx1 = 1 - 3x^2 / a^2 + 2x^3 / a^3
Ny1 = 1 - 3y^2 / b^2 + 2y^3 / b^3
Nx2 = x - 2x^2 / a + x^3 / a^2
Ny2 = y - 2y^2 / b + y^3 / b^2


 #       w11              θx11          θy11            w12
Nwθ = [Nx1 * Ny1        Nx1 * Ny2       Ny1 * Nx2       Nx2 * Ny1   ] 



Nv = [N_o3_x_1 * N_o1_y_1       N_o3_x_1 * N_o1_y_2     N_o3_x_2 * N_o1_y_1     N_o3_x_2 * N_o1_y_2]

Nw = [N_o3_x_1 * N_o3_y_1       N_o3_x_1 * N_o3_y_2     N_o3_x_2 * N_o3_y_1     N_o3_x_2 * N_o3_y_2]

Nθx = []