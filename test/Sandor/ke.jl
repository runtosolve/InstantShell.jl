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



