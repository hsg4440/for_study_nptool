import numpy as np

Nx = 12288
Ny = 24320

Sx = 700 * 110/3
Sy = 20*960

Stot = (97220)**2

ratio = 0.0273


S_alu = Nx*Sy + Ny*Sx

expected = ratio*Stot
#print(S_alu/expected)
#print(S_alu/Stot)


r = (2*(20*960+40*660))/(1040*700)
print(r)
