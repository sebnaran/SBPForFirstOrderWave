from scipy.sparse import csr_matrix,lil_matrix
from scipy.sparse.linalg import norm, inv
from Functions import *
import math

N            = 1000 
a            = -1
b            = 1
eps          = 1
mu           = 1
CFL          = 0.5
T            = 1
xp,xm,dx     = StaggeredSpatialMesh(N,a,b)
dt           = CFL*dx
TN           = math.ceil(T/dt)
Dm           = ConstructDm(dx,N)
Dp           = ConstructDp(dx,N)

Ppinv, pminv = ConstructPpminv(dx,N)
Pp           = inv(Ppinv)
Pm           = inv(pminv)

E            = InitE(xp)
H            = InitH(xm)

Ener         = eps*E.dot(Pp.dot(E))+mu*H.dot(Pm.dot(H))
print(Ener)
for i in range(TN):
    Ener = eps*E.dot(Pp.dot(E))+mu*H.dot(Pm.dot(H))
    print(Ener)
    E    = E+dt*Dp.dot(H)
    H    = H+dt*Dm.dot(E)


