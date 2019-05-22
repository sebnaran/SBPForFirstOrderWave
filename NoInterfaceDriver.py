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

#This time integration uses an Euler Step.
#for i in range(TN):
#    Ener = eps*E.dot(Pp.dot(E))+mu*H.dot(Pm.dot(H))
#    print(Ener)
#    E    = E+dt*Dp.dot(H)
#    H    = H+dt*Dm.dot(E)


#This time integration uses an RK4 method.
for i in range(TN):
    Ener = eps*E.dot(Pp.dot(E))+mu*H.dot(Pm.dot(H))
    print(Ener)
    
    Ek1 = dt*Dp.dot(H)
    Hk1 = dt*Dm.dot(E)
    
    Ek2 = dt*Dp.dot(H+Hk1/2)
    Hk2 = dt*Dm.dot(E+Ek1/2)

    Ek3 = dt*Dp.dot(H+Hk2/2)
    Hk3 = dt*Dm.dot(E+Ek2/2)

    Ek4 = dt*Dp.dot(H+Hk3)
    Hk4 = dt*Dm.dot(E+Ek3)

    E   = E+(Ek1+2*Ek2+2*Ek3+Ek4)/6
    H   = H+(Hk1+2*Hk2+2*Hk3+Hk4)/6





#print('Boundary Values')
#print(E[0])
#print(E[len(E)-1])
#print(H[0])
#print(H[len(H)-1])

