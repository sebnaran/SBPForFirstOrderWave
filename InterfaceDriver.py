from scipy.sparse import csr_matrix,lil_matrix
from scipy.sparse.linalg import norm, inv
from Functions import *
import math

N            = 10 
I            = [-1,0,1]
eps          = [1,2]
mu           = [1,2]
CFL          = 0.5
T            = 1

xp,xm,dx     = InterfaceMesh(N,I)
dt           = CFL*dx
TN           = math.ceil(T/dt)
NI           = len(I)-1


Dm           = ConstructDm(dx,N)
Dp           = ConstructDp(dx,N)
Ppinv, Pminv = ConstructPpminv(dx,N)

Dm           = [ Dm ]*NI
Dp           = [ Dm ]*NI
Pp           = [ inv(Ppinv) ]*NI
Pm           = [ inv(Pminv) ]*NI
Ppinv        = [ Ppinv ]*NI
Pminv        = [ Pminv ]*NI
AD1          = [ csr_matrix( ([-0.5],([N],[0]) ),shape=(N+1,1))  ]*NI
AD2          = [ csr_matrix( ([ 0.5],([0],[0]) ),shape=(N+2,1))  ]*NI

E            = [ [] ]*NI
H            = [ [] ]*NI
Ener         = 0
for i in range(NI):

    E[i]               = InitE(xp[i])
    H[i]               = InitH(xm[i])
    Ener               = Ener+eps[i]*E[i].dot( Pp[i].dot(E[i]) )\
                             +mu[i]*H[i].dot(  Pm[i].dot(H[i]) )

print(Ener)



This time integration uses an Euler Step.
for j in range(TN):
    
    E    = E+dt*Dp.dot(H)
    H    = H+dt*Dm.dot(E)
    
    Ener = 0
    
    for i in range(IN):
        Ener = Ener+eps[i]*E[i].dot( Pp[i].dot(E[i]) )\
                             +mu[i]*H[i].dot(  Pm[i].dot(H[i]) )

    print(Ener)


#This time integration uses an RK4 method.
#for i in range(TN):
#    Ener = eps*E.dot(Pp.dot(E))+mu*H.dot(Pm.dot(H))
#    print(Ener)
#    
#    Ek1 = dt*Dp.dot(H)
#    Hk1 = dt*Dm.dot(E)
#    
#    Ek2 = dt*Dp.dot(H+Hk1/2)
#    Hk2 = dt*Dm.dot(E+Ek1/2)
#
#    Ek3 = dt*Dp.dot(H+Hk2/2)
#    Hk3 = dt*Dm.dot(E+Ek2/2)
#
#    Ek4 = dt*Dp.dot(H+Hk3)
#    Hk4 = dt*Dm.dot(E+Ek3)
#
#    E   = E+(Ek1+2*Ek2+2*Ek3+Ek4)/6
#    H   = H+(Hk1+2*Hk2+2*Hk3+Hk4)/6





#print('Boundary Values')
#print(E[0])
#print(E[len(E)-1])
#print(H[0])
#print(H[len(H)-1])
