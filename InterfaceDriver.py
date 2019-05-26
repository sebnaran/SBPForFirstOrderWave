from scipy.sparse import csr_matrix,lil_matrix
from scipy.sparse.linalg import norm, inv
from Functions import *
import math

N            = 1000 
I            = [-1,0,1]
eps          = [1,2]
mu           = [1,2]
CFL          = 0.1
T            = 1

xp,xm,dx     = InterfaceMesh(N,I)
dt           = CFL*dx
TN           = math.ceil(T/dt)
NI           = len(I)-1


Dm           = ConstructDm(dx,N)
Dp           = ConstructDp(dx,N)
Ppinv, Pminv = ConstructPpminv(dx,N)

Dm           = [ Dm ]*NI
Dp           = [ Dp ]*NI
Pp           = [ inv(Ppinv) ]*NI
Pm           = [ inv(Pminv) ]*NI
Ppinv        = [ Ppinv ]*NI
Pminv        = [ Pminv ]*NI


AD11         = np.zeros(N+1)
AD11[N]      = -0.5
AD12         = np.zeros(N+2)
AD12[N+1]    = -0.5

AD21         = np.zeros(N+1)
AD21[0]      = 0.5
AD22         = np.zeros(N+2)
AD22[0]      = 0.5


#AD1          = [ csr_matrix( ([-0.5],([N],[0]) ),shape=(N+1,1))  ]*NI
#AD2          = [ csr_matrix( ([ 0.5],([0],[0]) ),shape=(N+2,1))  ]*NI

E            = [ [] ]*NI
H            = [ [] ]*NI
Ener         = 0
for i in range(NI):

    E[i]               = InitE(xp[i])
    H[i]               = InitH(xm[i])
    Ener               = Ener+eps[i]*E[i].dot( Pp[i].dot(E[i]) )\
                             +mu[i]*H[i].dot(  Pm[i].dot(H[i]) )

print(Ener)



#This time integration uses an Euler Step.
for j in range(TN):
    
    HN   = H[0][len(H[0])-1]
    H0   = H[1][0]

    EN   = E[0][len(E[0])-1]
    E0   = E[1][0]

    
    E[0] = E[0]+dt*( Dp[0].dot(H[0])+\
                     Ppinv[0].dot(AD11)*( HN-H0 ) )\
                     /eps[0]
            
    H[0] = H[0]+dt*( Dm[0].dot(E[0])+\
                   ( EN-E0 )*Pminv[0].dot(AD12) )\
                    /mu[0]      
    
    E[1] = E[1]+dt*( Dp[1].dot(H[1])+\
                    ( H0-HN )*Ppinv[1].dot(AD21) )\
                    /eps[1]

    H[1] = H[1]+dt*( Dm[1].dot(E[1])+\
                   ( E0-EN )*Pminv[0].dot(AD22) )\
                    /mu[1]
    
    Ener = 0
    for i in range(NI):
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
