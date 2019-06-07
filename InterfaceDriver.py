from scipy.sparse import csr_matrix,lil_matrix
from scipy.sparse.linalg import norm, inv
from scipy.integrate import odeint
from Functions import *
import math

N            = 10000 
I            = [-1,0,1]
eps          = [1,2]
mu           = [1,1]
CFL          = 0.05
T            = 0.1 

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

    E[i]               = IntInitE(xp[i])
    H[i]               = IntInitH(xm[i])
    Ener               = Ener+eps[i]*E[i].dot( Pp[i].dot(E[i]) )\
                              +mu[i]*H[i].dot( Pm[i].dot(H[i]) )

print(Ener)


#for j in range(TN):
#    
#    HN   = H[0][len(H[0])-1]
#    H0   = H[1][0]
#
#    EN   = E[0][len(E[0])-1]
#    E0   = E[1][0]
#
#    
#    E[0] = E[0]+(dt/eps[0])*( Dp[0].dot(H[0])+\
#                             ( HN-H0 )*Ppinv[0].dot(AD11) )
#            
#    H[0] = H[0]+(dt/mu[0])*( Dm[0].dot(E[0])+\
#                   ( EN-E0 )*Pminv[0].dot(AD12) )     
#    
#    E[1] = E[1]+(dt/eps[1])*( Dp[1].dot(H[1])+\
#                    ( H0-HN )*Ppinv[1].dot(AD21) )
#
#    H[1] = H[1]+(dt/mu[1])*( Dm[1].dot(E[1])+\
#                   ( E0-EN )*Pminv[0].dot(AD22) )
# 
#    for i in range(NI):
#        Ener = Ener+eps[i]*E[i].dot( Pp[i].dot(E[i]) )\
#                             +mu[i]*H[i].dot(  Pm[i].dot(H[i]) )
#
#    print(Ener)









EH = np.zeros(4*N+6)

EH = EoSet(EH,N,E[0])
EH = EtSet(EH,N,E[1])
EH = HoSet(EH,N,H[0])
EH = HtSet(EH,N,H[1])


def Func(EH,t):
    Eo = EoRetrieve(EH,N)
    Et = EtRetrieve(EH,N)
    Ho = HoRetrieve(EH,N)
    Ht = HtRetrieve(EH,N) 

    HN   = Ho[len(Ho)-1]
    H0   = Ht[0]
    EN   = Eo[len(Eo)-1]
    E0   = Et[0]
    
    TEo = (Dp[0].dot(Ho)+Ppinv[0].dot(AD11)*(HN-H0))/eps[0]
    THo = (Dm[0].dot(Eo)+Pminv[0].dot(AD12)*(EN-E0))/mu[0]
    TEt = (Dp[1].dot(Ht)+Ppinv[1].dot(AD21)*(H0-HN))/eps[1]
    THt = (Dm[1].dot(Et)+Pminv[1].dot(AD22)*(E0-EN))/mu[1]

    EH = EoSet(EH,N,TEo)
    EH = EtSet(EH,N,TEt)
    EH = HoSet(EH,N,THo)
    EH = HtSet(EH,N,THt)
    return EH

#for j in range(TN):
#    EH = EH+dt*Func(EH,0)
#
#    Eo = EoRetrieve(EH,N)
#    Et = EtRetrieve(EH,N)
#    Ho = HoRetrieve(EH,N)
#    Ht = HtRetrieve(EH,N)
#
#    Ener = eps[0]*Eo.dot( Pp[0].dot(Eo) )+mu[0]*Ho.dot( Pm[0].dot(Ho) )+\
#       eps[1]*Et.dot( Pp[1].dot(Et) )+mu[1]*Ht.dot( Pm[1].dot(Ht) )
#    print(Ener)


   
#    Ener = 0
#    for i in range(NI):
#        Ener = Ener+eps[i]*E[i].dot( Pp[i].dot(E[i]) )\
#                             +mu[i]*H[i].dot(  Pm[i].dot(H[i]) )
#
#    print(Ener)



#print(err)

#for j in range(TN):
#
#    k1 = dt*Func(EH,0)
#    k2 = dt*Func(EH+0.5*k1,0)
#    k3 = dt*Func(EH+0.5*k2,0)
#    k4 = dt*Func(EH+k3,0)
#    
#    EH = EH + (k1+2*k2+2*k3+k4)/6
#    Eo = EoRetrieve(EH,N)
#    Et = EtRetrieve(EH,N)
#    Ho = HoRetrieve(EH,N)
#    Ht = HtRetrieve(EH,N)
#
#    Ener = eps[0]*Eo.dot( Pp[0].dot(Eo) )+mu[0]*Ho.dot( Pm[0].dot(Ho) )+\
#       eps[1]*Et.dot( Pp[1].dot(Et) )+mu[1]*Ht.dot( Pm[1].dot(Ht) )
#    print(Ener)
#

















t  = np.linspace(0,T,TN)
#
EH = odeint(Func,EH,t)
#

EH = EH[TN-1,:]
Eo = EoRetrieve(EH,N)
Et = EtRetrieve(EH,N)
Ho = HoRetrieve(EH,N)
Ht = HtRetrieve(EH,N)
#
Ener = eps[0]*Eo.dot( Pp[0].dot(Eo) )+mu[0]*Ho.dot( Pm[0].dot(Ho) )+\
       eps[1]*Et.dot( Pp[1].dot(Et) )+mu[1]*Ht.dot( Pm[1].dot(Ht) )
print(Ener)












#
#
#
##This time integration uses an Euler Step.
#for j in range(TN):
#    
#    HN   = H[0][len(H[0])-1]
#    H0   = H[1][0]
#
#    EN   = E[0][len(E[0])-1]
#    E0   = E[1][0]
#
#    
#    E[0] = E[0]+dt*( Dp[0].dot(H[0])+\
#                     Ppinv[0].dot(AD11)*( HN-H0 ) )\
#                     /eps[0]
#            
#    H[0] = H[0]+dt*( Dm[0].dot(E[0])+\
#                   ( EN-E0 )*Pminv[0].dot(AD12) )\
#                    /mu[0]      
#    
#    E[1] = E[1]+dt*( Dp[1].dot(H[1])+\
#                    ( H0-HN )*Ppinv[1].dot(AD21) )\
#                    /eps[1]
#
#    H[1] = H[1]+dt*( Dm[1].dot(E[1])+\
#                   ( E0-EN )*Pminv[0].dot(AD22) )\
#                    /mu[1]
#    
#    Ener = 0
#    for i in range(NI):
#        Ener = Ener+eps[i]*E[i].dot( Pp[i].dot(E[i]) )\
#                             +mu[i]*H[i].dot(  Pm[i].dot(H[i]) )
#
#    print(Ener)


ExE1 = IntExactE(xp[0],T)
ExE2 = IntExactE(xp[1],T)
ExH1 = IntExactH(xm[0],T)
ExH2 = IntExactH(xm[1],T)

E1err = ExE1-E[0]
E2err = ExE2-E[1]
H1err = ExH1-H[0]
H2err = ExH2-H[1]

err = E1err.dot(Pp[0].dot(E1err))+E2err.dot(Pp[1].dot(E2err))+\
      H1err.dot(Pm[0].dot(H1err))+H2err.dot(Pm[1].dot(H2err))
err = math.sqrt(err)

print(err)




















