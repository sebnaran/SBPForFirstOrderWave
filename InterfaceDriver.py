from scipy.sparse import csr_matrix,lil_matrix
from scipy.sparse.linalg import norm, inv
from scipy.integrate import odeint
from Functions import *
import math

N            = 4
I            = [-1,0,1]
eps          = [1,1]
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

Pp           = inv(Ppinv)
Pm           = inv(Pminv)

A1      = np.zeros(N+1)
A1[N]   = -0.5
B1      = np.zeros(N+1)
B1[0]   = 0.5

A2      = np.zeros(N+2)
A2[N+1] = -0.5
B2      = np.zeros(N+2)
B2[0]   = 0.5

A3      = np.zeros(N+1)
A3[0]   = 0.5
B3      = np.zeros(N+1)
B3[N]   = -0.5

A4      = np.zeros(N+2)
A4[0]   = 0.5
B4      = np.zeros(N+2)
B4[N+1] = -0.5

Eo = IntInitE(xp[0])
Ho = IntInitH(xm[0])
Et = IntInitE(xp[1])
Ht = IntInitH(xm[1])


Ener = (eps[0]*Eo.dot(Pp.dot(Eo))+mu[0]*Ho.dot(Pm.dot(Ho))+\
        eps[1]*Et.dot(Pp.dot(Et))+mu[1]*Ht.dot(Pm.dot(Ht)))/2 

print(Ener)

EH = np.zeros(4*N+6)

EH = EoSet(EH,N,Eo)
EH = EtSet(EH,N,Et)
EH = HoSet(EH,N,Ho)
EH = HtSet(EH,N,Ht)

def tDeriv(EH,t):
    
    Eo = EoRetrieve(EH,N)
    Et = EtRetrieve(EH,N)
    Ho = HoRetrieve(EH,N)
    Ht = HtRetrieve(EH,N)

    HoN   = Ho[len(Ho)-1]
    Ht0   = Ht[0]
    EoN   = Eo[len(Eo)-1]
    Et0   = Et[0]

    Ho0   = Ho[0]
    HtN   = Ht[len(Ht)-1]
    Eo0   = Eo[0]
    EtN   = Et[len(Et)-1]

    TEo = (Dp.dot(Ho)+\
           Ppinv.dot(A1)*(HoN-Ht0)+\
           Ppinv.dot(B1)*(Ho0-HtN))/eps[0]
    THo = (Dm.dot(Eo)+\
           Pminv.dot(A2)*(EoN-Et0)+\
           Pminv.dot(B2)*(Eo0-EtN))/mu[0]
    TEt = (Dp.dot(Ht)+\
           Ppinv.dot(A3)*(Ht0-HoN)+\
           Ppinv.dot(B3)*(HtN-Ho0))/eps[1]
    THt = (Dm.dot(Et)+\
           Pminv.dot(A4)*(Et0-EoN)+\
           Pminv.dot(B4)*(EtN-Eo0))/mu[1]

    EH = EoSet(EH,N,TEo)
    EH = EtSet(EH,N,TEt)
    EH = HoSet(EH,N,THo)
    EH = HtSet(EH,N,THt)
    return EH

#RungeKutta4
#for j in range(TN):
#
#    k1 = dt*tDeriv(EH,0)
#    k2 = dt*tDeriv(EH+0.5*k1,0)
#    k3 = dt*tDeriv(EH+0.5*k2,0)
#    k4 = dt*tDeriv(EH+k3,0)
#
#    EH = EH + (k1+2*k2+2*k3+k4)/6

#Euler Step

#for j in range(TN):
#    EH = EH+dt*tDeriv(EH,0)




#This uses a time integrator from scypy.
t  = np.linspace(0,T,TN)
EH = odeint(tDeriv,EH,t)
EH = EH[TN-1,:]


Eo = EoRetrieve(EH,N)
Et = EtRetrieve(EH,N)
Ho = HoRetrieve(EH,N)
Ht = HtRetrieve(EH,N)
###
Ener = eps[0]*Eo.dot( Pp.dot(Eo) )+mu[0]*Ho.dot( Pm.dot(Ho) )+\
       eps[1]*Et.dot( Pp.dot(Et) )+mu[1]*Ht.dot( Pm.dot(Ht) )/2
print(Ener)



ExEo = IntExactE(xp[0],T)
ExEt = IntExactE(xp[1],T)
ExHo = IntExactH(xm[0],T)
ExHt = IntExactH(xm[1],T)

Eoerr = ExEo-Eo
Eterr = ExEt-Et
Hoerr = ExHo-Ho
Hterr = ExHt-Ht

err = Eoerr.dot(Pp.dot(Eoerr))+Eterr.dot(Pp.dot(Eterr))+\
      Hoerr.dot(Pm.dot(Hoerr))+Hterr.dot(Pm.dot(Hterr))

err = math.sqrt(err)

print(err)

