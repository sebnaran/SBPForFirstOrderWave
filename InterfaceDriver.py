from scipy.sparse import csr_matrix,lil_matrix
from scipy.sparse.linalg import norm, inv
from scipy.integrate import odeint
from Functions import *
import math

N            = 1280
I            = [-1,0,1]
eps          = [1,4]
mu           = [1,1]
CFL          = 0.5
T            = 1 

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


EnerI = (eps[0]*Eo.dot(Pp.dot(Eo))+mu[0]*Ho.dot(Pm.dot(Ho))+\
        eps[1]*Et.dot(Pp.dot(Et))+mu[1]*Ht.dot(Pm.dot(Ht)))/2 



#EH = np.zeros(4*N+6)

#EH = EoSet(EH,N,Eo)
#EH = EtSet(EH,N,Et)
#EH = HoSet(EH,N,Ho)
#EH = HtSet(EH,N,Ht)

#RungeKutta4

def PtDeriv(Eo,Et,Ho,Ht):
    TEo,TEt,THo,THt = tDeriv(Eo,Et,Ho,Ht,0,eps,mu,Dp,Dm,Ppinv,Pminv,A1,A2,A3,A4,B1,B2,B3,B4,N)
    return TEo,TEt,THo,THt

for j in range(TN):
     
    tEo,tEt,tHo,tHt = PtDeriv(Eo,Et,Ho,Ht)

    Eok1 = dt*tEo
    Hok1 = dt*tHo
    Etk1 = dt*tEt   
    Htk1 = dt*tHt   
    
    tEo,tEt,tHo,tHt = PtDeriv(Eo+0.5*Eok1,Et+0.5*Etk1,Ho+0.5*Hok1,Ht+0.5*Htk1)

    Eok2 = dt*tEo
    Hok2 = dt*tHo
    Etk2 = dt*tEt 
    Htk2 = dt*tHt 

    tEo,tEt,tHo,tHt = PtDeriv(Eo+0.5*Eok2,Et+0.5*Etk2,Ho+0.5*Hok2,Ht+0.5*Htk2)

    Eok3 = dt*tEo
    Hok3 = dt*tHo
    Etk3 = dt*tEt
    Htk3 = dt*tHt

    tEo,tEt,tHo,tHt = PtDeriv(Eo+Eok3,Et+Etk3,Ho+Hok3,Ht+Htk3)

    Eok4 = dt*tEo
    Hok4 = dt*tHo
    Etk4 = dt*tEt
    Htk4 = dt*tHt


    Eo = Eo+(Eok1+2*Eok2+2*Eok3+Eok4)/6 
    Ho = Ho+(Hok1+2*Hok2+2*Hok3+Hok4)/6
    Et = Et+(Etk1+2*Etk2+2*Etk3+Etk4)/6
    Ht = Ht+(Htk1+2*Htk2+2*Htk3+Htk4)/6

#    Ener = eps[0]*Eo.dot( Pp.dot(Eo) )+mu[0]*Ho.dot( Pm.dot(Ho) )+\
#           eps[1]*Et.dot( Pp.dot(Et) )+mu[1]*Ht.dot( Pm.dot(Ht) )/2
#    print(Ener)




#for j in range(TN):
#    tEo,tEt,tHo,tHt = tDeriv(Eo,Et,Ho,Ht,0,eps,mu,Dp,Dm,Ppinv,Pminv,A1,A2,A3,A4,B1,B2,B3,B4,N)
#
#    Eo = Eo+dt*tEo
#    Ho = Ho+dt*tHo
#    Et = Et+dt*tEt
#    Ht = Ht+dt*tHt
##
#    Ener = eps[0]*Eo.dot( Pp.dot(Eo) )+mu[0]*Ho.dot( Pm.dot(Ho) )+\
#       eps[1]*Et.dot( Pp.dot(Et) )+mu[1]*Ht.dot( Pm.dot(Ht) )/2
#    print(Ener)





#This uses a time integrator from scypy.
#t  = np.linspace(0,T,TN)
#EH = odeint(tDeriv,EH,t)
#EH = EH[TN-1,:]


#Eo = EoRetrieve(EH,N)
#Et = EtRetrieve(EH,N)
#Ho = HoRetrieve(EH,N)
#Ht = HtRetrieve(EH,N)


EnerF = eps[0]*Eo.dot(Pp.dot(Eo))+mu[0]*Ho.dot(Pm.dot(Ho)) +\
       eps[1]*Et.dot(Pp.dot(Et))+mu[1]*Ht.dot(Pm.dot(Ht))/2
print(EnerF)

ExEo = IntExactE(xp[0],T)
ExEt = IntExactE(xp[1],T)
ExHo = IntExactH(xm[0],T)
ExHt = IntExactH(xm[1],T)

Eoerr = ExEo-Eo
Eterr = ExEt-Et
Hoerr = ExHo-Ho
Hterr = ExHt-Ht

Eerr = Eoerr.dot(Pp.dot(Eoerr))+Eterr.dot(Pp.dot(Eterr))
Herr = Hoerr.dot(Pm.dot(Hoerr))+Hterr.dot(Pm.dot(Hterr))

Eerr = math.sqrt(Eerr)
Herr = math.sqrt(Herr)

print(Eerr)
print(Herr)


