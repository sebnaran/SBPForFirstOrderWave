from scipy.sparse import csr_matrix
from scipy.sparse import lil_matrix
from scipy.integrate import ode
import numpy as np
import math
from math import cos,sin,pi
##These functions relate to the initiation of the code and mesh
##building.

def StaggeredSpatialMesh(N,a,b):
    #This function creates a pair of staggered mesh of [-1,1]
    #It takes the number of subpartitions and returns the 
    #two grids alongside with the mesh size dx.
    if b < a:
        print('The interval is not well defined b < a')
        return 
    dx = (b-a)/N
    xp = [a+dx*i for i in range(N+1)]
    xm = [a+dx/2+dx*i for i in range(N)]
    xm = [a] + xm + [b]
    return xp,xm,dx

def InterfaceMesh(N,I):
    #It is assumed that a and b are lists of the points where there is an 
    #interface. This function will return two lists of meshes, one for each 
    #different material. These meshes will be collocated at the material 
    #interfaces and staggered in the interior of each material.

    n = len(I)

    
    xp = [[]]*(n-1)
    xm = [[]]*(n-1)
    for i in range(n-1):
        xp[i], xm[i], dx = StaggeredSpatialMesh(N,I[i],I[i+1])
    return xp, xm, dx        

    

def ConstructDm(dx,N):
    #Here we will build the finite difference operator from xplus    #grid to the xminus grid.
    Dm          = lil_matrix((N+2,N+1))
    Dm[0,0]     = -1
    Dm[0,1]     = 1
    Dm[1,0]     = -1
    Dm[1,1]     = 1
    Dm[2,0]     = -0.2
    Dm[2,1]     = -0.6
    Dm[2,2]     = 0.8

    Dm[N-1,N]   = 0.2
    Dm[N-1,N-1] = 0.6    
    Dm[N-1,N-2] = -0.8
    Dm[N,N]     = 1
    Dm[N,N-1]   = -1
    Dm[N+1,N]   = 1
    Dm[N+1,N-1] = -1

    for row in range(3,N-1):
        Dm[row,row-1] = -1
        Dm[row,row]   = 1

    Dm = Dm/dx
    return Dm.tocsr()
    
def ConstructDp(dx,N):
    #Here we will build the finite difference operator from 
    #xminus to the xplus grid.
    Dp = lil_matrix((N+1,N+2))
    
    Dp[0,0] = -1
    Dp[0,1] = 0.5
    Dp[0,2] = 0.5
    Dp[1,0] = -0.5
    Dp[1,1] = -0.25
    Dp[1,2] = 0.75
    for row in range(2,N-1):
        Dp[row,row]   = -1
        Dp[row,row+1] = 1
    Dp[N,N+1]   = 1
    Dp[N,N]     = -0.5
    Dp[N,N-1]   = -0.5
    Dp[N-1,N+1] = 0.5
    Dp[N-1,N]   = 0.25
    Dp[N-1,N-1] = -0.75
    Dp = Dp/dx
    return Dp.tocsr()



def ConstructPpminv(dx,N):
    #Here we build the inner product matrix on the xplus grid.
    Ppinv          = lil_matrix((N+1,N+1))
    Ppinv.setdiag(1)
    Ppinv[0,0]     = 2
    Ppinv[N,N]     = 2

    #And the inner product matrix on the xminus grid.
    Pminv          = lil_matrix((N+2,N+2))
    Pminv.setdiag(1)
    Pminv[0,0]     = 2
    Pminv[1,1]     = 4
    Pminv[2,2]     = 0.8
    Pminv[N+1,N+1] = 2
    Pminv[N,N]     = 4
    Pminv[N-1,N-1] = 0.8

    Ppinv  = Ppinv/dx
    Pminv  = Pminv/dx

    return Ppinv.tocsr(), Pminv.tocsr()

def NoIntInitE(x):
    return np.array([math.sin(math.pi*y) for y in x])
    
def NoIntInitH(x):
    return np.array([-math.cos(math.pi*y) for y in x])

def NoIntExactE(x,t):
    sol = [math.sin(math.pi*y)*(math.cos(math.pi*t)+\
                                math.sin(math.pi*t)) for y in x]
    return np.array(sol)

def NoIntExactH(x,t):
    sol = [math.cos(math.pi*y)*(math.sin(math.pi*t)-\
                                math.cos(math.pi*t)) for y in x]
    return np.array(sol)

def IntInitE(x):
    def f(xi):
        if xi<0:
            return cos(2*pi*xi)+2*sin(2*pi*xi)
        else:
            return cos(4*pi*xi)+sin(4*pi*xi)
    return np.array([f(y) for y in x])
            
def IntInitH(x):
    def f(xi):
        if xi<0:
            return -2*cos(2*pi*xi)+sin(2*pi*xi)
        else:
            return -2*cos(4*pi*xi)+2*sin(4*pi*xi)
    return np.array([f(y) for y in x])



def IntExactE(x,t):
    def f(xi,t):
        if xi<0:
            eval = cos(2*pi*t)*cos(2*pi*xi)+\
                 2*cos(2*pi*t)*sin(2*pi*xi)+\
                   sin(2*pi*t)*cos(2*pi*xi)+\
                 2*sin(2*pi*t)*sin(2*pi*xi)
            return eval
        else:
            eval = cos(2*pi*t)*cos(4*pi*xi)+\
                   cos(2*pi*t)*sin(4*pi*xi)+\
                   sin(2*pi*t)*cos(4*pi*xi)+\
                   sin(2*pi*t)*sin(4*pi*xi)
            return eval
    return np.array([f(y,t) for y in x])

def IntExactH(x,t):
    def f(xi,t):
        if xi<0:
            eval = -2*cos(2*pi*t)*cos(2*pi*xi)+\
                      cos(2*pi*t)*sin(2*pi*xi)+\
                    2*sin(2*pi*t)*cos(2*pi*xi)+\
                   -1*sin(2*pi*t)*sin(2*pi*xi)
            return eval
        else:
            eval = -2*cos(2*pi*t)*cos(4*pi*xi)+\
                    2*cos(2*pi*t)*sin(4*pi*xi)+\
                    2*sin(2*pi*t)*cos(4*pi*xi)+\
                   -2*sin(2*pi*t)*sin(4*pi*xi)
            return eval
    return np.array([f(y,t) for y in x])


def EoRetrieve(EH,N):
    return EH[0:N+1]

def EtRetrieve(EH,N):
    return EH[N+1:2*N+2]

def HoRetrieve(EH,N):
    return EH[2*N+2:3*N+4]

def HtRetrieve(EH,N):
    return EH[3*N+4:4*N+6]

def EoSet(EH,N,Eo):
    if len(Eo) == N+1:
        EH[0:N+1] = Eo
        return EH
    else:
        print('wrong size of Eo')
    return

def EtSet(EH,N,Et):
    if len(Et) == N+1:
        EH[N+1:2*N+2] = Et
        return EH
    else:
        print('wrong size of Et')
    return

def HoSet(EH,N,Ho):
    if len(Ho) == N+2:
        EH[2*N+2:3*N+4] = Ho
        return EH
    else:
        print('wrong size of Ho')
    return

def HtSet(EH,N,Ht):
    if len(Ht) == N+2:
        EH[3*N+4:4*N+6] = Ht
        return EH
    else:
        print('wrong size of Ht')
    return

def tDeriv(Eo,Et,Ho,Ht,t,eps,mu,Dp,Dm,Ppinv,Pminv,A1,A2,A3,A4,B1,B2,B3,B4,N):

    #Eo = EoRetrieve(EEH,N)
    #Et = EtRetrieve(EEH,N)
    #Ho = HoRetrieve(EEH,N)
    #Ht = HtRetrieve(EEH,N)

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

    #EEH = EoSet(EEH,N,TEo)
    #EEH = EtSet(EEH,N,TEt)
    #EEH = HoSet(EEH,N,THo)
    #EEH = HtSet(EEH,N,THt)
    return TEo,TEt,THo,THt

















#def ConstructPinv(dx,N):
#    #Here we build the matrix P
#    Pplusinv,Pminusinv = ConstructPpm(dx,N)
#    
#    Pinv = lil_matrix((2*N+3,2*N+3))
#    Pinv[0:N+2,0:N+2] = Pminusinv
#    Pinv[N+2:2*N+3,N+2:2*N+3] = Pplusinv
#    Pinv = Pinv.tocsr()
#    return Pinv
#
#
#def ConstructAtilde(N,dx):
#    #In this function we combine the two matrices Dplus and Dminus
#    Dminus              = ConstructDm(dx,N)
#    Dplus               = ConstructDp(dx,N)
#    A                   = lil_matrix((2*N+3,2*N+3))
#    A[0:N+1, N+1:2*N+3] = Dminus
#    A[N+1:2*N+3, 0:N+1] = Dplus
#    A                   = A.tocsr()  #Here we change the data structure to one thatis suited for efficient computation.
#    return A






