from scipy.sparse import csr_matrix,lil_matrix
from scipy.sparse.linalg import norm, inv
import numpy as np
import math
from Functions import *

#These functions test the SBP construction of the SBP matrices and the SBP
#Property
def test_ConstructDpDm():
    #This routine will test that, at least when N = 4, 
    #the finite difference operators are constructed well.
    N   = 5
    dx  = 2/N

    Dm  = ConstructDm(dx,N)
    TDm = csr_matrix([[-1,1,0,0,0,0],\
                      [-1,1,0,0,0,0],\
                      [-0.2,-0.6,0.8,0,0,0],\
                      [0,0,-1,1,0,0],\
                      [0,0,0,-0.8,0.6,0.2],\
                      [0,0,0,0,-1,1],\
                      [0,0,0,0,-1,1]])/dx
    
    merr = norm(Dm-TDm)
    
    Dp  = ConstructDp(dx,N)
    TDp = csr_matrix([[-1,0.5,0.5,0,0,0,0],\
                      [-0.5,-0.25,0.75,0,0,0,0],\
                      [0,0,-1,1,0,0,0],\
                      [0,0,0,-1,1,0,0],\
                      [0,0,0,0,-0.75,0.25,0.5],\
                      [0,0,0,0,-0.5,-0.5,1]])/dx
    perr = norm(Dp-TDp)

    eps = perr+merr
    assert eps < 0.0001

def test_ConstructPpminv():
    #This routine will test if the inverse of the
    #inner product matrices are constructed correctly
    N            = 4
    dx           = 2/N
    Ppinv, Pminv = ConstructPpminv(dx,N)
    TPpinv       = csr_matrix([[2,0,0,0,0],\
                               [0,1,0,0,0],\
                               [0,0,1,0,0],\
                               [0,0,0,1,0],\
                               [0,0,0,0,2]])/dx

    TPminv       = csr_matrix([[2,0,0,0,0,0],\
                               [0,4,0,0,0,0],\
                               [0,0,0.8,0,0,0],\
                               [0,0,0,0.8,0,0],\
                               [0,0,0,0,4,0],\
                               [0,0,0,0,0,2]])/dx
    
    eps1          = norm(Ppinv-TPpinv)+norm(Pminv-TPminv)
    assert eps1 < 0.0001
    
    Pp  = inv(Ppinv)
    TPp = csr_matrix([[0.5,0,0,0,0],\
                      [0,1,0,0,0],\
                      [0,0,1,0,0],\
                      [0,0,0,1,0],\
                      [0,0,0,0,0.5]])*dx
    
    Pm  = inv(Pminv)
    TPm = csr_matrix([[0.5,0,0,0,0,0],\
                      [0,0.25,0,0,0,0],\
                      [0,0,1.25,0,0,0],\
                      [0,0,0,1.25,0,0],\
                      [0,0,0,0,0.25,0],\
                      [0,0,0,0,0,0.5]])*dx
    
    eps2 = norm(Pp-TPp)+norm(Pm-TPm)
    assert eps2 < 0.0001


def test_SBPProterty():
    #This routine will check that the finite difference 
    #operators and the inner product matrices satisfy the
    #SBP property at least when N = 4

    N  = 4 
    dx = 2/N

    Dp           = ConstructDp(dx,N)
    Dm           = ConstructDm(dx,N)
    Ppinv, Pminv = ConstructPpminv(dx,N)
    Pp           = inv(Ppinv)
    Pm           = inv(Pminv)

##    Pm           = csr_matrix([[0,0,0,0,0,0],\
##                               [0,1,0,0,0,0],\
##                               [0,0,1,0,0,0],\
##                               [0,0,0,1,0,0],\
##                               [0,0,0,0,1,0],\
##                               [0,0,0,0,0,0]])*dx
## 
    Qp = Pp.dot(Dp)
    Qm = Pm.dot(Dm)

    Q = Qp + Qm.transpose()
    
    FAlist = [-3,-2,-1,0,1,2,3]
    LAlist = [-2,-1,0,1,2,3,-3]

    FBlist = [-6,-4,-2,0,2,4,6]
    LBlist = [6,-6,-4,-2,0,2,4]

    for i in range(6):
        FA = FAlist[i] 
        FB = FBlist[i]
        LA = LAlist[i]
        LB = LBlist[i]
        
        A = np.array([FA,0,0,0,LA])
        B = np.array([FB,0,0,0,0,LB])
        #A = np.array([FA,0,0,0,0,0,0,0,LA])
        #B = np.array([FB,0,0,0,0,0,0,0,0,LB])
        
        err = abs(LB*LA-FB*FA-A.dot(Q.dot(B)))
        assert err < 0.0001

    for i in range(6):
        FA = FAlist[i] 
        FB = FBlist[i]
        LA = LAlist[i]
        LB = LBlist[i]
        
        A = np.array([FA,6,8,2,LA])
        B = np.array([FB,-3,-9,2,3,LB])
        #A = np.array([FA,6,8,2,2,2,2,2,LA])
        #B = np.array([FB,-3,-9,2,3,3,5,7,9,LB])
        
        err = abs(LB*LA-FB*FA-A.dot(Q.dot(B)))
        assert err < 0.0001

#Here we test the construction of the mesh

def test_InterfaceMesh():
    N          = 2
    I          = [-1,0,1]
    xp, xm, dx = InterfaceMesh(N,I)
    
    Txp        = np.array([[-1,-0.5,0],[0,0.5,1]])
    Txm        = np.array([[-1,-0.75,-0.25,0],[0,0.25,0.75,1]])
    xp         = np.array(xp)
    xm         = np.array(xm)
    
    eps        = np.linalg.norm(xp-Txp) + np.linalg.norm(xm-Txm)
    #print(xp-Txp)
    #print(xm-Txm)
    #print(xp)
    #print(xm)
    assert eps < 0.0001

def test_SettingRetrieving():

    N  = 5
    EH = np.zeros(4*N+6)

    Eo = np.array([0,1,2,3,4,5])
    Et = np.array([6,7,8,9,10,11])

    Ho = np.array([12,13,14,15,16,17,18])
    Ht = np.array([19,20,21,22,23,24,25])
    
    EH = EoSet(EH,N,Eo)
    EH = EtSet(EH,N,Et)
    EH = HoSet(EH,N,Ho)
    EH = HtSet(EH,N,Ht)

    TEo = EoRetrieve(EH,N)
    TEt = EtRetrieve(EH,N)
    THo = HoRetrieve(EH,N)
    THt = HtRetrieve(EH,N)
    
    err = np.linalg.norm(Eo-TEo)+\
          np.linalg.norm(Et-TEt)+\
          np.linalg.norm(Ho-THo)+\
          np.linalg.norm(Ht-THt)
    assert err<0.0001
   

def test_Aij():
    N         = 5
    AD11      = np.zeros(N+1)
    AD11[N]   = -0.5
    AD12      = np.zeros(N+2)
    AD12[N+1] = -0.5

    AD21         = np.zeros(N+1)
    AD21[0]      = 0.5
    AD22         = np.zeros(N+2)
    AD22[0]      = 0.5
    
    
    E = np.array([1,2,3,4,5,6])
    H = np.array([10,11,12,13,14,15,16])

    assert E.dot(AD11)+3   < 0.001
    assert E.dot(AD21)-0.5 < 0.001
    assert H.dot(AD12)+8   < 0.001
    assert H.dot(AD22)-5   < 0.001

def test_TimeDeriv():

    N            = 100  
    I            = [-1,0,1]
    eps          = [1,2]
    mu           = [1,1]

    xp,xm,dx     = InterfaceMesh(N,I)
    NI           = len(I)-1


    Dm           = ConstructDm(dx,N)
    Dp           = ConstructDp(dx,N)
    Ppinv, Pminv = ConstructPpminv(dx,N)

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


    E            = [ [] ]*NI
    H            = [ [] ]*NI
    for i in range(NI):

        E[i]               = IntInitE(xp[i])
        H[i]               = IntInitH(xm[i])

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

    EH = np.zeros(4*N+6)

    EH = EoSet(EH,N,E[0])
    EH = EtSet(EH,N,E[1])
    EH = HoSet(EH,N,H[0])
    EH = HtSet(EH,N,H[1])
    
    TEo = EoRetrieve(EH,N)
    TEt = EtRetrieve(EH,N)
    THo = HoRetrieve(EH,N)
    THt = HtRetrieve(EH,N)
    
    print( Dm[0].dot(TEo) )
    print( Dm[0].dot(E[0]) )
    EH = Func(EH,0)

    TEo = EoRetrieve(EH,N)
    TEt = EtRetrieve(EH,N)
    THo = HoRetrieve(EH,N)
    THt = HtRetrieve(EH,N)

    HN   = H[0][len(H[0])-1]
    H0   = H[1][0]
    EN   = E[0][len(E[0])-1]
    E0   = E[1][0]
    
    
    NEo = (Dp[0].dot(H[0])+Ppinv[0].dot(AD11)*(HN-H0))/eps[0]
    NHo = (Dm[0].dot(E[0])+Pminv[0].dot(AD12)*(EN-E0))/mu[0]
    NEt = (Dp[1].dot(H[1])+Ppinv[1].dot(AD21)*(H0-HN))/eps[1]
    NHt = (Dm[1].dot(E[1])+Pminv[1].dot(AD22)*(E0-EN))/mu[1]
    


    errEo = np.linalg.norm(TEo-NEo)
    errEt = np.linalg.norm(TEt-NEt)
    errHo = np.linalg.norm(THo-NHo)
    errHt = np.linalg.norm(THt-NHt)
    
    
    assert errEo < 0.001
    assert errEt < 0.001
    assert errHo < 0.001
    assert errHt < 0.001



           







