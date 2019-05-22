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


