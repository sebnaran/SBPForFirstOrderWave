from scipy.sparse import csr_matrix,lil_matrix
from scipy.sparse.linalg import norm, inv
import numpy as np
import math
from Functions import ConstructDp, ConstructDm, ConstructPpminv


def test_ConstructDpDm():
    #This routine will test that, at least when N = 4, 
    #the finite difference operators are constructed well.
    N   = 4
    dx  = 2/N

    Dm  = ConstructDm(dx,N)
    TDm = csr_matrix([[-1,1,0,0,0],\
                      [-1,1,0,0,0],\
                      [0,-1,1,0,0],\
                      [0,0,-1,1,0],\
                      [0,0,0,-1,1],\
                      [0,0,0,-1,1]])/dx
    
    merr = norm(Dm-TDm)
    
    Dp  = ConstructDp(dx,N)
    TDp = csr_matrix([[-1,1,0,0,0,0],\
                      [0,-1,1,0,0,0],\
                      [0,0,-1,1,0,0],\
                      [0,0,0,-1,1,0],\
                      [0,0,0,0,-1,1]])/dx
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
    
    eps          = norm(Ppinv-TPpinv)+norm(Pminv-TPminv)
    assert eps < 0.0001


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
        
        err = abs(LB*LA-FB*FA-A.dot(Q.dot(B)))
        assert err < 0.0001

    for i in range(6):
        FA = FAlist[i] 
        FB = FBlist[i]
        LA = LAlist[i]
        LB = LBlist[i]
        
        A = np.array([FA,6,8,2,LA])
        B = np.array([FB,-3,-9,2,3,LB])
        
        err = abs(LB*LA-FB*FA-A.dot(Q.dot(B)))
        assert err < 0.0001



