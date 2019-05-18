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

    Dp  = ConstructDp(dx,N)
    TDp = csr_matrix([[-1,1,0,0,0],\
                      [-1,1,0,0,0],\
                      [0,-1,1,0,0],\
                      [0,0,-1,1,0],\
                      [0,0,0,-1,1],\
                      [0,0,0,-1,1]])/dx
    
    perr = norm(Dp-TDp)
    
    Dm  = ConstructDm(dx,N)
    TDm = csr_matrix([[-1,1,0,0,0,0],\
                      [0,-1,1,0,0,0],\
                      [0,0,-1,1,0,0],\
                      [0,0,0,-1,1,0],\
                      [0,0,0,0,-1,1]])/dx
    merr = norm(Dm-TDm)

    eps = perr+merr
    assert eps < 0.001

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
    assert eps < 0.001


#def CheckSBPProperty():
#    #This routine will check that the finite difference 
#    #operators and the inner product matrices satisfy the
#    #SBP property at least when N = 4
#
#    N  = 4
#    dx = 2/N
#
#    Dp           = ConstructDp(dx,N)
#    Dm           = ConstructDm(dx,N)
#    Ppinv, Pminv = ConstructPpminv(dx,N)
#    Pp           = inv(Ppinv)
#    Pm           = inv(Pminv)
#
#    Qp = Pp.dot(Dp)
#    Qm = Pm.dot(Dm)
#
#    Q = Qp + Qm.transpose()
#    assert 1 < 2
