from scipy.sparse import csr_matrix,lil_matrix
from scipy.sparse.linalg import norm
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
    
    pluserror = norm(Dp-TDp)
    
    Dm  = ConstructDm(dx,N)
    TDm = csr_matrix([[-1,1,0,0,0,0],\
                      [0,-1,1,0,0,0],\
                      [0,0,-1,1,0,0],\
                      [0,0,0,-1,1,0],\
                      [0,0,0,0,-1,1]])/dx
    minuserror = norm(Dm-TDm)

    eps = pluserror+minuserror
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

#test_ConstructDp():
#   #    Dp = ConstructDp(dx,4)
#    Dm = ConstructDm(dx,4)
#
#    Pp = inv(PPlusinv)
#    Pm = inv(Pminusinv)
