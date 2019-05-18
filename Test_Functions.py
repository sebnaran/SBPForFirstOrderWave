from scipy.sparse import csr_matrix,lil_matrix
from scipy.sparse.linalg import norm
import numpy as np
import math
from Functions import ConstructDp, ConstructDm, ConstructPpminv


def test_ConstructDpDm():
    dx = 0.5
    Dp = ConstructDp(dx,4)
    Dm = ConstructDm(dx,4)
    TDp = csr_matrix([[-1,1,0,0,0],\
                      [-1,1,0,0,0],\
                      [0,-1,1,0,0],\
                      [0,0,-1,1,0],\
                      [0,0,0,-1,1],\
                      [0,0,0,-1,1]])
    TDp = TDp/dx
    eps = norm(Dp-TDp)
    print(eps)
    assert eps < 0.001



#test_ConstructDp():
#    dx = 2/N
#    Dp = ConstructDp(dx,4)
#    Dm = ConstructDm(dx,4)
#    PPlusinv, Pminusinv = ConstructPpminv(dx,4)
#    Pp = inv(PPlusinv)
#    Pm = inv(Pminusinv)
