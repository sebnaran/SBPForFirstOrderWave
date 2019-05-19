from scipy.sparse import csr_matrix
from scipy.sparse import lil_matrix
import numpy as np
import math
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
    xm = [a] + xminus + [b]
    return xplus,xminus,dx

def ConstructDm(dx,N):
    #Here we will build the finite difference operator from xplus    #grid to the xminus grid.
    Dm          = lil_matrix((N+2,N+1))
    Dm[0,0]     = -1
    Dm[0,1]     = 1
    Dm[N+1,N]   = 1
    Dm[N+1,N-1] = -1

    for row in range(1,N+1):
        Dm[row,row-1] = -1
        Dm[row,row]   = 1

    Dm = Dm/dx
    return Dm.tocsr()
    
def ConstructDp(dx,N):
    #Here we will build the finite difference operator from 
    #xminus to the xplus grid.
    Dp = lil_matrix((N+1,N+2))
    
    Dp[0,0] = -2
    Dp[0,1] = 2
    for row in range(1,N):
        Dp[row,row]   = -1
        Dp[row,row+1] = 1
    Dp[N,N+1] = 2
    Dp[N,N]   = -2
    Dp = Dp/dx
    return Dp.tocsr()

def ConstructAtilde(N,dx):
    #In this function we combine the two matrices Dplus and Dminus
    Dminus              = ConstructDm(dx,N)
    Dplus               = ConstructDp(dx,N)
    A = lil_matrix((2*N+3,2*N+3))
    A[0:N+1, N+1:2*N+3] = Dminus
    A[N+1:2*N+3, 0:N+1] = Dplus
    A = A.tocsr()  #Here we change the data structure to one thatis suited for efficient computation.
    return A


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

def ConstructPinv(dx,N):
    #Here we build the matrix P
    Pplusinv,Pminusinv = ConstructPpm(dx,N)
    
    Pinv = lil_matrix((2*N+3,2*N+3))
    Pinv[0:N+2,0:N+2] = Pminusinv
    Pinv[N+2:2*N+3,N+2:2*N+3] = Pplusinv
    Pinv = Pinv.tocsr()
    return Pinv





























