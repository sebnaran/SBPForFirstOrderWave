from scipy.sparse import csr_matrix
from scipy.sparse import lil_matrix
import numpy as np
import math
##These functions relate to the initiation of the code and mesh
##building.

def StaggeredSpatialMesh(N):
    #This function creates a pair of staggered mesh of [-1,1]
    #It takes the number of subpartitions and returns the 
    #two grids alongside with the mesh size dx.
    dx = 2/N
    xplus = [-1+dx*i for i in range(N+1)]
    xminus = [-1+dx/2+dx*i for i in range(N)]
    xminus = [-1] + xminus + [1]
    return xplus,xminus,dx

def ConstructDp(dx,N):
    #Here we will build the finite difference operator from xplus    #grid to the xminus grid.
    Dplus = lil_matrix((N+2,N+1))
    Dplus[0,0] = -1
    Dplus[0,1] = 1
    Dplus[N+1,N] = 1
    Dplus[N+1,N-1] = -1

    for row in range(1,N+1):
        Dplus[row,row-1] = -1
        Dplus[row,row] = 1
    return (1/dx)*Dplus
    
def ConstructDm(dx,N):
    #Here we will build the finite difference operator from 
    #xminus to the xplus grid.
    Dminus = lil_matrix((N+1,N+2))

    for row in range(N+1):
        Dminus[row,row] = -1
        Dminus[row,row+1] = 1

    return (1/dx)*Dminus

def ConstructAtilde(N,dx):
    #In this function we combine the two matrices Dplus and Dminus
    Dminus  = ConstructDm(dx,N)
    Dplus = ConstructDp(dx,N)
    A = lil_matrix((2*N+3,2*N+3))
    A[0:N+1, N+1:2*N+3] = Dminus
    A[N+1:2*N+3, 0:N+1] = Dplus
    A = A.tocsr()  #Here we change the data structure to one thatis suited for efficient computation.
    return A


def ConstructPpminv(dx,N):
    #Here we build the inner product matrix on the xplus grid.
    Pplusinv = lil_matrix((N+1,N+1))
    Pplusinv.setdiag(1)
    Pplusinv[0,0] = 2
    Pplusinv[N+1,N+1] = 2
    #And the inner product matrix on the xminus grid.
    Pminusinv = lil_matrix((N+2,N+2))
    Pminusinv.setdiag(1)
    Pminusinv[0,0] = 2
    Pminusinv[1,1] = 4
    Pminusinv[2,2] = 0.8
    Pminusinv[N+1,N+1] = 0.8
    Pminusinv[N,N] = 4
    Pminusinv[N-1,N-1] = 2
    return (1/dx)*Pplusinv, (1/dx)*Pminusinv

def ConstructPinv(dx,N):
    #Here we build the matrix P
    Pplusinv,Pminusinv = ConstructPpm(dx,N)
    
    Pinv = lil_matrix((2*N+3,2*N+3))
    Pinv[0:N+2,0:N+2] = Pminusinv
    Pinv[N+2:2*N+3,N+2:2*N+3] = Pplusinv
    Pinv = Pinv.tocsr()
    return Pinv





























