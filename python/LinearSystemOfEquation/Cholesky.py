from scipy import array, linalg, dot
# Cholesky 分解
import numpy as np
A = np.array([[1,2,4],[2,13,23],[4,23,77]])
L0 = linalg.cholesky(A, lower=True)
n = 3
L = np.zeros((n,n))# 修正
d = np.zeros((n))


def Cholesky():
    # A = {L}{L^T}
    v = np.zeros((n))
    for j in range(0,n):
        for i in range(j,n):
            v[i] = A[i,j]
            for k in range(0,j):
                v[i] -= L[j,k]*L[i,k]
        L[i,j] = v[i] / np.sqrt(v[j])

def Modified1():
    # A = {L}{d}{L^T}
    d[0] = A[0,0]
    L[0,0] = 1
    for i in range(1,n):
        
        for j in range(0,i):
            lld = A[i,j]
            for k in range(0,j):
                lld -= L[i,k]*L[j,k]*d[k]
            L[i,j] = 1 / d[j] * lld
            
        ld = A[i,i]
        for k in range(0,i):
            ld -= L[i,k]*L[i,k]*d[k]
        d[i] = ld
        L[i,i] = 1
        
def Modified2():
    # A = {L}{d}{L^T}
    L[0,0] = A[0,0]
    d[0] = 1 / L[0,0]
    for i in range(1,n):
        for j in range(0,i+1):
            lld = A[i,j]
            for k in range(0,j):
                lld -= L[i,k]*L[j,k]*d[k]
            L[i,j] = lld
        d[i] = 1 / L[i,i]
        
def IncompleteCholesky():
    d[0] = A[0,0]
    L[0,0] = 1
    for i in range(1,n):
        for j in range(0,i):
            if(abs(A[i,j]) < 1e-10):
                continue
            lld = A[i,j]
            for k in range(0,j):
                lld -= L[i,k]*L[j,k]*d[k]
            L[i,j] = 1 / d[j]
        
Modified2()
        

        