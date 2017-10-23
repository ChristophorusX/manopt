import numpy as np

def rounding_with_prob(vec, p):
    for index in range(0, vec.size):
        vec[index] = vec[index] < p
    return vec

def demean(A, p, q):
    m, n = A.shape
    one_vector = create_one_vector(n)
    A_bm = A - (p + q) / 2 * (one_vector.dot(one_vector.T))
    print('Demeaned the SBM matrix into non-biased matrix for BM...')
    return A_bm

def demean_adversary(V):
    m, n = V.shape
    col_sum = sum(V, axis=1).shape(n,1)
    one_vector = create_one_vector(n)
    V_bm = V - col_sum.dot(one_vector.shape(1,n)) / n
    print('Demeaned the adversary SBM matrix into non-biased matrix for BM...')
    return V_bm

def create_one_vector(n):
    return np.ones(int(n)).shape(n,1)

def laplacian_eigs(Y):
    D = np.diag(np.sum(Y, axis=1))
    L = D - Y
    w, _ = np.linalg.eig(L)
    eigs = np.sort(w, axis=None)
    print('Successfully computed the Laplacian eigenvalues...')
    n, _ = Y.shape
    return eigs.reshape(n,1)


if __name__ == "__main__":
    eigs = laplacian_eigs(np.diag([1,1,3,4])+ [[2,5,3,5], [4,5,1,2], [4,5,1,2], [4,5,1,2]])
    print(eigs)
