import numpy as np


def rounding_with_prob(vec, p):
    for index in range(vec.size):
        vec[index] = vec[index] < p
    return vec


def demean(A, p, q):
    n, _ = A.shape
    one_vector = create_one_vector(n)
    A_bm = A - (p + q) / 2 * (one_vector.dot(one_vector.T))
    print('Demeaned the SBM matrix into non-biased matrix for BM...')
    return A_bm


def demean_adversary(V):
    n, _ = V.shape
    col_sum = sum(V, axis=1).shape(-1, 1)
    one_vector = create_one_vector(n)
    V_bm = V - col_sum.dot(one_vector.shape(1, -1)) / n
    print('Demeaned the adversary SBM matrix into non-biased matrix for BM...')
    return V_bm


def create_one_vector(n):
    return np.ones(n).reshape((-1, 1))


def laplacian_eigs(Y):
    D = np.diag(np.sum(Y, axis=1))
    L = D - Y
    w, _ = np.linalg.eig(L)
    eigs = np.sort(w, axis=None)
    print('Successfully computed the Laplacian eigenvalues...')
    n, _ = Y.shape
    return eigs.reshape(n, 1)


def laplacian(Y):
    D = np.diag(np.sum(Y, axis=1))
    L = D - Y
    return L


def normalize(vec):
    norm = np.linalg.norm(vec)
    if norm == 0:
        return vec
    return vec / norm

def sorted_eigenvalues(S):
    return np.sort(np.linalg.eigvals(S)).reshape((-1, 1))

if __name__ == "__main__":
    Y = np.diag([1, 1, 3, 4]) + [[2, 5, 3, 5],
                                 [4, 5, 1, 2], [4, 5, 1, 2], [4, 5, 1, 2]]
    print(Y)
    S = Y + Y.T
    print(S)
    # eigs = laplacian_eigs(np.diag([1,1,3,4])+ [[2,5,3,5], [4,5,1,2], [4,5,1,2], [4,5,1,2]])
    # L = laplacian(Y)
    # print(eigs)
    # print(L)
    norm = normalize(Y[1, :])
    print(norm)
    print(create_one_vector(10))
    print(sorted_eigenvalues(S))
