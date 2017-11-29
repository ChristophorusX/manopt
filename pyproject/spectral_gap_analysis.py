import numpy as np


def spectral_gap(A):
    return np.sort(np.linalg.eigvals(A))[1]


if __name__ == '__main__':
    matrix = [[1, 2, 5], [3, 4, 6], [7, 8, 9]]
    matrix = np.array(matrix)
    print(spectral_gap(matrix))
