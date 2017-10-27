import numpy as np
import matplotlib.pyplot as plt
from scipy import optimize as opt
import sbm_generator as gensbm
import sync_generator as gensync
import aux


def augmented_lagrangian(Y, k):
    n, _ = Y.shape
    A = np.empty((n, n, n))
    for i in range(n):
        A[i] = np.diag(_basis_vector(n, i))
    y = np.ones(n).reshape((-1, 1))
    R = np.random.random_sample((n, k))
    penalty = 1
    gamma = 10
    eta = .25
    target = .01
    vec = _constraint_term_vec(A, R)
    v = vec.reshape((1, -1)).dot(vec)
    v_best = v
    counter = 0
    while v > target and counter < 10:
        counter += 1
        Rv = _R_to_Rv(R, k)
        # print(_jacobian(Rv, Y, A, y, penalty, k))
        print('Starting L-BFGS-B on augmented Lagrangian...')
        optimizer = opt.minimize(lambda R_vec: _augmented_lagrangian_func(
            A, R_vec, y, penalty, n, k), Rv, jac=lambda R_vec: _jacobian(R_vec, Y, A, y, penalty, k), method="L-BFGS-B")
        print('Finishing L-BFGS-B on augmented Lagrangian...')
        R = _Rv_to_R(optimizer.x.reshape((-1, 1)), n, k)
        X = R.dot(R.T)
        vec = _constraint_term_vec(A, X)
        v = vec.reshape((1, -1)).dot(vec)
        print('Finish updating variables...')
        # print(R)
        # plt.scatter(R.T[0], R.T[1], alpha=0.15, color='orange')
        # circle = plt.Circle((0, 0), 1, fill=False, linestyle='--')
        # plt.gcf().gca().add_artist(circle)
        # plt.show()
        # print('v: ', v)
        # print('v_best: ', v_best)
        print(_Rv_to_R(optimizer.x))
        if v < eta * v_best:
            y = y - penalty * vec
            v_best = v
        else:
            penalty = gamma * penalty
    print('Augmented Lagrangian terminated.')


def _basis_vector(size, index):
    vec = np.zeros(size)
    vec[index] = 1
    return vec


def _A_trace_vec(A, X):
    n, _, _ = A.shape
    vec = np.empty(n)
    for i in range(n):
        vec[i] = np.trace(A[i].dot(X))
    return vec.reshape((-1, 1))


def _constraint_term_vec(A, X):
    n, _, _ = A.shape
    vec = _A_trace_vec(A, X).reshape((1, -1))
    # for general case, all one vector should be replaced by a vector b
    constraint = vec - np.ones(n)
    return constraint.reshape((-1, 1))


def _augmented_lagrangian_func(A, Rv, y, penalty, n, k):
    Rv = Rv.reshape((-1, 1))
    R = _Rv_to_R(Rv, n, k)
    X = R.dot(R.T)
    vec = _constraint_term_vec(A, X)
    # print('-np.trace(Y.dot(X)): ', -np.trace(Y.dot(X)))
    # print('- y.reshape((1, -1)).dot(vec): ', - y.reshape((1, -1)).dot(vec))
    # print('vec.reshape((1, -1)).dot(vec): ', vec.reshape((1, -1)).dot(vec))
    # print('penalty / 2 * vec.reshape((1, -1)).dot(vec): ', penalty / 2 * vec.reshape((1, -1)).dot(vec))
    return -np.trace(Y.dot(X)) - y.reshape((1, -1)).dot(vec) + penalty / 2 * vec.reshape((1, -1)).dot(vec)


def _Rv_to_R(Rv, n, k):
    Rv = Rv.reshape((-1, 1))
    U = Rv[0:n, 0].reshape((-1, 1))
    for i in range(1, k):
        vec = Rv[i * n:(i + 1) * n, 0].reshape((-1, 1))
        U = np.hstack((U, vec))
    return U


def _R_to_Rv(R, k):
    u = R[:, 0].reshape((-1, 1))
    for i in range(1, k):
        u = np.vstack((u, R[:, i].reshape((-1, 1))))
    return u


def _jacobian(Rv, Y, A, y, penalty, k):
    m, n, _ = A.shape
    R = _Rv_to_R(Rv, n, k)
    X = R.dot(R.T)
    y_row = y.reshape((1, -1)).ravel()
    b_row = np.ones(n)
    vec_trace_A = np.array([np.trace(A[l].dot(X))
                            for l in range(m)])
    vec_second_part = [y_row[l] * A[l].dot(R) for l in range(m)]
    vec_third_part = [(vec_trace_A[l] - b_row[l]) * A[l].dot(R)
                      for l in range(m)]
    jacobian = -2 * Y.dot(R) - 2 * np.sum(vec_second_part,
                                          axis=0) + 2 * penalty * np.sum(vec_third_part, axis=0)
    jac_vec = _R_to_Rv(jacobian, k)
    return jac_vec.reshape((1, -1)).ravel()


if __name__ == "__main__":
    Y, z = gensbm.sbm_linear(500, 10, 2)
    Y = aux.demean(Y, 10, 2)
    # Y, z = gensync.synchronization_usual(50, .5, 20)
    # print(Y)
    augmented_lagrangian(Y, 2)
