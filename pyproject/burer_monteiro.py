import numpy as np
import matplotlib.pyplot as plt
from scipy import optimize as opt
import sbm_generator as gensbm
import sync_generator as gensync
import aux
import time


def augmented_lagrangian(Y, k):
    start_pre = time.time()
    n, _ = Y.shape
    y = np.ones(n).reshape((-1, 1))
    R = np.random.random_sample((n, k))
    penalty = 1
    gamma = 10
    eta = .25
    target = .01
    vec = _constraint_term_vec(n, R)
    v = vec.reshape((1, -1)).dot(vec)
    v_best = v
    end_pre = time.time()
    start_loop = time.time()
    while v > target:
        Rv = _R_to_Rv(R, k)
        # print(_jacobian(Rv, Y, A, y, penalty, k))
        print('Starting L-BFGS-B on augmented Lagrangian...')
        optimizer = opt.minimize(lambda R_vec: _augmented_lagrangian_func(
            R_vec, y, penalty, n, k), Rv, jac=lambda R_vec: _jacobian(R_vec, Y, n, y, penalty, k), method="L-BFGS-B")
        print('Finishing L-BFGS-B on augmented Lagrangian...')
        R = _Rv_to_R(optimizer.x, n, k)
        vec = _constraint_term_vec(n, R)
        v = vec.reshape((1, -1)).dot(vec)
        print('Finish updating variables...')
        # print(R)
        _plot_R(R)
        # print('v: ', v)
        # print('v_best: ', v_best)
        # print(optimizer)
        if v < eta * v_best:
            y = y - penalty * vec
            v_best = v
        else:
            penalty = gamma * penalty
    end_loop = time.time()
    print('Timing pre part: ', end_pre - start_pre)
    print('Timing loop part: ', end_loop - start_loop)
    print('Augmented Lagrangian terminated.')


def _basis_vector(size, index):
    vec = np.zeros(size)
    vec[index] = 1
    return vec


def _A_trace_vec(n, R):
    # X = R.dot(R.T)
    vec = np.empty(n)
    for i in range(n):
        # vec[i] = np.trace(A[i].dot(X))
        vec[i] = R[i, :].dot(R[i, :])
    return vec.reshape((-1, 1))


def _constraint_term_vec(n, R):
    vec = _A_trace_vec(n, R)
    # for general case, all one vector should be replaced by a vector b
    constraint = vec - np.ones(n).reshape((-1, 1))
    return constraint


def _augmented_lagrangian_func(Rv, y, penalty, n, k):
    # print('Start computing objective function...')
    R = _Rv_to_R(Rv, n, k)
    # X = R.dot(R.T)
    vec = _constraint_term_vec(n, R)
    # print('-np.trace(Y.dot(X)): ', -np.trace(Y.dot(X)))
    # print('- y.reshape((1, -1)).dot(vec): ', - y.reshape((1, -1)).dot(vec))
    # print('vec.reshape((1, -1)).dot(vec): ', vec.reshape((1, -1)).dot(vec))
    # print('penalty / 2 * vec.reshape((1, -1)).dot(vec): ', penalty / 2 * vec.reshape((1, -1)).dot(vec))
    objective = -np.trace(Y.dot(R.dot(R.T))) - y.reshape((1, -1)
                                                         ).dot(vec) + penalty / 2 * vec.reshape((1, -1)).dot(vec)
    # print('Finish computing objective function!!!')
    return objective


def _Rv_to_R(Rv, n, k):
    U = Rv.reshape((n, k))
    # Rv = Rv.reshape((-1, 1))
    # U = Rv[0:n, 0].reshape((-1, 1))
    # for i in range(1, k):
    #     vec = Rv[i * n:(i + 1) * n, 0].reshape((-1, 1))
    #     U = np.hstack((U, vec))
    return U


def _R_to_Rv(R, k):
    u = R.reshape((1, -1)).ravel()
    # u = R[:, 0].reshape((-1, 1))
    # for i in range(1, k):
    #     u = np.vstack((u, R[:, i].reshape((-1, 1))))
    return u

# def _take_one_row(Z, R, l):
#     Z[l,:] = R[l,:]
#     return Z


def _jacobian(Rv, Y, n, y, penalty, k):
    # print('Start computing Jacobian...')
    R = _Rv_to_R(Rv, n, k)
    # X = R.dot(R.T)
    # y_row = y.reshape((1, -1)).ravel()
    # vec_trace_A = np.array([np.trace(A[l].dot(X))
    # for l in range(m)])
    vec_trace_A = _A_trace_vec(n, R).ravel()
    # vec_second_part = [y.reshape((1, -1)).ravel()
    #                    [l] * _take_one_row(np.zeros((n, k)), R, l) for l in range(n)]
    # vec_third_part = [(vec_trace_A[l] - 1) * _take_one_row(np.zeros((n, k)), R, l)
    #                   for l in range(n)]
    vec_second_part = R.copy()
    for l in range(n):
        vec_second_part[l,:] *= y.ravel()[l]
    vec_third_part = R.copy()
    for l in range(n):
        vec_third_part[l,:] *= (vec_trace_A[l] - 1)
    jacobian = -2 * Y.dot(R) - 2 * vec_second_part + 2 * penalty * vec_third_part
    jac_vec = _R_to_Rv(jacobian, k)
    # print('Finish computing Jacobian!!!')
    return jac_vec.reshape((1, -1)).ravel()

def _plot_R(R):
    plt.scatter(R.T[0], R.T[1], alpha=0.25, label='Rows of R')
    circle = plt.Circle((0, 0), 1, fill=False, linestyle='--', color='xkcd:blue violet')
    plt.gcf().gca().add_artist(circle)
    plt.legend()
    plt.show()


if __name__ == "__main__":
    Y, z = gensbm.sbm_linear(500, 10, 2)
    Y = aux.demean(Y, 10, 2)
    # Y, z = gensync.synchronization_usual(1000, .5, 2)
    # print(Y)
    # print(Y)
    # Y_re = Y.reshape((1, -1)).ravel()
    # print(Y_re)
    # Y_back = Y_re.reshape((10,10))
    # print(Y_back)
    augmented_lagrangian(Y, 2)
