import numpy as np
import cvxpy as cvx
import sbm_generator as gen


def sdp_relaxation(Y, z):
    print('Solving sdp relaxation problem...')
    n, _ = Y.shape
    X = cvx.Semidef(n)
    A = X * Y
    objective = cvx.Maximize(cvx.trace(A))
    constraints = [cvx.diag(X) == 1]
    problem = cvx.Problem(objective, constraints)
    problem.solve()
    print('Status: ' + problem.status)
    print('Optimal value: \n', problem.value)
    print('Verifying optimality (dual value): \n',
          np.sum(constraints[0].dual_value))
    print('Optimal X: \n', X.value)
    print('Optimal dual D (only diagonal entries): \n',
          constraints[0].dual_value)
    return problem.value, X, constraints[0].dual_value


if __name__ == "__main__":
    Y, z = gen.sbm_linear(10, 9, 2)
    sdp_relaxation(Y, z)
