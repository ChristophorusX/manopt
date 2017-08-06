# Autogenerated with SMOP 
from smop.core import *
# dominant_invariant_subspace_complex.m

    
@function
def dominant_invariant_subspace_complex(A=None,p=None,*args,**kwargs):
    varargin = dominant_invariant_subspace_complex.varargin
    nargin = dominant_invariant_subspace_complex.nargin

    # Returns a unitary basis of the dominant invariant p-subspace of A.
    
    # function X = dominant_invariant_subspace(A, p)
    
    # Input: A complex, Hermitian matrix A of size nxn and an integer p < n.
# Output: A complex, unitary matrix X of size nxp such that trace(X'*A*X)
#         is maximized. That is, the columns of X form a unitary basis
#         of a dominant subspace of dimension p of A.
    
    # The optimization is performed on the complex Grassmann manifold, since
# only the space spanned by the columns of X matters.
    
    # See dominant_invariant_subspace for more details in the real case.
    
    # See also: dominant_invariant_subspace grassmanncomplexfactory
    
    # This file is part of Manopt and is copyrighted. See the license file.
    
    # Main author: Nicolas Boumal, June 30, 2015
# Contributors:
    
    # Change log:
    
    # Generate some random data to test the function
    if logical_not(exist('A','var')) or isempty(A):
        A=randn(128) + dot(1j,randn(128))
# dominant_invariant_subspace_complex.m:27
        A=(A + A.T) / 2
# dominant_invariant_subspace_complex.m:28
    
    if logical_not(exist('p','var')) or isempty(p):
        p=3
# dominant_invariant_subspace_complex.m:31
    
    
    # Make sure the input matrix is Hermitian
    n=size(A,1)
# dominant_invariant_subspace_complex.m:35
    assert_(size(A,2) == n,'A must be square.')
    assert_(norm(A - A.T,'fro') < dot(n,eps),'A must be Hermitian.')
    assert_(p <= n,'p must be smaller than n.')
    
    Gr=grassmanncomplexfactory(n,p)
# dominant_invariant_subspace_complex.m:41
    problem.M = copy(Gr)
# dominant_invariant_subspace_complex.m:42
    problem.cost = copy(lambda X=None: - real(trace(dot(dot(X.T,A),X))))
# dominant_invariant_subspace_complex.m:43
    problem.egrad = copy(lambda X=None: dot(dot(- 2,A),X))
# dominant_invariant_subspace_complex.m:44
    problem.ehess = copy(lambda X=None,H=None: dot(dot(- 2,A),H))
# dominant_invariant_subspace_complex.m:45
    
    # These can be commented out.
    # checkgradient(problem);
    # pause;
    # checkhessian(problem);
    # pause;
    
    # Issue a call to a solver. A random initial guess will be chosen and
    # default options are selected except for the ones we specify here.
    options.Delta_bar = copy(dot(8,sqrt(p)))
# dominant_invariant_subspace_complex.m:56
    X,costX,info,options=trustregions(problem,[],options,nargout=4)
# dominant_invariant_subspace_complex.m:57
    
    
    fprintf('Options used:\\n')
    disp(options)
    
    # Riemannian Hessian on the tangent space at (any) X. Computing the
    # spectrum at the solution gives us some idea of the conditioning of
    # the problem. If we were to implement a preconditioner for the
    # Hessian, this would also inform us on its performance.
    
    # Notice that (typically) all eigenvalues of the Hessian at the
    # solution are positive, i.e., we find an isolated minimizer. If we
    # replace the Grassmann manifold by the Stiefel manifold, hence still
    # optimizing over orthonormal matrices but ignoring the invariance
    # cost(XQ) = cost(X) for all Q orthogonal, then we see
    # dim O(p) = p(p-1)/2 zero eigenvalues in the Hessian spectrum, making
    # the optimizer not isolated anymore.
    if Gr.dim() < 512:
        evs=hessianspectrum(problem,X)
# dominant_invariant_subspace_complex.m:76
        stairs(sort(evs))
        title(cat('Eigenvalues of the Hessian of the cost function ','at the solution'))
        xlabel('Eigenvalue number (sorted)')
        ylabel('Value of the eigenvalue')
    
    return X,info
    
if __name__ == '__main__':
    pass
    