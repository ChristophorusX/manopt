# Autogenerated with SMOP 
from smop.core import *
# tangentorthobasis.m

    
@function
def tangentorthobasis(M=None,x=None,n=None,*args,**kwargs):
    varargin = tangentorthobasis.varargin
    nargin = tangentorthobasis.nargin

    # Returns an orthonormal basis of tangent vectors in the Manopt framework.
    
    # function orthobasis = tangentorthobasis(M, x)
# function orthobasis = tangentorthobasis(M, x, n)
    
    # M is a Manopt manifold structure obtained from a factory.
# x is a point on the manifold M.
# n (optional) is the dimension of the random subspace to span; by default,
#   n = M.dim() so that the returned basis spans the whole tangent space.
    
    # orthobasis is a cell of n tangent vectors at x.
# With high probability, they form an orthonormal basis of the tangent
# space at x. If necessary, this can be checked by calling
#   G = grammatrix(orthobasis)
# and verifying that G is invertible.
    
    # See also: grammatrix orthogonalize lincomb plotprofile
    
    # This file is part of Manopt: www.manopt.org.
# Original author: Nicolas Boumal, April 28, 2016.
# Contributors: 
# Change log:
    
    dim=M.dim()
# tangentorthobasis.m:25
    if logical_not(exist('n','var')) or isempty(n):
        n=copy(dim)
# tangentorthobasis.m:27
    
    assert_(n >= 0 and n <= dim and n == round(n),'n must be an integer between 0 and M.dim().')
    basis=cell(n,1)
# tangentorthobasis.m:32
    
    # are linearly independent.
    for k in arange(1,n).reshape(-1):
        basis[k]=M.randvec(x)
# tangentorthobasis.m:37
    
    
    # The Gram-Schmidt process transforms any n linearly independent
    # vectors into n orthonormal vectors spanning the same subspace.
    orthobasis=orthogonalize(M,x,basis)
# tangentorthobasis.m:42
    return orthobasis
    
if __name__ == '__main__':
    pass
    