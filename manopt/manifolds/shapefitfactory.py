# Autogenerated with SMOP 
from smop.core import *
# euclidean/shapefitfactory.m

    
@function
def shapefitfactory(VJt=None,*args,**kwargs):
    varargin = shapefitfactory.varargin
    nargin = shapefitfactory.nargin

    # Linear manifold structure for optimization over the ShapeFit search space
    
    # function M = shapefitfactory(VJt)
    
    # Input: VJt is a matrix of size dxn, such that VJt * ones(n, 1) = 0.
    
    # Returns M, a structure describing the Euclidean space of d-by-n matrices
# equipped with the standard Frobenius distance and associated trace inner
# product, as a manifold for Manopt. Matrices on M, denoted by T, have size
# dxn and obey T*ones(n, 1) = 0 (centered columns) and <VJt, T> = 1, where
# <A, B> = Trace(A' * B).
    
    # See this paper: http://arxiv.org/abs/1506.01437
# ShapeFit: Exact location recovery from corrupted pairwise directions, 2015
# Paul Hand, Choongbum Lee, Vladislav Voroninski
    
    # See also: shapefit_smoothed
    
    # This file is part of Manopt: www.manopt.org.
# Original author: Nicolas Boumal, June 18, 2015.
# Contributors: 
# Change log:
    
    d,n=size(VJt,nargout=2)
# euclidean/shapefitfactory.m:25
    M.name = copy(lambda : sprintf('ShapeFit space of size %d x %d',d,n))
# euclidean/shapefitfactory.m:27
    M.dim = copy(lambda : dot(d,n) - d - 1)
# euclidean/shapefitfactory.m:29
    M.inner = copy(lambda x=None,d1=None,d2=None: dot(ravel(d1).T,ravel(d2)))
# euclidean/shapefitfactory.m:31
    M.norm = copy(lambda x=None,d=None: norm(d,'fro'))
# euclidean/shapefitfactory.m:33
    M.dist = copy(lambda x=None,y=None: norm(x - y,'fro'))
# euclidean/shapefitfactory.m:35
    M.typicaldist = copy(lambda : sqrt(dot(d,n)))
# euclidean/shapefitfactory.m:37
    M.proj = copy(lambda T=None,U=None: projection(U))
# euclidean/shapefitfactory.m:39
    VJt_normed=VJt / norm(VJt,'fro')
# euclidean/shapefitfactory.m:40
    
@function
def projection(U=None,*args,**kwargs):
    varargin = projection.varargin
    nargin = projection.nargin

    # Center the columns
    PU=bsxfun(minus,U,mean(U,2))
# euclidean/shapefitfactory.m:43
    
    # Note: these two actions can be executed separately, without
        # interference, owing to VJt having centered columns itself.
    PU=PU - dot((dot(ravel(VJt_normed).T,ravel(U))),VJt_normed)
# euclidean/shapefitfactory.m:47
    return PU
    
if __name__ == '__main__':
    pass
    
    
    M.egrad2rgrad = copy(M.proj)
# euclidean/shapefitfactory.m:50
    M.ehess2rhess = copy(lambda x=None,eg=None,eh=None,d=None: projection(eh))
# euclidean/shapefitfactory.m:52
    M.tangent = copy(lambda x=None,d=None: d)
# euclidean/shapefitfactory.m:54
    M.exp = copy(exp)
# euclidean/shapefitfactory.m:56
    
@function
def exp(x=None,d=None,t=None,*args,**kwargs):
    varargin = exp.varargin
    nargin = exp.nargin

    if nargin == 3:
        y=x + dot(t,d)
# euclidean/shapefitfactory.m:59
    else:
        y=x + d
# euclidean/shapefitfactory.m:61
    
    return y
    
if __name__ == '__main__':
    pass
    
    
    M.retr = copy(M.exp)
# euclidean/shapefitfactory.m:65
    M.log = copy(lambda x=None,y=None: y - x)
# euclidean/shapefitfactory.m:67
    M.hash = copy(lambda x=None: cat('z',hashmd5(ravel(x))))
# euclidean/shapefitfactory.m:69
    M.randvec = copy(lambda x=None: randvec())
# euclidean/shapefitfactory.m:71
    
@function
def randvec(*args,**kwargs):
    varargin = randvec.varargin
    nargin = randvec.nargin

    u=projection(randn(d,n))
# euclidean/shapefitfactory.m:73
    u=u / norm(u,'fro')
# euclidean/shapefitfactory.m:74
    return u
    
if __name__ == '__main__':
    pass
    
    
    # We exploit the fact that VJt_normed belongs to the manifold
    M.rand = copy(lambda : VJt_normed + dot(randn(1),randvec()))
# euclidean/shapefitfactory.m:78
    M.lincomb = copy(matrixlincomb)
# euclidean/shapefitfactory.m:80
    M.zerovec = copy(lambda x=None: zeros(d,n))
# euclidean/shapefitfactory.m:82
    M.transp = copy(lambda x1=None,x2=None,d=None: d)
# euclidean/shapefitfactory.m:84
    M.pairmean = copy(lambda x1=None,x2=None: dot(0.5,(x1 + x2)))
# euclidean/shapefitfactory.m:86
    M.vec = copy(lambda x=None,u_mat=None: ravel(u_mat))
# euclidean/shapefitfactory.m:88
    M.mat = copy(lambda x=None,u_vec=None: reshape(u_vec,cat(d,n)))
# euclidean/shapefitfactory.m:89
    M.vecmatareisometries = copy(lambda : true)
# euclidean/shapefitfactory.m:90
    return u
    
if __name__ == '__main__':
    pass
    