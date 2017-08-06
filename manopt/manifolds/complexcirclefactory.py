# Autogenerated with SMOP 
from smop.core import *
# complexcircle/complexcirclefactory.m

    
@function
def complexcirclefactory(n=None,*args,**kwargs):
    varargin = complexcirclefactory.varargin
    nargin = complexcirclefactory.nargin

    # Returns a manifold struct to optimize over unit-modulus complex numbers.
    
    # function M = complexcirclefactory()
# function M = complexcirclefactory(n)
    
    # Description of vectors z in C^n (complex) such that each component z(i)
# has unit modulus. The manifold structure is the Riemannian submanifold
# structure from the embedding space R^2 x ... x R^2, i.e., the complex
# circle is identified with the unit circle in the real plane.
    
    # By default, n = 1.
    
    # See also spherecomplexfactory
    
    # This file is part of Manopt: www.manopt.org.
# Original author: Nicolas Boumal, Dec. 30, 2012.
# Contributors: 
# Change log:
    
    #   July 7, 2014 (NB): Added ehess2rhess function.
    
    #   Sep. 3, 2014 (NB): Correction to the dist function (extract real part).
    
    #   April 13, 2015 (NB): Fixed logarithm.
    
    #   Oct. 8, 2016 (NB)
#       Code for exponential was simplified to only treat the zero vector
#       as a particular case.
    
    if logical_not(exist('n','var')):
        n=1
# complexcircle/complexcirclefactory.m:32
    
    M.name = copy(lambda : sprintf('Complex circle (S^1)^%d',n))
# complexcircle/complexcirclefactory.m:35
    M.dim = copy(lambda : n)
# complexcircle/complexcirclefactory.m:37
    M.inner = copy(lambda z=None,v=None,w=None: real(dot(v.T,w)))
# complexcircle/complexcirclefactory.m:39
    M.norm = copy(lambda x=None,v=None: norm(v))
# complexcircle/complexcirclefactory.m:41
    M.dist = copy(lambda x=None,y=None: norm(acos(real(multiply(conj(x),y)))))
# complexcircle/complexcirclefactory.m:43
    M.typicaldist = copy(lambda : dot(pi,sqrt(n)))
# complexcircle/complexcirclefactory.m:45
    M.proj = copy(lambda z=None,u=None: u - multiply(real(multiply(conj(u),z)),z))
# complexcircle/complexcirclefactory.m:47
    M.tangent = copy(M.proj)
# complexcircle/complexcirclefactory.m:49
    
    # Riemannian gradient amounts to an orthogonal projection.
    M.egrad2rgrad = copy(M.proj)
# complexcircle/complexcirclefactory.m:53
    M.ehess2rhess = copy(ehess2rhess)
# complexcircle/complexcirclefactory.m:55
    
@function
def ehess2rhess(z=None,egrad=None,ehess=None,zdot=None,*args,**kwargs):
    varargin = ehess2rhess.varargin
    nargin = ehess2rhess.nargin

    rhess=M.proj(z,ehess - multiply(real(multiply(z,conj(egrad))),zdot))
# complexcircle/complexcirclefactory.m:57
    return rhess
    
if __name__ == '__main__':
    pass
    
    
    M.exp = copy(exponential)
# complexcircle/complexcirclefactory.m:60
    
@function
def exponential(z=None,v=None,t=None,*args,**kwargs):
    varargin = exponential.varargin
    nargin = exponential.nargin

    
    if nargin == 2:
        # t = 1;
        tv=copy(v)
# complexcircle/complexcirclefactory.m:65
    else:
        tv=dot(t,v)
# complexcircle/complexcirclefactory.m:67
    
    y=zeros(n,1)
# complexcircle/complexcirclefactory.m:70
    nrm_tv=abs(tv)
# complexcircle/complexcirclefactory.m:72
    
    mask=(nrm_tv > 0)
# complexcircle/complexcirclefactory.m:75
    y[mask]=multiply(z[mask],cos(nrm_tv[mask])) + multiply(tv[mask],(sin(nrm_tv[mask]) / nrm_tv[mask]))
# complexcircle/complexcirclefactory.m:76
    y[logical_not(mask)]=z[logical_not(mask)]
# complexcircle/complexcirclefactory.m:78
    return y
    
if __name__ == '__main__':
    pass
    
    
    M.retr = copy(retraction)
# complexcircle/complexcirclefactory.m:82
    
@function
def retraction(z=None,v=None,t=None,*args,**kwargs):
    varargin = retraction.varargin
    nargin = retraction.nargin

    if nargin == 2:
        # t = 1;
        tv=copy(v)
# complexcircle/complexcirclefactory.m:86
    else:
        tv=dot(t,v)
# complexcircle/complexcirclefactory.m:88
    
    y=z + tv
# complexcircle/complexcirclefactory.m:90
    y=y / abs(y)
# complexcircle/complexcirclefactory.m:91
    return y
    
if __name__ == '__main__':
    pass
    
    M.log = copy(logarithm)
# complexcircle/complexcirclefactory.m:94
    
@function
def logarithm(x1=None,x2=None,*args,**kwargs):
    varargin = logarithm.varargin
    nargin = logarithm.nargin

    v=M.proj(x1,x2 - x1)
# complexcircle/complexcirclefactory.m:96
    di=real(acos(real(multiply(conj(x1),x2))))
# complexcircle/complexcirclefactory.m:97
    nv=abs(v)
# complexcircle/complexcirclefactory.m:98
    factors=di / nv
# complexcircle/complexcirclefactory.m:99
    factors[di <= 1e-06]=1
# complexcircle/complexcirclefactory.m:100
    v=multiply(v,factors)
# complexcircle/complexcirclefactory.m:101
    return v
    
if __name__ == '__main__':
    pass
    
    
    M.hash = copy(lambda z=None: cat('z',hashmd5(cat([real(ravel(z))],[imag(ravel(z))]))))
# complexcircle/complexcirclefactory.m:104
    M.rand = copy(random)
# complexcircle/complexcirclefactory.m:106
    
@function
def random(*args,**kwargs):
    varargin = random.varargin
    nargin = random.nargin

    z=randn(n,1) + dot(1j,randn(n,1))
# complexcircle/complexcirclefactory.m:108
    z=z / abs(z)
# complexcircle/complexcirclefactory.m:109
    return z
    
if __name__ == '__main__':
    pass
    
    
    M.randvec = copy(randomvec)
# complexcircle/complexcirclefactory.m:112
    
@function
def randomvec(z=None,*args,**kwargs):
    varargin = randomvec.varargin
    nargin = randomvec.nargin

    # i*z(k) is a basis vector of the tangent vector to the k-th circle
    v=multiply(randn(n,1),(dot(1j,z)))
# complexcircle/complexcirclefactory.m:115
    v=v / norm(v)
# complexcircle/complexcirclefactory.m:116
    return v
    
if __name__ == '__main__':
    pass
    
    
    M.lincomb = copy(matrixlincomb)
# complexcircle/complexcirclefactory.m:119
    M.zerovec = copy(lambda x=None: zeros(n,1))
# complexcircle/complexcirclefactory.m:121
    M.transp = copy(lambda x1=None,x2=None,d=None: M.proj(x2,d))
# complexcircle/complexcirclefactory.m:123
    M.pairmean = copy(pairmean)
# complexcircle/complexcirclefactory.m:125
    
@function
def pairmean(z1=None,z2=None,*args,**kwargs):
    varargin = pairmean.varargin
    nargin = pairmean.nargin

    z=z1 + z2
# complexcircle/complexcirclefactory.m:127
    z=z / abs(z)
# complexcircle/complexcirclefactory.m:128
    return z
    
if __name__ == '__main__':
    pass
    
    M.vec = copy(lambda x=None,u_mat=None: cat([real(u_mat)],[imag(u_mat)]))
# complexcircle/complexcirclefactory.m:131
    M.mat = copy(lambda x=None,u_vec=None: u_vec[1:n] + dot(1j,u_vec[(n + 1):end()]))
# complexcircle/complexcirclefactory.m:132
    M.vecmatareisometries = copy(lambda : true)
# complexcircle/complexcirclefactory.m:133
    return z
    
if __name__ == '__main__':
    pass
    