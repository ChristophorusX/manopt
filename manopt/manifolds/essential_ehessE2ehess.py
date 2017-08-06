# Autogenerated with SMOP 
from smop.core import *
# essential/essential_ehessE2ehess.m

    
@function
def essential_ehessE2ehess(X=None,egradE=None,ehessE=None,S=None,*args,**kwargs):
    varargin = essential_ehessE2ehess.varargin
    nargin = essential_ehessE2ehess.nargin

    # Converts the Hessian in essential matrix E to the Hessian in X.
    
    # function ehess = essential_ehessE2ehess(X, egradE, ehessE, S)
    
    # egradE is the function handle for the gradient in E.
# ehessE is the function handle for the Hessian in E.
# S is the search direction in the space of X.
    
    # The output is a matrix in the space of X.
    
    # See also: essential_costE2cost essential_egradE2egrad
    
    # This file is part of Manopt: www.manopt.org.
# Original author: Roberto Tron, Aug. 8, 2014
# Contributors: Bamdev Mishra, May 22, 2015.
    
    e3hat=matlabarray(cat([0,- 1,0],[1,0,0],[0,0,0]))
# essential/essential_ehessE2ehess.m:19
    RA=X[:,1:3,:]
# essential/essential_ehessE2ehess.m:21
    RB=X[:,4:6,:]
# essential/essential_ehessE2ehess.m:22
    E=multiprod(multiprod(multitransp(RA),e3hat),RB)
# essential/essential_ehessE2ehess.m:23
    
    G=egradE[E]
# essential/essential_ehessE2ehess.m:24
    V=essential_sharp(multiprod(essential_flat(X),essential_flat(S)))
# essential/essential_ehessE2ehess.m:26
    VA=V[:,1:3,:]
# essential/essential_ehessE2ehess.m:27
    VB=V[:,4:6,:]
# essential/essential_ehessE2ehess.m:28
    dE=multiprod(multiprod(multitransp(RA),e3hat),VB) + multiprod(multiprod(multitransp(VA),e3hat),RB)
# essential/essential_ehessE2ehess.m:30
    dG=ehessE[E,dE]
# essential/essential_ehessE2ehess.m:32
    
    ehess=multiprod(e3hat,cat(2,multiprod(VB,multitransp(G)) + multiprod(RB,multitransp(dG)),- multiprod(VA,G) - multiprod(RA,dG)))
# essential/essential_ehessE2ehess.m:35
    return ehess
    
if __name__ == '__main__':
    pass
    