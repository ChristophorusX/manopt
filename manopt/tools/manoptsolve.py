# Autogenerated with SMOP 
from smop.core import *
# manoptsolve.m

    
@function
def manoptsolve(problem=None,x0=None,options=None,*args,**kwargs):
    varargin = manoptsolve.varargin
    nargin = manoptsolve.nargin

    # Gateway helper function to call a Manopt solver, chosen in the options.
    
    # function [x, cost, info, options] = manoptsolve(problem)
# function [x, cost, info, options] = manoptsolve(problem, x0)
# function [x, cost, info, options] = manoptsolve(problem, x0, options)
# function [x, cost, info, options] = manoptsolve(problem, [], options)
    
    # Depending on what is available in the Manopt problem structure, one of
# the Manopt solvers will be called and the outputs passed along. It is
# also possible to force the choice of a solver by specifying it in the
# options structure. For example:
    
    #    options.solver = @trustregions;
    
    # Simply specify a function handle to a Manopt solver.
    
    # See also: trustregions conjugategradient steepestdescent
    
    # This file is part of Manopt: www.manopt.org.
# Original author: Nicolas Boumal, Aug. 13, 2014.
# Contributors: 
# Change log:
    
    # At the very least, we need a cost function.
    if logical_not(canGetCost(problem)):
        error('The problem structure must specify a cost function.')
    
    
    # Depending on the number of differentials available, pick a different
    # default solver.
    if logical_not(canGetGradient(problem)):
        localdefaults.solver = copy(neldermead)
# manoptsolve.m:33
    else:
        if logical_not(canGetHessian(problem)):
            localdefaults.solver = copy(conjugategradient)
# manoptsolve.m:35
        else:
            localdefaults.solver = copy(trustregions)
# manoptsolve.m:37
    
    
    # Merge local defaults with user options, if any.
    if logical_not(exist('options','var')) or isempty(options):
        options=struct()
# manoptsolve.m:42
    
    options=mergeOptions(localdefaults,options)
# manoptsolve.m:44
    
    if logical_not(exist('x0','var')):
        x0=matlabarray([])
# manoptsolve.m:48
    
    
    # Issue the actual call.
    x,cost,info,options=options.solver(problem,x0,options,nargout=4)
# manoptsolve.m:52
    return x,cost,info,options
    
if __name__ == '__main__':
    pass
    