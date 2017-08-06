# Autogenerated with SMOP 
from smop.core import *
# steepestdescent/steepestdescent.m

    
@function
def steepestdescent(problem=None,x=None,options=None,*args,**kwargs):
    varargin = steepestdescent.varargin
    nargin = steepestdescent.nargin

    # Steepest descent (gradient descent) minimization algorithm for Manopt.
    
    # function [x, cost, info, options] = steepestdescent(problem)
# function [x, cost, info, options] = steepestdescent(problem, x0)
# function [x, cost, info, options] = steepestdescent(problem, x0, options)
# function [x, cost, info, options] = steepestdescent(problem, [], options)
    
    # Apply the steepest descent minimization algorithm to the problem defined
# in the problem structure, starting at x0 if it is provided (otherwise, at
# a random point on the manifold). To specify options whilst not specifying
# an initial guess, give x0 as [] (the empty matrix).
    
    # In most of the examples bundled with the toolbox (see link below), the
# solver can be replaced by the present one if need be.
    
    # The outputs x and cost are the best reached point on the manifold and its
# cost. The struct-array info contains information about the iterations:
#   iter : the iteration number (0 for the initial guess)
#   cost : cost value
#   time : elapsed time in seconds
#   gradnorm : Riemannian norm of the gradient
#   stepsize : norm of the last tangent vector retracted
#   linesearch : information logged by options.linesearch
#   And possibly additional information logged by options.statsfun.
# For example, type [info.gradnorm] to obtain a vector of the successive
# gradient norms reached.
    
    # The options structure is used to overwrite the default values. All
# options have a default value and are hence optional. To force an option
# value, pass an options structure with a field options.optionname, where
# optionname is one of the following and the default value is indicated
# between parentheses:
    
    #   tolgradnorm (1e-6)
#       The algorithm terminates if the norm of the gradient drops below this.
#   maxiter (1000)
#       The algorithm terminates if maxiter iterations have been executed.
#   maxtime (Inf)
#       The algorithm terminates if maxtime seconds elapsed.
#   minstepsize (1e-10)
#       The algorithm terminates if the linesearch returns a displacement
#       vector (to be retracted) smaller in norm than this value.
#   linesearch (@linesearch or @linesearch_hint)
#       Function handle to a line search function. The options structure is
#       passed to the line search too, so you can pass it parameters. See
#       each line search's documentation for info. Another available line
#       search in manopt is @linesearch_adaptive, in
#       /manopt/linesearch/linesearch_adaptive.m
#       If the problem structure includes a line search hint, then the
#       default line search used is @linesearch_hint.
#   statsfun (none)
#       Function handle to a function that will be called after each
#       iteration to provide the opportunity to log additional statistics.
#       They will be returned in the info struct. See the generic Manopt
#       documentation about solvers for further information.
#   stopfun (none)
#       Function handle to a function that will be called at each iteration
#       to provide the opportunity to specify additional stopping criteria.
#       See the generic Manopt documentation about solvers for further
#       information.
#   verbosity (3)
#       Integer number used to tune the amount of output the algorithm
#       generates during execution (mostly as text in the command window).
#       The higher, the more output. 0 means silent.
#   storedepth (2)
#       Maximum number of different points x of the manifold for which a
#       store structure will be kept in memory in the storedb. If the
#       caching features of Manopt are not used, this is irrelevant. For
#       the SD algorithm, a store depth of 2 should always be sufficient.
    
    
    # See also: conjugategradient trustregions manopt/solvers/linesearch manopt/examples
    
    # This file is part of Manopt: www.manopt.org.
# Original author: Nicolas Boumal, Dec. 30, 2012.
# Contributors: 
# Change log:
    
    #   April 3, 2015 (NB):
#       Works with the new StoreDB class system.
    
    
    # Verify that the problem description is sufficient for the solver.
    if logical_not(canGetCost(problem)):
        warning('manopt:getCost','No cost provided. The algorithm will likely abort.')
    
    if logical_not(canGetGradient(problem)) and logical_not(canGetApproxGradient(problem)):
        # Note: we do not give a warning if an approximate Hessian is
        # explicitly given in the problem description, as in that case the user
        # seems to be aware of the issue.
        warning('manopt:getGradient:approx',cat('No gradient provided. Using an FD approximation instead (slow).\\n','It may be necessary to increase options.tolgradnorm.\\n','To disable this warning: warning(\'off\', \'manopt:getGradient:approx\')'))
        problem.approxgrad = copy(approxgradientFD(problem))
# steepestdescent/steepestdescent.m:97
    
    
    # Set local defaults here
    localdefaults.minstepsize = copy(1e-10)
# steepestdescent/steepestdescent.m:101
    localdefaults.maxiter = copy(1000)
# steepestdescent/steepestdescent.m:102
    localdefaults.tolgradnorm = copy(1e-06)
# steepestdescent/steepestdescent.m:103
    
    # line-search algorithms, choose a default line-search that works on
    # its own (typical) or that uses the hint.
    if logical_not(canGetLinesearch(problem)):
        localdefaults.linesearch = copy(linesearch)
# steepestdescent/steepestdescent.m:109
    else:
        localdefaults.linesearch = copy(linesearch_hint)
# steepestdescent/steepestdescent.m:111
    
    
    # Merge global and local defaults, then merge w/ user options, if any.
    localdefaults=mergeOptions(getGlobalDefaults(),localdefaults)
# steepestdescent/steepestdescent.m:115
    if logical_not(exist('options','var')) or isempty(options):
        options=struct()
# steepestdescent/steepestdescent.m:117
    
    options=mergeOptions(localdefaults,options)
# steepestdescent/steepestdescent.m:119
    timetic=tic()
# steepestdescent/steepestdescent.m:121
    
    if logical_not(exist('x','var')) or isempty(x):
        x=problem.M.rand()
# steepestdescent/steepestdescent.m:125
    
    
    # Create a store database and get a key for the current x
    storedb=StoreDB(options.storedepth)
# steepestdescent/steepestdescent.m:129
    key=storedb.getNewKey()
# steepestdescent/steepestdescent.m:130
    
    cost,grad=getCostGrad(problem,x,storedb,key,nargout=2)
# steepestdescent/steepestdescent.m:133
    gradnorm=problem.M.norm(x,grad)
# steepestdescent/steepestdescent.m:134
    
    # At any point, iter is the number of fully executed iterations so far.
    iter=0
# steepestdescent/steepestdescent.m:138
    
    stats=savestats()
# steepestdescent/steepestdescent.m:141
    info[1]=stats
# steepestdescent/steepestdescent.m:142
    info[min(10000,options.maxiter + 1)].iter = copy([])
# steepestdescent/steepestdescent.m:143
    if options.verbosity >= 2:
        fprintf(' iter\\t               cost val\\t    grad. norm\\n')
    
    
    # Start iterating until stopping criterion triggers
    while true:

        # Display iteration information
        if options.verbosity >= 2:
            fprintf('%5d\\t%+.16e\\t%.8e\\n',iter,cost,gradnorm)
        # Start timing this iteration
        timetic=tic()
# steepestdescent/steepestdescent.m:158
        stop,reason=stoppingcriterion(problem,x,options,info,iter + 1,nargout=2)
# steepestdescent/steepestdescent.m:161
        if logical_not(stop) and stats.stepsize < options.minstepsize:
            stop=copy(true)
# steepestdescent/steepestdescent.m:166
            reason=sprintf(cat('Last stepsize smaller than minimum ','allowed; options.minstepsize = %g.'),options.minstepsize)
# steepestdescent/steepestdescent.m:167
        if stop:
            if options.verbosity >= 1:
                fprintf(cat(reason,'\\n'))
            break
        # Pick the descent direction as minus the gradient
        desc_dir=problem.M.lincomb(x,- 1,grad)
# steepestdescent/steepestdescent.m:180
        stepsize,newx,newkey,lsstats=options.linesearch(problem,x,desc_dir,cost,- gradnorm ** 2,options,storedb,key,nargout=4)
# steepestdescent/steepestdescent.m:183
        newcost,newgrad=getCostGrad(problem,newx,storedb,newkey,nargout=2)
# steepestdescent/steepestdescent.m:188
        newgradnorm=problem.M.norm(newx,newgrad)
# steepestdescent/steepestdescent.m:189
        storedb.purge()
        x=copy(newx)
# steepestdescent/steepestdescent.m:195
        key=copy(newkey)
# steepestdescent/steepestdescent.m:196
        cost=copy(newcost)
# steepestdescent/steepestdescent.m:197
        grad=copy(newgrad)
# steepestdescent/steepestdescent.m:198
        gradnorm=copy(newgradnorm)
# steepestdescent/steepestdescent.m:199
        iter=iter + 1
# steepestdescent/steepestdescent.m:202
        stats=savestats()
# steepestdescent/steepestdescent.m:205
        info[iter + 1]=stats
# steepestdescent/steepestdescent.m:206

    
    
    
    info=info[1:iter + 1]
# steepestdescent/steepestdescent.m:211
    if options.verbosity >= 1:
        fprintf('Total time is %f [s] (excludes statsfun)\\n',info[end()].time)
    
    
    
    
    # Routine in charge of collecting the current iteration stats
    
@function
def savestats(*args,**kwargs):
    varargin = savestats.varargin
    nargin = savestats.nargin

    stats.iter = copy(iter)
# steepestdescent/steepestdescent.m:222
    stats.cost = copy(cost)
# steepestdescent/steepestdescent.m:223
    stats.gradnorm = copy(gradnorm)
# steepestdescent/steepestdescent.m:224
    if iter == 0:
        stats.stepsize = copy(NaN)
# steepestdescent/steepestdescent.m:226
        stats.time = copy(toc(timetic))
# steepestdescent/steepestdescent.m:227
        stats.linesearch = copy([])
# steepestdescent/steepestdescent.m:228
    else:
        stats.stepsize = copy(stepsize)
# steepestdescent/steepestdescent.m:230
        stats.time = copy(info[iter].time + toc(timetic))
# steepestdescent/steepestdescent.m:231
        stats.linesearch = copy(lsstats)
# steepestdescent/steepestdescent.m:232
    
    stats=applyStatsfun(problem,x,storedb,key,options,stats)
# steepestdescent/steepestdescent.m:234
    return stats
    
if __name__ == '__main__':
    pass
    
    
    return stats
    
if __name__ == '__main__':
    pass
    