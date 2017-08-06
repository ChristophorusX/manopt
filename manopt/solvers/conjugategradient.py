# Autogenerated with SMOP 
from smop.core import *
# conjugategradient/conjugategradient.m

    
@function
def conjugategradient(problem=None,x=None,options=None,*args,**kwargs):
    varargin = conjugategradient.varargin
    nargin = conjugategradient.nargin

    # Conjugate gradient minimization algorithm for Manopt.
    
    # function [x, cost, info, options] = conjugategradient(problem)
# function [x, cost, info, options] = conjugategradient(problem, x0)
# function [x, cost, info, options] = conjugategradient(problem, x0, options)
# function [x, cost, info, options] = conjugategradient(problem, [], options)
    
    # Apply the conjugate gradient minimization algorithm to the problem
# defined in the problem structure, starting at x0 if it is provided
# (otherwise, at a random point on the manifold). To specify options whilst
# not specifying an initial guess, give x0 as [] (the empty matrix).
    
    # The outputs x and cost are the best reached point on the manifold and its
# cost. The struct-array info contains information about the iterations:
#   iter : the iteration number (0 for the initial guess)
#   cost : cost value
#   time : elapsed time in seconds
#   gradnorm : Riemannian norm of the gradient
#   stepsize : norm of the last tangent vector retracted
#   beta : value of the beta parameter (see options.beta_type)
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
#   beta_type ('H-S')
#       Conjugate gradient beta rule used to construct the new search
#       direction, based on a linear combination of the previous search
#       direction and the new (preconditioned) gradient. Possible values
#       for this parameter are:
#           'S-D', 'steep' for beta = 0 (preconditioned steepest descent)
#           'F-R' for Fletcher-Reeves's rule
#           'P-R' for Polak-Ribiere's modified rule
#           'H-S' for Hestenes-Stiefel's modified rule
#           'H-Z' for Hager-Zhang's modified rule
#       See Hager and Zhang 2006, "A survey of nonlinear conjugate gradient
#       methods" for a description of these rules in the Euclidean case and
#       for an explanation of how to adapt them to the preconditioned case.
#       The adaption to the Riemannian case is straightforward: see in code
#       for details. Modified rules take the max between 0 and the computed
#       beta value, which provides automatic restart, except for H-Z which
#       uses a different modification.
#   orth_value (Inf)
#       Following Powell's restart strategy (Math. prog. 1977), restart CG
#       (that is, make a -preconditioned- gradient step) if two successive
#       -preconditioned- gradients are "too" parallel. See for example
#       Hager and Zhang 2006, "A survey of nonlinear conjugate gradient
#       methods", page 12. An infinite value disables this strategy. See in
#       code formula for the specific criterion used.
#   linesearch (@linesearch_adaptive or @linesearch_hint)
#       Function handle to a line search function. The options structure is
#       passed to the line search too, so you can pass it parameters. See
#       each line search's documentation for info. Another available line
#       search in manopt is @linesearch, in /manopt/linesearch/linesearch.m
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
#       the CG algorithm, a store depth of 2 should always be sufficient.
    
    
    # In most of the examples bundled with the toolbox (see link below), the
# solver can be replaced by the present one if need be.
    
    # See also: steepestdescent trustregions manopt/solvers/linesearch manopt/examples
    
    # An explicit, general listing of this algorithm, with preconditioning,
# can be found in the following paper:
#     @Article{boumal2015lowrank,
#       Title   = {Low-rank matrix completion via preconditioned optimization on the {G}rassmann manifold},
#       Author  = {Boumal, N. and Absil, P.-A.},
#       Journal = {Linear Algebra and its Applications},
#       Year    = {2015},
#       Pages   = {200--239},
#       Volume  = {475},
#       Doi     = {10.1016/j.laa.2015.02.027},
#     }
    
    # This file is part of Manopt: www.manopt.org.
# Original author: Bamdev Mishra, Dec. 30, 2012.
# Contributors: Nicolas Boumal
# Change log:
    
    #   March 14, 2013, NB:
#       Added preconditioner support : see Section 8 in
#       https://www.math.lsu.edu/~hozhang/papers/cgsurvey.pdf
#    
#   Sept. 13, 2013, NB:
#       Now logging beta parameter too.
#    
#	Nov. 7, 2013, NB:
#       The search direction is no longer normalized before it is passed
#       to the linesearch. This way, it is up to the designers of the
#       linesearch to decide whether they want to use the norm of the
#       search direction in their algorithm or not. There are reasons
#       against it, but practical evidence that it may help too, so we
#       allow it. The default linesearch_adaptive used does exploit the
#       norm information. The base linesearch does not. You may select it
#       by setting options.linesearch = @linesearch;
    
    #	Nov. 29, 2013, NB:
#       Documentation improved: options are now explicitly described.
#       Removed the Daniel rule for beta: it was not appropriate for
#       preconditioned CG and I could not find a proper reference for it.
    
    #   April 3, 2015 (NB):
#       Works with the new StoreDB class system.
    
    # Verify that the problem description is sufficient for the solver.
    if logical_not(canGetCost(problem)):
        warning('manopt:getCost','No cost provided. The algorithm will likely abort.')
    
    if logical_not(canGetGradient(problem)) and logical_not(canGetApproxGradient(problem)):
        warning('manopt:getGradient:approx',cat('No gradient provided. Using an FD approximation instead (slow).\\n','It may be necessary to increase options.tolgradnorm.\\n','To disable this warning: warning(\'off\', \'manopt:getGradient:approx\')'))
        problem.approxgrad = copy(approxgradientFD(problem))
# conjugategradient/conjugategradient.m:151
    
    # Set local defaults here
    localdefaults.minstepsize = copy(1e-10)
# conjugategradient/conjugategradient.m:155
    localdefaults.maxiter = copy(1000)
# conjugategradient/conjugategradient.m:156
    localdefaults.tolgradnorm = copy(1e-06)
# conjugategradient/conjugategradient.m:157
    localdefaults.storedepth = copy(20)
# conjugategradient/conjugategradient.m:158
    # Changed by NB : H-S has the "auto restart" property.
# See Hager-Zhang 2005/2006 survey about CG methods.
# The auto restart comes from the 'max(0, ...)', not so much from the
# reason stated in Hager-Zhang I think. P-R also has auto restart.
    localdefaults.beta_type = copy('H-S')
# conjugategradient/conjugategradient.m:163
    localdefaults.orth_value = copy(Inf)
# conjugategradient/conjugategradient.m:164
    
    
    # Depending on whether the problem structure specifies a hint for
# line-search algorithms, choose a default line-search that works on
# its own (typical) or that uses the hint.
    if logical_not(canGetLinesearch(problem)):
        localdefaults.linesearch = copy(linesearch_adaptive)
# conjugategradient/conjugategradient.m:171
    else:
        localdefaults.linesearch = copy(linesearch_hint)
# conjugategradient/conjugategradient.m:173
    
    # Merge global and local defaults, then merge w/ user options, if any.
    localdefaults=mergeOptions(getGlobalDefaults(),localdefaults)
# conjugategradient/conjugategradient.m:177
    if logical_not(exist('options','var')) or isempty(options):
        options=struct()
# conjugategradient/conjugategradient.m:179
    
    options=mergeOptions(localdefaults,options)
# conjugategradient/conjugategradient.m:181
    # For convenience
    inner=problem.M.inner
# conjugategradient/conjugategradient.m:184
    lincomb=problem.M.lincomb
# conjugategradient/conjugategradient.m:185
    timetic=tic()
# conjugategradient/conjugategradient.m:187
    # If no initial point x is given by the user, generate one at random.
    if logical_not(exist('x','var')) or isempty(x):
        x=problem.M.rand()
# conjugategradient/conjugategradient.m:191
    
    # Create a store database and generate a key for the current x
    storedb=StoreDB(options.storedepth)
# conjugategradient/conjugategradient.m:195
    key=storedb.getNewKey()
# conjugategradient/conjugategradient.m:196
    # Compute cost-related quantities for x
    cost,grad=getCostGrad(problem,x,storedb,key,nargout=2)
# conjugategradient/conjugategradient.m:199
    gradnorm=problem.M.norm(x,grad)
# conjugategradient/conjugategradient.m:200
    Pgrad=getPrecon(problem,x,grad,storedb,key)
# conjugategradient/conjugategradient.m:201
    gradPgrad=inner[x,grad,Pgrad]
# conjugategradient/conjugategradient.m:202
    # Iteration counter (at any point, iter is the number of fully executed
# iterations so far)
    iter=0
# conjugategradient/conjugategradient.m:206
    # Save stats in a struct array info and preallocate.
    stats=savestats()
# conjugategradient/conjugategradient.m:209
    info[1]=stats
# conjugategradient/conjugategradient.m:210
    info[min(10000,options.maxiter + 1)].iter = copy([])
# conjugategradient/conjugategradient.m:211
    if options.verbosity >= 2:
        fprintf(' iter\\t               cost val\\t    grad. norm\\n')
    
    # Compute a first descent direction (not normalized)
    desc_dir=lincomb[x,- 1,Pgrad]
# conjugategradient/conjugategradient.m:219
    # Start iterating until stopping criterion triggers
    while true:

        # Display iteration information
        if options.verbosity >= 2:
            fprintf('%5d\\t%+.16e\\t%.8e\\n',iter,cost,gradnorm)
        # Start timing this iteration
        timetic=tic()
# conjugategradient/conjugategradient.m:231
        stop,reason=stoppingcriterion(problem,x,options,info,iter + 1,nargout=2)
# conjugategradient/conjugategradient.m:234
        if logical_not(stop) and abs(stats.stepsize) < options.minstepsize:
            stop=copy(true)
# conjugategradient/conjugategradient.m:238
            reason=sprintf(cat('Last stepsize smaller than minimum ','allowed; options.minstepsize = %g.'),options.minstepsize)
# conjugategradient/conjugategradient.m:239
        if stop:
            if options.verbosity >= 1:
                fprintf(cat(reason,'\\n'))
            break
        # The line search algorithms require the directional derivative of the
    # cost at the current point x along the search direction.
        df0=inner[x,grad,desc_dir]
# conjugategradient/conjugategradient.m:254
        # negative gradient. Equivalent to resetting the CG direction to a
    # steepest descent step, which discards the past information.
        if df0 >= 0:
            # Or we switch to the negative gradient direction.
            if options.verbosity >= 3:
                fprintf(cat('Conjugate gradient info: got an ascent direction ','(df0 = %2e), reset to the (preconditioned) ','steepest descent direction.\\n'),df0)
            # Reset to negative gradient: this discards the CG memory.
            desc_dir=lincomb[x,- 1,Pgrad]
# conjugategradient/conjugategradient.m:268
            df0=- gradPgrad
# conjugategradient/conjugategradient.m:269
        # Execute line search
        stepsize,newx,newkey,lsstats=options.linesearch(problem,x,desc_dir,cost,df0,options,storedb,key,nargout=4)
# conjugategradient/conjugategradient.m:275
        newcost,newgrad=getCostGrad(problem,newx,storedb,newkey,nargout=2)
# conjugategradient/conjugategradient.m:280
        newgradnorm=problem.M.norm(newx,newgrad)
# conjugategradient/conjugategradient.m:281
        Pnewgrad=getPrecon(problem,newx,newgrad,storedb,newkey)
# conjugategradient/conjugategradient.m:282
        newgradPnewgrad=inner[newx,newgrad,Pnewgrad]
# conjugategradient/conjugategradient.m:283
        # This paper https://www.math.lsu.edu/~hozhang/papers/cgsurvey.pdf
	# by Hager and Zhang lists many known beta rules. The rules defined
    # here can be found in that paper (or are provided with additional
    # references), adapted to the Riemannian setting.
	#
        if strcmpi(options.beta_type,'steep') or strcmpi(options.beta_type,'S-D'):
            beta=0
# conjugategradient/conjugategradient.m:296
            desc_dir=lincomb[x,- 1,Pnewgrad]
# conjugategradient/conjugategradient.m:297
        else:
            oldgrad=problem.M.transp(x,newx,grad)
# conjugategradient/conjugategradient.m:301
            orth_grads=inner[newx,oldgrad,Pnewgrad] / newgradPnewgrad
# conjugategradient/conjugategradient.m:302
            # survey on conjugate gradient methods, for example)
            if abs(orth_grads) >= options.orth_value:
                beta=0
# conjugategradient/conjugategradient.m:307
                desc_dir=lincomb[x,- 1,Pnewgrad]
# conjugategradient/conjugategradient.m:308
            else:
                desc_dir=problem.M.transp(x,newx,desc_dir)
# conjugategradient/conjugategradient.m:312
                if 'F-R' == upper(options.beta_type):
                    beta=newgradPnewgrad / gradPgrad
# conjugategradient/conjugategradient.m:317
                else:
                    if 'P-R' == upper(options.beta_type):
                        # vector grad(new) - transported grad(current)
                        diff=lincomb[newx,1,newgrad,- 1,oldgrad]
# conjugategradient/conjugategradient.m:321
                        ip_diff=inner[newx,Pnewgrad,diff]
# conjugategradient/conjugategradient.m:322
                        beta=ip_diff / gradPgrad
# conjugategradient/conjugategradient.m:323
                        beta=max(0,beta)
# conjugategradient/conjugategradient.m:324
                    else:
                        if 'H-S' == upper(options.beta_type):
                            diff=lincomb[newx,1,newgrad,- 1,oldgrad]
# conjugategradient/conjugategradient.m:327
                            ip_diff=inner[newx,Pnewgrad,diff]
# conjugategradient/conjugategradient.m:328
                            beta=ip_diff / inner[newx,diff,desc_dir]
# conjugategradient/conjugategradient.m:329
                            beta=max(0,beta)
# conjugategradient/conjugategradient.m:330
                        else:
                            if 'H-Z' == upper(options.beta_type):
                                diff=lincomb[newx,1,newgrad,- 1,oldgrad]
# conjugategradient/conjugategradient.m:333
                                Poldgrad=problem.M.transp(x,newx,Pgrad)
# conjugategradient/conjugategradient.m:334
                                Pdiff=lincomb[newx,1,Pnewgrad,- 1,Poldgrad]
# conjugategradient/conjugategradient.m:335
                                deno=inner[newx,diff,desc_dir]
# conjugategradient/conjugategradient.m:336
                                numo=inner[newx,diff,Pnewgrad]
# conjugategradient/conjugategradient.m:337
                                numo=numo - dot(dot(2,inner[newx,diff,Pdiff]),inner[newx,desc_dir,newgrad]) / deno
# conjugategradient/conjugategradient.m:338
                                beta=numo / deno
# conjugategradient/conjugategradient.m:340
                                desc_dir_norm=problem.M.norm(newx,desc_dir)
# conjugategradient/conjugategradient.m:343
                                eta_HZ=- 1 / (dot(desc_dir_norm,min(0.01,gradnorm)))
# conjugategradient/conjugategradient.m:344
                                beta=max(beta,eta_HZ)
# conjugategradient/conjugategradient.m:345
                            else:
                                error(cat('Unknown options.beta_type. ','Should be steep, S-D, F-R, P-R, H-S or H-Z.'))
                desc_dir=lincomb[newx,- 1,Pnewgrad,beta,desc_dir]
# conjugategradient/conjugategradient.m:352
        # Make sure we don't use too much memory for the store database
        storedb.purge()
        x=copy(newx)
# conjugategradient/conjugategradient.m:362
        key=copy(newkey)
# conjugategradient/conjugategradient.m:363
        cost=copy(newcost)
# conjugategradient/conjugategradient.m:364
        grad=copy(newgrad)
# conjugategradient/conjugategradient.m:365
        Pgrad=copy(Pnewgrad)
# conjugategradient/conjugategradient.m:366
        gradnorm=copy(newgradnorm)
# conjugategradient/conjugategradient.m:367
        gradPgrad=copy(newgradPnewgrad)
# conjugategradient/conjugategradient.m:368
        iter=iter + 1
# conjugategradient/conjugategradient.m:371
        stats=savestats()
# conjugategradient/conjugategradient.m:374
        info[iter + 1]=stats
# conjugategradient/conjugategradient.m:375

    
    info=info[1:iter + 1]
# conjugategradient/conjugategradient.m:380
    if options.verbosity >= 1:
        fprintf('Total time is %f [s] (excludes statsfun)\\n',info[end()].time)
    
    # Routine in charge of collecting the current iteration stats
    
@function
def savestats(*args,**kwargs):
    varargin = savestats.varargin
    nargin = savestats.nargin

    stats.iter = copy(iter)
# conjugategradient/conjugategradient.m:389
    stats.cost = copy(cost)
# conjugategradient/conjugategradient.m:390
    stats.gradnorm = copy(gradnorm)
# conjugategradient/conjugategradient.m:391
    if iter == 0:
        stats.stepsize = copy(nan)
# conjugategradient/conjugategradient.m:393
        stats.time = copy(toc(timetic))
# conjugategradient/conjugategradient.m:394
        stats.linesearch = copy([])
# conjugategradient/conjugategradient.m:395
        stats.beta = copy(0)
# conjugategradient/conjugategradient.m:396
    else:
        stats.stepsize = copy(stepsize)
# conjugategradient/conjugategradient.m:398
        stats.time = copy(info[iter].time + toc(timetic))
# conjugategradient/conjugategradient.m:399
        stats.linesearch = copy(lsstats)
# conjugategradient/conjugategradient.m:400
        stats.beta = copy(beta)
# conjugategradient/conjugategradient.m:401
    
    stats=applyStatsfun(problem,x,storedb,key,options,stats)
# conjugategradient/conjugategradient.m:403
    return stats
    
if __name__ == '__main__':
    pass
    
    return stats
    
if __name__ == '__main__':
    pass
    