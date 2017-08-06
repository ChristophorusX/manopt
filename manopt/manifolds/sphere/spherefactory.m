function M = spherefactory(n, m)
% Returns a manifold struct to optimize over unit-norm vectors or matrices.
%
% function M = spherefactory(n)
% function M = spherefactory(n, m)
%
% Manifold of n-by-m real matrices of unit Frobenius norm.
% By default, m = 1, which corresponds to the unit sphere in R^n. The
% metric is such that the sphere is a Riemannian submanifold of the space
% of nxm matrices with the usual trace inner product, i.e., the usual
% metric.
% 
% See also: obliquefactory spherecomplexfactory

% This file is part of Manopt: www.manopt.org.
% Original author: Nicolas Boumal, Dec. 30, 2012.
% Contributors: 
% Change log: 
%
%   Oct. 8, 2016 (NB)
%       Code for exponential was simplified to only treat the zero vector
%       as a particular case.
%
%   Oct. 22, 2016 (NB)
%       Distance function dist now significantly more accurate for points
%       within 1e-7 and less from each other.

    
    if ~exist('m', 'var')
        m = 1;
    end

    if m == 1
        M.name = @() sprintf('Sphere S^%d', n-1);
    else
        M.name = @() sprintf('Unit F-norm %dx%d matrices', n, m);
    end
    
    M.dim = @() n*m-1;
    
    M.inner = @(x, d1, d2) d1(:).'*d2(:);
    
    M.norm = @(x, d) norm(d, 'fro');
    
    M.dist = @dist;
    function d = dist(x, y)
        
        % This computation cannot be accurate below an output of 2e-8.
        % The reason is: if two unit-norm vectors x and y are very close to
        % one another, their inner product is about 1. The machine
        % precision at 1 is eps(1) = 2e-16. The correct value for
        % acos(1-eps(1)) is about 2e-8. This can be checked with the
        % syms toolbox: syms x; f = acos(1-x); vpa(subs(f, x, eps(1)), 32)
        % Thus, if x and y are actually closer to each other than 2e-8,
        % their inner product will be even closer to 1, but that cannot be
        % represented in IEEE arithmetic. Thus, their inner product will be
        % rounded to either 1 (giving 0 distance) or to 1-eps(1), which
        % gives a distance of 2e-8, or to something even further from 1. No
        % distance between 0 and 2e-8 can thus be computed this way.
        d = real(acos(x(:).'*y(:)));
        
        % Hence, if the distance proves dangerously small so that it is
        % possible that we suffered from round-off, we compute the distance
        % in the embedding space instead. At this scale, this is quite
        % accurate.
        if d < 1e-7
            d = norm(x-y, 'fro');
        end
    end
    
    M.typicaldist = @() pi;
    
    M.proj = @(x, d) d - x*(x(:).'*d(:));
    
    M.tangent = M.proj;
	
    % For Riemannian submanifolds, converting a Euclidean gradient into a
    % Riemannian gradient amounts to an orthogonal projection.
	M.egrad2rgrad = M.proj;
	
	M.ehess2rhess = @ehess2rhess;
	function rhess = ehess2rhess(x, egrad, ehess, u)
        rhess = M.proj(x, ehess) - (x(:)'*egrad(:))*u;
	end
    
    M.exp = @exponential;
    
    M.retr = @retraction;

    M.log = @logarithm;
    function v = logarithm(x1, x2)
        v = M.proj(x1, x2 - x1);
        di = M.dist(x1, x2);
        % If the two points are "far apart", correct the norm.
        if di > 1e-6
            nv = norm(v, 'fro');
            v = v * (di / nv);
        end
    end
    
    M.hash = @(x) ['z' hashmd5(x(:))];
    
    M.rand = @() random(n, m);
    
    M.randvec = @(x) randomvec(n, m, x);
    
    M.lincomb = @matrixlincomb;
    
    M.zerovec = @(x) zeros(n, m);
    
    M.transp = @(x1, x2, d) M.proj(x2, d);
    
    M.pairmean = @pairmean;
    function y = pairmean(x1, x2)
        y = x1+x2;
        y = y / norm(y, 'fro');
    end

    M.vec = @(x, u_mat) u_mat(:);
    M.mat = @(x, u_vec) reshape(u_vec, [n, m]);
    M.vecmatareisometries = @() true;

end

% Exponential on the sphere
function y = exponential(x, d, t)

    if nargin == 2
        % t = 1
        td = d;
    else
        td = t*d;
    end
    
    nrm_td = norm(td, 'fro');
    
    % Former versions of Manopt avoided the computation of sin(a)/a for
    % small a, but further investigations suggest this computation is
    % well-behaved numerically.
    if nrm_td > 0
        y = x*cos(nrm_td) + td*(sin(nrm_td)/nrm_td);
    else
        y = x;
    end

end

% Retraction on the sphere
function y = retraction(x, d, t)

    if nargin == 2
        % t = 1;
        td = d;
    else
        td = t*d;
    end
    
    y = x + td;
    y = y / norm(y, 'fro');

end

% Uniform random sampling on the sphere.
function x = random(n, m)

    x = randn(n, m);
    x = x / norm(x, 'fro');

end

% Random normalized tangent vector at x.
function d = randomvec(n, m, x)

    d = randn(n, m);
    d = d - x*(x(:).'*d(:));
    d = d / norm(d, 'fro');

end
