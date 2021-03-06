function [m, nextr] = adaptRobust(prevr, Uconst, Lconst, keps, A, B)
% [m, nextr] = adaptRobust(prevr, Uconst, Lconst, keps, A, B)
%   Returns next mesh point location and its reach to the right
%   Inputs are previous r location, const for Upper bound (should be -theta + V - log L), const for Lower
%   bound (should be -theta -W + log U); we return s.t. within keps=k.epsilon
%   Assumes we're working from left (low) to right (high)
%   With left bound A, right boung 1-B
%
%   Use a proper nonlinear root-finder.

% For now use loose upper bound, i.e. assume integral is max val * interval
% No seems not much harder to try integral so we'll do that
% ^^ Great idea, Adrian!

global MESH_TOL

fudge = 1e-6;

insurance = 1e-4;

if (1 - B - A) < MESH_TOL
    %warning('Interval size 1 - B - A is already smaller than fzero tolerance. Setting m = nextr = 1 - B.');
    m = 1 - B;
    nextr = m;
    return;
end

C=Uconst-Lconst;
if C<0
    error('Error: Upper bound constant < Lower bound constant\n');    
end
if A<0 || (1-B)>1 || A>(1-B)
    error('Error: invalid A, 1-B bounds\n');    
end
if prevr<A || prevr>1-B
    error('Error: Previous r not in [A,1-B]\n');    
end
if Uconst+log(prevr/(1-prevr)) < 0 - fudge
    error('Error: Upper bound<0 at prevr\n');    
end
if Lconst+log((1-B)/B) > 0 + fudge
    error('Error: Lower bound>0 at right edge 1-B\n');    
end

p = prevr;

% opts = struct('TolX', MESH_TOL, 'TolFun', MESH_TOL);


%der   = @(t) log(t) - log1p(-t);

%ubInt = @(t) intp2m(Uconst, t, p) - keps;
if intNew(1 - B, Uconst, -keps, p) < 0
%if ubInt(p) < 0 && ubInt(1 - B) < 0
    % Both endpoints of the feasible interval are negative, so fzero won't
    % find a root. But this means the next mesh point to be added goes over
    % the end (1 - B). So just add that endpoint and be done.
    m     = 1 - B;
    nextr = 1 - B;
    return;
else
    % NOTE: Matlab's fzero dominates 'bisection', even though it does other
    % junk.
    %
    % Evaluating the ubInt function is actually a bit expensive. The
    % adaptive and fdm could be moved to a mexFile... ick, and the mesh
    % construction too.
    
%      m     = fzero(ubInt, [p 1-B], opts);
%     m = fzero_fast([p 1-B], MESH_TOL, Uconst, -keps, p);
    m = fzero_fast(@intNew, [p 1-B], MESH_TOL, Uconst, -keps, p);    
%     mSlow = fzero(@intNew, [p 1-B], opts, Uconst, -keps, p);
%     assert(m == mSlow);

    %m = steffensen(ubInt, p, MESH_TOL);
    %m = newton(ubInt, der, p, MESH_TOL);
    %m = bisection(ubInt, p, 1 - B, MESH_TOL);    
end

%lbInt = @(t) intp2m(Lconst, t, m) + keps;
%if lbInt(m) > 0 && lbInt(1 - B) > 0
if intNew(1 - B, Lconst, keps, m) > 0
    % Both endpoints of the feasible interval are negative, so fzero won't
    % find a root. But this means the mesh point m reaches over the end, so
    % we'll just set the reach to be exactly the end and be done.
    nextr = 1 - B;
else
%     nextr = fzero(lbInt, [m 1 - B], opts);
%     nextr = fzero_fast([m 1-B], MESH_TOL, Lconst, keps, m);
    nextr = fzero_fast(@intNew, [m 1-B], MESH_TOL, Lconst, keps, m);
    
%     nextrSlow = fzero(@intNew, [m 1-B], opts, Lconst, keps, m);
%     assert(nextr == nextrSlow);
    
    %nextr = steffensen(lbInt, m, MESH_TOL);
    %nextr = newton(lbInt, der, m, MESH_TOL);
    %nextr = bisection(lbInt, m, 1 - B, MESH_TOL);
end

%fprintf('adaptRobust prevr = %g, m = %g, nextr = %g\n', prevr, m, nextr);

assert(m >= prevr && nextr >= m, 'integrals were not right');
%assert(nextr > prevr, 'No progress made!');

end

