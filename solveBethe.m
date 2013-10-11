function [logZ, oneMargs, twoMargs, meshSize, misc] = solveBethe(theta, W, varargin)
% [mu, xi, twoMargs, energy, complexity] = solveBetheExact(theta, W, epsilon)
%
%   Find (global?) minimum of the Bethe free energy.
%
%   Returns the singleton potentials in mu (mu(i) == P(x_i = 1)) and
%   pairwise (2x2) potentials in twoMargs, with xi(i,j) = P(x_i = x_j = 1)
%
%   Right now, parameterizes an ILP using our bounds. But this will seem
%   silly. We should use gradient information. What about alpha-BB?
%
%   complexity is a structure with fields solverTime, sumN, prodN, thisN
%   count the sum, product, and vector of mesh points.
%
%   NOTE: Energy is -log Z. But the outer (dual decomposition) loop
%   introduces some constant shift, which it takes care of.
%

    % FILL IN THE DETAILS;    
    
    p = inputParser;
    p.addRequired('theta', @isnumeric);
    p.addRequired('W', @isnumeric);    
    p.addParamValue('epsilon', 0.1); % It's gonna be slow...
        
    p.parse(theta, W, varargin{:});    
    o = p.Results;
    
    W = sparse(triu(W, 1) + triu(W, 1)');
        
    [gams, meshSize] = bestMesh(theta, W, o.epsilon);  
    [logZ, oneMargs, twoMargs, misc] = BetheGams_mex(theta, W, gams);

end

