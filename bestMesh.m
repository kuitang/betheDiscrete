function [gam, complexity] = bestMesh(theta, W, epsilon, verbose)
% [gams, complexity] = bestMesh(theta, W, epsilon)
%
%   NOTE: JUST USE ADAPTIVEMINSUM. For our regime (not "too" many nodes,
%   this is the best anyways).
%
%   Try all our mesh methods (1st deriv, second deriv, Lagrangian, simple)
%   and return the mesh grid for one with smallest total points.
%
%   theta is an N-vector, W is a weight matrix (symmetric), epsilon is the
%   accuracy threshold.
%
%   gams is an N-cell array, each element is the ordered vector of mesh
%   grids in that variable; complexity is a structue with fields sumN,
%   thisN, prodN.

    if nargin == 3
        verbose = false;
    end
        
    % Get the gams from somewhere...
    %[A, B] = bpbound(theta, W);  
    
    %[A, B] = BBPNew(theta, W);        
    [A, B, alpha, L, U] = BBPNew(theta, W);
    
    %fmMethods = { 'simple', 'minsum', 'adaptivesimple', 'adaptiveminsum' };
    
    %fmMethods = { 'adaptiveminsum' };
    fmMethods = { 'minsum' };
    

    
    [gams, sumN, prodN, thisN] = cellfun(@(m) fdm(theta, W, A, B, epsilon, m, L, U), fmMethods, 'UniformOutput', false);
    [~, i] = min(cell2mat(sumN));
    %i = 4;
    gam = gams{i};
    complexity = struct('sumN', sumN{i}, 'thisN', thisN{i}, 'prodN', prodN{i});
        
    if thisN{i} >= 1000
        warning('Problem may have high complexity:');
        complexity
    end
    
    if verbose
        disp(['Methods: ' fmMethods]);
        disp(['sumN: ' sumN]);
        disp(['prodN: ' prodN]);
        fprintf(1, 'Best method: %s, sumN = %d\n', fmMethods{i}, sumN{i});
    end        

end

