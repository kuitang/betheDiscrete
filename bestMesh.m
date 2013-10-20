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

    fudge = 10 * eps;

    if nargin == 3
        verbose = false;
    end
        
    % Get the gams from somewhere...
%     [Amk, Bmk] = MKNew(theta, W);      

        
    [Amk, Bmk] = MKNew_mex(theta, sparse(W));
%     fprintf('A 1-norm deviation from MK and MK MEX: %g\n', norm(Amk - Amex, 1));
%     fprintf('B 1-norm deviation from MK and MK MEX: %g\n', norm(Bmk - Bmex, 1));   
    Smk = 1 - Bmk - Amk;
    
%     [Abp, Bbp] = BBPNew(theta, W);
%     Sbp = 1 - Bbp - Abp;
    
    assert(all(Smk > -fudge), 'Smk was less than 0 (over fudge)!');
%     assert(all(Sbp > -fudge), 'Sbp was less than 0 (over fudge)!');
    
    % Set the flipped quantities to the upper bound.
    mkFlip = Smk < 0;
    Amk(mkFlip) = 1 - Bmk(mkFlip);
    
%     bpFlip = Sbp < 0;
%     Abp(bpFlip) = 1 - Bbp(bpFlip);
    
%     if sum(Sbp) < sum(Smk)
%         warning('Strangely, Sbp = %g < Smk = %g');
%         A = Abp; B = Bbp;        
%     else
%         A = Amk; B = Bmk;        
%     end
    A = Amk;
    B = Bmk;
    N = length(A);
        
    [gams, sumN, prodN, thisN] = fdm(theta, W, A, B, epsilon, 'minsum');        
    
%    gam = cell(N, 1);
%    for n = 1:N
%        lb = A(n);
%        ub = 1 - B(n);
%        gam{n} = [lb:gams(n):ub ub];
%        assert(all(diff(gam{n}) <= gams(n) + 10*eps), 'Did not cover!');
%    end
%
    [gams, sumN, prodN, thisN] = fdm(theta, W, A, B, epsilon, 'adaptiveminsum');  

    gam = gams;
    
    complexity = struct('sumN', sumN, 'thisN', thisN, 'prodN', prodN);        
    if sumN >= 1000
        warning('Problem may have high complexity: sumN = %d', sumN);
    end        

end

