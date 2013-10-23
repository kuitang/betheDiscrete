function [gam, cpx, epsilon] = bestMesh(theta, W, minPoints, maxPoints)
% [gams, complexity, epsilon] = bestMesh(theta, W, minPoints, maxPoints)
%
%   Return the best mesh and corresponding epsilon (doubling and halving)
%   to attain a mesh with size in [minPoints, maxPoints].
%
%   This just uses adaptiveminsum.

    global MESH_TOL MESH_MAX_POINTS MIN_EPS
    MESH_TOL = 1e6;
    MIN_EPS  = 1e-4;    
    
    fudge = 1e-10;
    
    assert(minPoints < maxPoints, 'minPoints is not less than maxPoints');
            
    [A, B] = MKNew_mex(theta, sparse(W));
    % HACK: Remove zeros and 1s. Not really justified, but what the hell,
    % we decided the model anyways.
    A(A == 0) = eps;
    A(A == 1) = 1 - eps;
    B(B == 0) = eps;
    B(B == 1 ) = eps;
    S = 1 - B - A;
    assert(all(S > -fudge), 'Smk was less than 0 (over fudge)!');
    N = length(A);
    
%     figure;
%     hist(S);
    
    if maxPoints < N
        warning('bestMesh:maxPoints', 'maxPoints < N, a contradiction. Setting maxPoints = N.');
        maxPoints = N;
    end
    MESH_MAX_POINTS = maxPoints;
        
    % Set the flipped quantities to the upper bound.
    mkFlip = S < 0;
    A(mkFlip) = 1 - B(mkFlip);

    
    % Bisectio with infeasible start.
    
    hiLoRatio = 100;
    loEps = 0.1;
    hiEps = hiLoRatio * loEps;
    
    [loGam, loCpx] = fdm(theta, W, A, B, loEps, 'adaptiveminsum');
    [hiGam, hiCpx] = fdm(theta, W, A, B, hiEps, 'adaptiveminsum');
            
    nBisects = 0;
    while true
        fprintf('nBisects = %d; loEps = %g; hiEps = %g\n', nBisects, loEps, hiEps);
        if loEps <= MIN_EPS        
            epsilon = loEps;
            gam = loGam;
            cpx = loCpx;
            return;        
        elseif loCpx.sumN < minPoints && hiCpx.sumN < minPoints % both too coarse; need lower epsilon.
            % Switch the low to the high.
            hiEps = loEps;
            hiGam = loGam;
            hiCpx = loCpx;
            
            loEps = loEps / hiLoRatio;
            [loGam, loCpx] = fdm(theta, W, A, B, loEps, 'adaptiveminsum');

        elseif isempty(loGam) && isempty(hiGam) % both too fine; need higher epsilon.
            % Switch the high to the low.
            loEps = hiEps;
            loGam = hiGam;
            loCpx = hiCpx;
            
            hiEps = hiLoRatio * hiEps;
            [hiGam, hiCpx] = fdm(theta, W, A, B, hiEps, 'adaptiveminsum'); 
        elseif hiCpx.sumN < minPoints && isempty(loGam) % hi is too coarse but lo is too fine.
            % Halve the upswitch.
            hiLoRatio = hiLoRatio / 2.0;
                        
            % Decrease hiEps from its current value, since we halved
            % hiLoRatio.            
            hiEps = hiLoRatio * loEps;
            assert(loEps < hiEps, 'You messed up hi/lo.');
            [hiGam, hiCpx] = fdm(theta, W, A, B, hiEps, 'adaptiveminsum');
        elseif loCpx.sumN >= minPoints && loCpx.sumN <= maxPoints
            % feasible!
            epsilon = loEps;
            gam = loGam;
            cpx = loCpx;
            return;
        elseif hiCpx.sumN >= minPoints && hiCpx.sumN <= maxPoints
            % feasible!
            epsilon = hiEps;
            gam = hiGam;
            cpx = hiCpx;
            return;
        else
            error('Unreachable!');
        end     
        
        nBisects = nBisects + 1;
    end
end

