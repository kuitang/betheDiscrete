dbstop if error;

N = 12;
Ntrials = 100;

%% barrier
for t = 1:Ntrials
    %% One trial
    W     = rand(N, N);
    W(1:4:end) = 0;
    W = triu(W, 1) + triu(W, 1)';

    theta = randn(N, 1) - sum(W, 2) / 2;
    
    nEdges = N * (N - 1) / 2;
    pots = cell(nEdges, 1);
    pW = zeros(N, N);    
    idx = 1;
    for j = 1:N
        for i = 1:(j-1)
            pW(i,j) = idx;
            pots{idx} = [0 0 ; 0 -W(i,j)];
            idx = idx + 1;
        end
    end
    
    thetaFull = [zeros(N,1) -theta];
    
    [energy, xmap, misc] = submodularMAP_mex(theta, sparse(W))    
    
    thetaFull    
    pots{:}
    [energyFull, xmapFull, miscFull] = submodularMAPFull_mex(thetaFull, sparse(pW), pots);
    thetaFull    
    pots{:}
    
    assertElementsAlmostEqual(energy, energyFull);
    

    % Brute force check
    xs = enumerate(2 * ones(N, 1)) - 1;
    nxs = size(xs, 1);

    e(nxs) = 0;
    for n = 1:nxs
        e(n) = mrfEnergy(theta, W, xs(n,:)');
    end

    [eMin, iMin] = min(e);
    
    assertElementsAlmostEqual(eMin, energy);
end
