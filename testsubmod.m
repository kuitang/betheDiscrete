dbstop if error;

N = 12;
Ntrials = 100;

for t = 1:Ntrials
    %% One trial
    W     = rand(N, N);
    W(1:4:end) = 0;
    W = triu(W, 1) + triu(W, 1)';

    theta = randn(N, 1) - sum(W, 2) / 2;
    
    %% Sub
    theta
    W

    [energy, xmap, misc] = submodularMAP_mex(theta, sparse(W))

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
