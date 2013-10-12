function [ ] = mySubmod(N)

    % Standard constructions
    W = rand(N, N)
    W(1:(N+1):end) = 0;
    theta = randn(N, 1) - sum(W, 2) / 2
    
    % Upper bound
    maxB = -min(theta) + 0.1;
    % Lower bound -- would just set all zeros to saturate. Though this
    % isn't exactly the right answer...
    minB = -max(theta) - 1;
    
    Nbs = 1000;
    bs = linspace(minB, maxB, Nbs);
    minCuts = zeros(Nbs, 1);
    nCount  = zeros(Nbs, 1);
    
    xSeq{1} = zeros(N, 1);
    
    for n = 1:Nbs        
        [e, x, m] = submodularMAP_mex(theta + bs(n), sparse(W));
        x = double(x);
        minCuts(n) = m.maxFlow;
        assert(all(x - xSeq{end} >= 0), 'not monotonic');        
            
        xSeq{end+1} = x;                    
        nCounts(n) = sum(x);
    end
    
    %nCounts  ./ max(minCuts); % rescale for display
    plot(bs, minCuts, bs, nCounts);
    
end
