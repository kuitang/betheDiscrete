function [ xs, fill_n ] = enumerate( limits, level, prefix, fill_n, xs)
% xs = enumerate(limits) generates all x values between 1 and limits.
%
% limits = vector of upper limits, one per dimension.
%
% The rest are recursion arguments; ignore them.

    D = length(limits);

    if nargin == 1
        N = prod(limits);
        xs = zeros(N, D);       
        level = 1;
        fill_n = 0;
    end
    
    % Leaf
    if level == D
        nnew = limits(end);
        next_fill_n = fill_n + nnew;
        idxs = (fill_n+1):next_fill_n;
        xs(idxs, 1:(D-1)) = repmat(prefix, nnew, 1);
        xs(idxs, D) = (1:nnew)';
        
        fill_n = next_fill_n;            
    else    
        for i = 1:limits(level)
            prefix(level) = i;
            [xs, fill_n] = enumerate(limits, level + 1, prefix, fill_n, xs);
        end
    end

end

