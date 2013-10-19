function [ A, B, varargout ] = BBPNew(theta, W, varargin)
% BBP Bethe belief propagation, modified better suit first derivative bound
%
%   [A, B, /alpha/, /L/, /U/] = BBPNew(theta, W, 'thresh', 0.002, 'maxIter', 
%       Inf, 'A', Ainit, 'B', Binit) where /x/ denotes optional
%       output. A is the vector of lower bounds; B is the vector of
%       complementary upper bounds.
%
%       Note that you can use Ainit and Binit with the MK values to obtain
%       a (convergent) L, U.

    nNodes = size(W, 1);
    
    posW = zeros(nNodes, 1);
    negW = zeros(nNodes, 1);
    
    for j = 1:nNodes
        col = W(W(:,j) ~= 0,j);
        posW(j) = sum(col(col > 0));
        negW(j) = -sum(col(col < 0));
    end

    p = inputParser;    
    p.addRequired('theta');
    p.addRequired('W');
    p.addParamValue('thresh', 1e-6);
    p.addParamValue('maxIter', Inf);
    p.addParamValue('A', sigmoid(theta - negW));
    p.addParamValue('B', 1 - sigmoid(theta + posW));            
    p.parse(theta, W, varargin{:});
    
    o = p.Results;    
        
    alpha = exp(abs(W)) - 1;    
    
    converged = false;
            
    iter = 0;
    A = o.A;
    B = o.B;
    L = ones(nNodes, 1);
    U = ones(nNodes, 1);
    
    while ~converged && iter <= o.maxIter
        oldA = A;
        oldB = B;
        
        for j = 1:nNodes
            % We can't vectorize this because a new update (j) may access
            % our updated variable.
            L(j) = 1;
            U(j) = 1;
            
            for i = find(W(:,j))'
                a = alpha(i,j);
                if W(i,j) > 0
                    L(i) = L(i) * (1 + a*A(j) / (1 + a*(1 - B(i))*(1 - A(j))));
                    U(i) = U(i) * (1 + a*B(j) / (1 + a*(1 - A(i))*(1 - B(j))));
                else
                    L(i) = L(i) * (1 + a*B(j) / (1 + a*(1 - B(i))*(1 - B(j))));
                    U(i) = U(i) * (1 + a*A(j) / (1 + a*(1 - A(i))*(1 - A(j))));
                end
            end
                    
            A(j) = 1 / (1 + exp(-theta(j) + negW(j)) / L(j));
            B(j) = 1 / (1 + exp(theta(j)  + posW(j)) / U(j));            
        end
                
        % Lemma 9: At every iteration, each element of A, B
        % monotonically increase 
        %
        % NOTE: may fail due to fudge factor
%         assert(all(A >= oldA));
%         assert(all(B >= oldB));
        
        dA = abs(A - oldA);
        dB = abs(B - oldB);
        if all(dA < o.thresh) && all(dB < o.thresh)
            converged = true;
        elseif iter >= o.maxIter
            warning('BBPNew:converge', ...
                'Convergence threshold %d not reached after %d iterations.', ...
                o.thresh, o.maxIter);                
        end
        
        iter = iter + 1;
    end
    
    if nargout >= 3
        varargout{1} = alpha;
    end
    if nargout >= 4
        varargout{2} = L;
    end
    if nargout >= 5
        varargout{3} = U;
    end        
end

