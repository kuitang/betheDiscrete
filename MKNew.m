function [A,B] = MKNew(theta, W, thresh)
% function [A,B] = MKNew(N,theta,W,maxiter)
% Improved bound calculation algorithm
% section IIIE of http://cs.ru.nl/~jorism/articles/04385778.pdf

% Initialize bounds on cavity fields

if nargin == 2
    thresh = 1e6;
end

N = length(theta);

% Convert from {0, 1} to {-1, 1}
% A bit more cumbersome than just inverting Mooij's equations; just write
% out the grid notation
J = 0.25 * W;
eta = 0.5*theta + sum(J, 2);

etaLo = -Inf * spones(W);
etaHi = Inf  * spones(W);

[iVec, jVec, wVec] = find(W);
nI = length(iVec);

% Propagate bounds on cavity fields

nIter = 0;

JT = transpose(J);

etaLoNew = spones(W);
etaHiNew = spones(W);

while true        
    for ne = 1:nI
        i = iVec(ne);
        j = jVec(ne);
        
        etaLoNew(i,j) = eta(i);
        etaHiNew(i,j) = eta(i);
        
        kidxs = full(JT(:,i) ~= 0);
        kkidxs = full(J(:,i) ~= 0);
        assert(all(kidxs == kkidxs));        
        kidxs(j) = 0;
        
        a = atanh(tanh(J(kidxs,i)) .* tanh(etaLo(kidxs,i)));
        b = atanh(tanh(J(kidxs,i)) .* tanh(etaHi(kidxs,i)));        

%         fprintf('DEBUG: (i,j) = (%d,%d); kidxs = %s\n', i, j, num2str(find(kidxs)));
        
        dLo = sum(min(a, b));
        dHi = sum(max(a, b));
        
        etaLoNew(i,j) = etaLoNew(i,j) + dLo;
        etaHiNew(i,j) = etaHiNew(i,j) + dHi;
                        
%         for k = find(J(i,:))
%             if k ~= j
%                 a = atanh(tanh(J(k,i)) * tanh(etaLo(k,i)));
%                 b = atanh(tanh(J(k,i)) * tanh(etaHi(k,i)));
%                 etaLoNew(i,j) = etaLoNew(i,j) + min([a b]);
%                 etaHiNew(i,j) = etaHiNew(i,j) + max([a b]);                
%             end
%         end        
    end
    
%     fprintf('etaLo, etaLoNew\n');
%     [full(etaLo(:)) full(etaLoNew(:))]
%     
%     fprintf('etaHi, etaHiNew\n');
%     [full(etaHi(:)) full(etaHiNew(:))]
    
    % Yeah, this is slow, but it's the correct convergence criterion.
    conv = sum(abs(etaLo - etaLoNew)) + sum(abs(etaHi - etaHiNew))
    if conv < thresh
        break;
    end
    
    etaLo = etaLoNew;
    etaHi = etaHiNew;
    
    %fprintf(1, 'bpbound: iter %d, max(conv) = %g; thresh = %g\n', nIter, max(conv), thresh);
    nIter = nIter + 1;
end

% Intermedite check
[nonzeros(etaLo) nonzeros(etaHi)]

% Calculate beliefs in tanh parameterization
betaLo = eta;
betaHi = eta;

for ne = 1:nI
    i = iVec(ne);
    k = jVec(ne);
    
    a = atanh(tanh(J(k,i)) * tanh(etaLo(k,i)));
    b = atanh(tanh(J(k,i)) * tanh(etaHi(k,i)));
    
    betaLo(i) = betaLo(i) + min([a b]);
    betaHi(i) = betaHi(i) + max([a b]);      
end

[betaLo betaHi]

% Calculate A,B bounds on pseudomarginal
A = (tanh(betaLo) + 1) / 2;
B = 1 - ((tanh(betaHi) + 1) / 2);

if any(A == 0) || any(B == 0)
    warning('bpbound:zero', 'Got zero bounds: %d in A and %d in B', sum(A == 0), sum(B == 0));
end

return
