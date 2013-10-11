function [A,B] = bpbound(theta, W, thresh)
% function [A,B] = bpbound(N,theta,W,maxiter)
% Improved bound calculation algorithm
% section IIIE of http://cs.ru.nl/~jorism/articles/04385778.pdf

% Initialize bounds on cavity fields

if nargin == 2
    thresh = 1e-6;
end

N = length(theta);

eta_lo = zeros(N,N);
eta_hi = zeros(N,N);

% Convert from {0, 1} to {-1, 1}
% A bit more cumbersome than just inverting Mooij's equations; just write
% out the grid notation
J = 0.25 * W;
eta = 0.5*theta + sum(J, 2);

for i=1:N
  for j=1:N
    if J(i,j)
      eta_lo(i,j) = -Inf;
      eta_hi(i,j) = Inf;
    end
  end
end

% Propagate bounds on cavity fields

while true
  for i=1:N
    for j=1:N
      if J(i,j)
        eta_lo_new(i,j) = eta(i);
        eta_hi_new(i,j) = eta(i);
        for k=1:N
          if J(i,k) && k ~= j
            a = atanh(tanh(J(k,i)) * tanh(eta_lo(k,i)));
            b = atanh(tanh(J(k,i)) * tanh(eta_hi(k,i)));
            eta_lo_new(i,j) = eta_lo_new(i,j) + min([a b]);
            eta_hi_new(i,j) = eta_hi_new(i,j) + max([a b]);                                   
          end
        end
      end
    end
  end
  
  if sum(abs(eta_lo - eta_lo_new)) + sum(abs(eta_hi - eta_hi_new)) < thresh
      break;
  else    
      eta_lo = eta_lo_new;
      eta_hi = eta_hi_new;
  end
end

% Calculate beliefs in tanh parameterization
beta_lo = eta;
beta_hi = eta;
for i=1:N
  for k=1:N
    if J(i,k)
      a = atanh(tanh(J(k,i)) * tanh(eta_lo(k,i)));
      b = atanh(tanh(J(k,i)) * tanh(eta_hi(k,i)));
      beta_lo(i) = beta_lo(i) + min([a b]);
      beta_hi(i) = beta_hi(i) + max([a b]);
    end
  end
end

% Calculate A,B bounds on pseudomarginal
A = (tanh(beta_lo) + 1) / 2;
B = 1 - ((tanh(beta_hi) + 1) / 2);

return
