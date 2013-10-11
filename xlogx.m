function [ f ] = xlogx( x )
%[ f ] = xlogx( x )
%   Returns x * log x, using natural log and returns 0 if x=0

% log(1) = 0 still

f = x .* log(x);
f(x == 0) = 0;

end

