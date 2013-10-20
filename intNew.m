function [ integral ] = intNew(m, Uconst, keps, p)
%[ integral ] = intp2m( Uconst,m,p )
%   Returns integral of Uconst + log[q/(1-q)] from m to p
% integralSlow = Uconst*(m-p) + xlogx(m) + xlogx(1-m) - xlogx(p) - xlogx(1-p) + keps;

% xlogx is a bottleneck; optimize it here.
% x log x + x log(1 - x) is 0 iff x = 0 or x = 1.

if m == 0 || m == 1
    mlog = 0;
else
    mlog = m * log(m) + (1 - m) * log1p(-m);
end

if p == 0 || p == 1
    plog = 0;
else
    plog = p * log(p) + (1 - p) * log1p(-p);
end

integral = Uconst * (m - p) + mlog - plog + keps;

% assertElementsAlmostEqual(integral, integralSlow);

end

