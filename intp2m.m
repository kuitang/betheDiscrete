function [ integral ] = intp2m( Uconst,m,p )
%[ integral ] = intp2m( Uconst,m,p )
%   Returns integral of Uconst + log[q/(1-q)] from m to p
integral = Uconst*(m-p) + xlogx(m) + xlogx(1-m) - xlogx(p) - xlogx(1-p);
end

