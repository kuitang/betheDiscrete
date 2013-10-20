% function [b,fval,exitflag] = fzero_fast(FunFcn, x, tol, Uconst, keps, p)
function [b,fval,exitflag] = fzero_fast(FunFcn, x, tol, varargin)

% Initialization
fcount = 0;
iter = 0;
%intervaliter = 0;
exitflag = 1;
%procedure = ' ';

% Interval input
assert(numel(x) == 2, 'fast_fzero requires an initial interval');

a = x(1); savea=a;
b = x(2); saveb=b;

% faSlow = FunFcn(a,varargin{:});
% fbSlow = FunFcn(b,varargin{:});

fa = intNew(a, varargin{:});
fb = intNew(b, varargin{:});

% assert(fa == faSlow);
% assert(fb == fbSlow);

if any(~isfinite([fa fb])) || any(~isreal([fa fb]))
    error(message('MATLAB:fzero:ValuesAtEndPtsComplexOrNotFinite'))
end

savefa = fa; savefb = fb;

if ( fa == 0 )
    b = a;
    fval = fa;
    return
elseif ( fb == 0)
    % b = b;
    fval = fb;
    return
elseif (fa > 0) == (fb > 0)
    error(message('MATLAB:fzero:ValuesAtEndPtsSameSign'))
end
    
fc = fb;

% Main loop, exit from middle of the loop
while fb ~= 0 && a ~= b
    % Insure that b is the best result so far, a is the previous
    % value of b, and c is on the opposite side of the zero from b.
    if (fb > 0) == (fc > 0)
        c = a;  fc = fa;
        d = b - a;  e = d;
    end
    if abs(fc) < abs(fb)
        a = b;    b = c;    c = a;
        fa = fb;  fb = fc;  fc = fa;
    end
    
    % Convergence test and possible exit
    m = 0.5*(c - b);
    toler = 2.0*tol*max(abs(b),1.0);
    if (abs(m) <= toler) || (fb == 0.0) 
        break
    end
    
    % Choose bisection or interpolation
    if (abs(e) < toler) || (abs(fa) <= abs(fb))
        % Bisection
        d = m;  e = m;        
    else
        % Interpolation
        s = fb/fa;
        if (a == c)
            % Linear interpolation
            p = 2.0*m*s;
            q = 1.0 - s;
        else
            % Inverse quadratic interpolation
            q = fa/fc;
            r = fb/fc;
            p = s*(2.0*m*q*(q - r) - (b - a)*(r - 1.0));
            q = (q - 1.0)*(r - 1.0)*(s - 1.0);
        end;
        if p > 0, q = -q; else p = -p; end;
        % Is interpolated point acceptable
        if (2.0*p < 3.0*m*q - abs(toler*q)) && (p < abs(0.5*e*q))
            e = d;  d = p/q;
        else
            d = m;  e = m;            
        end;
    end % Interpolation
    
    % Next point
    a = b;
    fa = fb;
    if abs(d) > toler, b = b + d;
    elseif b > c, b = b - toler;
    else b = b + toler;
    end    
%     fbSlow = FunFcn(b,varargin{:});
    
    fb = intNew(b, varargin{:});
%     assert(fbSlow == fb);
    
    fcount = fcount + 1;
    iter = iter + 1;
end % Main loop

fval = fb; % b is the best value

if abs(fval) <= max(abs(savefa),abs(savefb))
    % noope
else
    exitflag = -5; 
end
