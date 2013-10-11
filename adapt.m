function [ m,nextr ] = adapt( prevr,Uconst,Lconst,keps,A,B,TOL )
%[ m,nextr ] = adapt( prevr,Uconst,Lconst,keps,A,B,TOL )
%   Returns next mesh point location and its reach to the right
%   Inputs are previous r location, const for Upper bound (should be -theta + V - log L), const for Lower
%   bound (should be -theta -W + log U); we return s.t. within keps=k.epsilon
%   Assumes we're working from left (low) to right (high)
%   With left bound A, right boung 1-B

% For now use loose upper bound, i.e. assume integral is max val * interval
% No seems not much harder to try integral so we'll do that
if nargin<7; TOL=0.99; end % We'll find m,nextr s.t. integrals are in range [LOTWOL,1] * keps

fudge = 1e-6;

C=Uconst-Lconst;
if C<0
    error('Error: Upper bound constant < Lower bound constant\n');    
end
if A<0 || (1-B)>1 || A>(1-B)
    error('Error: invalid A, 1-B bounds\n');    
end
if prevr<A || prevr>1-B
    error('Error: Previous r not in [A,1-B]\n');    
end
if Uconst+log(prevr/(1-prevr)) < 0 - fudge
    error('Error: Upper bound<0 at prevr\n');    
end
if Lconst+log((1-B)/B) > 0 + fudge
    error('Error: Lower bound>0 at right edge 1-B\n');    
end

p=prevr;
% Want to find m s.t. upper bd U(m) * (m-prevr) = keps
% No we'll try integral so we want m s.t. \int_p^m U(q) dq = keps
%   i.e. [Uc q + q log q + (1-q) log (1-q)]_p^m = keps
m=1-B; nextr=1-B;
if intp2m(Uconst,1-B,prevr)>keps % if <= then already done
    iter=0;
    curint=intp2m(Uconst,m,p);
    curfrac=curint/keps;
    bestlo=p; besthi=1-B;
    while (curfrac<TOL) || (curfrac>1)
        % Shall we use binary search or use curfrac? or use deriv at m?
        % Let's use deriv if poss else binary?
        % Actually here use binary down and deriv up (note deriv is
        % increasing monotonically so this undermoves up)
        iter=iter+1;
        if curfrac>1
            if m<besthi; besthi=m; end
            %m=p+(mlast-p)/2;
            m=p+(besthi-p)/(curfrac+0.01); % 0.01 to prevent numerical problems
%            if m<=bestlo || m>=besthi; m=(bestlo+besthi)/2; disp('Interpolated m'); end
            if m<=bestlo || m>=besthi; m=(bestlo+besthi)/2; end            
        else
            if m>bestlo; bestlo=m; end
            derivAtM=Uconst+log(m/(1-m));
            deltam=(keps-curint)/derivAtM *0.99;  % 0.99 is to prevent numerical overshoot
            m=m+deltam;
%            if m<=bestlo || m>=besthi; m=(bestlo+besthi)/2; disp('Interpolated m'); end
            if m<=bestlo || m>=besthi; m=(bestlo+besthi)/2; end            
        end
        curint=intp2m(Uconst,m,p);
        curfrac=curint/keps;
        %fprintf('m iteration #%d, m=%8.6f, curfrac=%8.6f\n',iter,m,curfrac);
    end
    %fprintf('Found m at %8.6f after %d iterations. curfrac=%8.6f\n',m,iter,curfrac);
    
    iter=0;
    r=1-B; bestlo=m; besthi=1-B;
    curint=intp2m(Lconst,m,r); % note this will be positive *
    if curint>keps % if <= then already done
        curfrac=curint/keps;
        while (curfrac<TOL) || (curfrac>1)
            iter=iter+1;
            if curfrac>1
                if r<besthi; besthi=r; end
                rlast=r; r=m+(rlast-m)/(curfrac+.01);
%                if r<=bestlo || r>=besthi; r=(bestlo+besthi)/2; disp('Interpolated r'); end
                if r<=bestlo || r>=besthi; r=(bestlo+besthi)/2; end               
            else
                if r>bestlo; bestlo=r; end
                derivAtR=-(Lconst+log(r/(1-r))); % positive
                deltar=(keps-curint)/derivAtR *0.99;
                r=r+deltar;
%                if r<=bestlo || r>=besthi; r=(bestlo+besthi)/2; disp('Interpolated r'); end
                if r<=bestlo || r>=besthi; r=(bestlo+besthi)/2; end               
            end
            curint=intp2m(Lconst,m,r); curfrac=curint/keps;
            %fprintf('r iteration #%d, r=%8.6f, curfrac=%8.6f\n',iter,r,curfrac);
        end
    end
    %fprintf('Found nextr at %8.6f after %d iterations. curfrac=%8.6f\n',r,iter,curfrac);
    nextr=r;
end

end

