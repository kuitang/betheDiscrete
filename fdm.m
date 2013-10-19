function [ gams,sumN,prodN,thisN ] = fdm(theta, W, A, B, epsilon, method, L, U)
%[ gams,sumN,prodN,thisN ] = fdm(theta, W, A, B, epsilon, method, L, U)
% First derivative mesh revised by AW.
% gams has all the mesh points
% sumN=total no. mesh points across all dims, prodN=\prod_i N_i
% thisN has N_i for each dimension
% methods: simple, minsum, adaptivesimple, adaptiveminsum (see SODA paper)

%   calculates the first-derivative mesh
%   to accuracy epsilon. gams is a N-vector of step sizes, one in each
%   variable. A, B, L, U are N-vectors of bounds; L and U come from BBP.
%   If not provided, they're calculated here based on the A, B vectors of
%   bounds provided.
%   Define \gamma as the *width* of the mesh (so that, e.g.
%   q_{i,n+1} = q_{i,n} + \gamma_i). (This is at variance with notes2,
%   where \gamma is the distance to the nearest mesh point. So this
%   function will sometimes use 2*\gamma).
%
%   This means first mesh point is at Ai+gamma_i/2 then keep
%   adding gamma_i to get the next point, until we reach >= Bi
%
%   method can be simple, lagrangian, adaptive. If adaptive, gams will be a
%   length-N cell array, and will contain an array of grid points instead
%   of just step sizes.
%
%   Calculate the maximum first derivative in each direction. The bound
%   is
%
%   -\theta_i - W_i + \log U_i + \log \frac{q_i}{1 - q_i}
%   \leq \frac{\partial F}{\partial q_i}
%   \leq -\theta_i + V_i - \log L_i + \log \frac{q_i}{1 - q_i}
%
%   To bound, we consider the maximum modulus of the bound on either side.
%   This amounts to checking the lower-left and upper-right corners.
%
%   Later we extend     - integrate rather than sum max
%                       - factor in direction, so e.g. if start from left
%                       Ai, the first mesh point should cover all possible
%                       global min q^ to the left, so only need use
%                       integral of upper bound of derivative (and can
%                       ignore the lower bound function). Then we compute
%                       the 'reach to the right' of the first mesh point,
%                       where for this we must use the lower bound
%                       derivative, which is used for the upper bound of
%                       the - derivative.  From there, again we use the
%                       upper bound only to compute the next mesh point,
%                       which must reach back to the left to the previous
%                       reach point, and continue...
%
%   Note the integral of \log \frac{q}{1 - q} is
%       q \log q + (1-q) \log (1-q) + const
%
%   We can raise the lower bound curve by \log U_i and lower the upper
%   bound curve by \log L_i.

n = length(theta);

Wpos =  sum(W .* (W > 0), 2); % W_i in notes
Wneg = -sum(W .* (W < 0), 2); % V_i in notes

% Check if L and U provided, if not then compute in one pass of BBP using
% A, B provided. Note the L and U obtained are the best we can easily
% obtain. If we used these to generate new A and B, they might be worse
% than the ones we were given. We assume here that we won't improve A and B
% since we assume something at least as good as BBP has already been run to
% completion to produce the A and B which were input.
if nargin<8    
    % Do one half iteration of BBP
    [~, ~, alpha, L, U] = BBPNew(theta, W, 'A', A, 'B', B, 'maxIter', 1);    
end
    
% Constants for the bound functions; exclude \log \frac{q_i}{1 - q_i}
lb   = -theta - Wpos + log(U);
ub   = -theta + Wneg - log(L);

ll   = lb + log(A) - log(1 - A);
ur   = ub + log(1 - B) - log(B);
D    = max(abs(ll), abs(ur));

S      = 1 - B - A; % gap sizes

assert(all(S >= -10 * eps), 'fudge!');
S(S < 0) = 0;

% Now calculate the gamma mesh.
if strcmp(method, 'simple')
    gams   = 2 * epsilon / n * (1./D);  % epsilon and n are scalars, D is a vector
elseif strcmp(method, 'minsum')
    sqrtSoverD = sqrt(S ./ D); sumsqrtSD=sum(sqrt(S .* D));
    
    sqrtSoverD(S == 0) = 0;
    
    gams   = 2 * epsilon * sqrtSoverD / sumsqrtSD;
elseif strcmp(method(1:8), 'adaptive')
    % This uses the whole bag of tricks. So adaptive, integrate,
    % directional. Though for the moment we keep the deltaF_i = epsilon/n
    % which can probably be improved similarly to lagrangian, need to think
    % For now I always go left->right placing points as extreme as
    % possible along the way. May need to intro tolerance for numerical
    % error *
    
    % lower prod. lowest?
    if strcmp(method(9:end), 'simple')
        k = 1/n * ones(n,1); keps = k*epsilon; % we'll make the max delta F_i = k epsilon 
    elseif strcmp(method(9:end), 'minsum')
        sqrtSD=sqrt(S .* D); sumsqrtSD=sum(sqrtSD); keps=epsilon/sumsqrtSD *sqrtSD; 
    else
        fprintf('adaptive but what type?\n'); beep; return
    end
    gams = cell(n,1); sumN=n; % N is total no. of mesh points in all dims
    thisN=ones(n,1);
    for i=1:n
        %fprintf('Starting computation of mesh points gams{i} for i=%d\n',i);
        Uconst=ub(i); Lconst=lb(i); prevr=A(i);
        gams{i}=[prevr]; % this is just to make it a vec, take it off at the end (unless A=1-B in which case leave it as is]
        if A(i)<1-B(i)
            while prevr<1-B(i)
%                 [m,nextr]=adapt(prevr,Uconst,Lconst,keps(i),A(i),B(i));
                [m,nextr]=adaptRobust(prevr,Uconst,Lconst,keps(i),A(i),B(i));

                gams{i}=[gams{i} m];
                prevr=nextr;
                sumN=sumN+1; thisN(i)=thisN(i)+1;
            end
            gams{i}=gams{i}(2:end); % strips off the initial point which was needed to make sure it would work as a vec
            sumN=sumN-1; thisN(i)=thisN(i)-1;
        end
    end  
end

if strcmp(method, 'simple') || strcmp(method, 'minsum') % calc N
    frac=(S - gams) ./ gams;
    frac(S == 0) = 0;
    thisN=round(frac+0.51)+1; % no. mesh points in each dimension
    sumN=sum(thisN);
end
prodN=prod(thisN);   

end

