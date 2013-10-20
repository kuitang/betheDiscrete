function p = steffensen(f,p0,tol)
% Stolen from http://en.wikipedia.org/wiki/Steffensen%27s_method.
% (Though for this problem you know the derivative so may as well use
% Newton; but maybe this isn't the bottleneck?)
%
% This function takes as inputs: a fixed point iteration function, f, 
% and initial guess to the fixed point, p0, and a tolerance, tol.
% The fixed point iteration function is assumed to be input as an
% inline function. 
% This function will calculate and return the fixed point, p, 
% that makes the expression f(x) = p true to within the desired 
% tolerance, tol. 
 
format compact % This shortens the output.
format long    % This prints more decimal places. 
 
for i=1:1000   % get ready to do a large, but finite, number of iterations.
               % This is so that if the method fails to converge, we won't
               % be stuck in an infinite loop.
    p1=f(p0);  % calculate the next two guesses for the fixed point.
    p2=f(p1);
    p=p0-(p1-p0)^2/(p2-2*p1+p0) % use Aitken's delta squared method to
                                % find a better approximation to p0.
    if abs(p-p0)<tol  % test to see if we are within tolerance.
        break         % if we are, stop the iterations, we have our answer.
    end
    p0=p;              % update p0 for the next iteration.
end
if abs(p-p0)>tol       % If we fail to meet the tolerance, we output a
                       % message of failure.
    error('steffensen:convergence', 'failed to converge in 1000 iterations.')
end
