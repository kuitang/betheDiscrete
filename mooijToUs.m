function [theta, W] = mooijToUs(eta, J)
% [theta, W] = mooijToUs(eta, J) change parameters from Mooij and Kappen
%
%   Their theta is eta here. theta is our theta.
    
    J = triu(J, 1);
    J = J + J';

    theta = 2*eta - 2*sum(J, 2);
    W = 4 * J;    

end

