function [J, eta] = usToMooij(W, theta)
% [J, eta] = usToMooij(W, theta) change parameters from us to Mooij

    J = 0.25 * W;
    eta = 0.5*theta + sum(J, 2);
    
end

