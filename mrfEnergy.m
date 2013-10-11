function [ energy, energyTheta, energyW ] = mrfEnergy(theta, W, x)
% energy = mrfEnergy(theta, W, x) computes energy of a binary MRF.
%   x and W must be column.
    
    energyTheta = -sum(theta .* x);
    energyMat   = bsxfun(@times, bsxfun(@times, x, W), x');
    energyW     = -sum(energyMat(:)) / 2;
    
    energy = energyTheta + energyW;

end

