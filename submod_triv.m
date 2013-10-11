theta = 10 * ones(3,1);
W     = zeros(3,3);

[energy, xmap, misc] = submodularMAP_mex(theta, sparse(W));
