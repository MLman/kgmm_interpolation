function gmm = randGMM(k, dim)
%RANDGMM generates a random k-GMM in R^dim.
%
%
%   See Also: GMDISTRIBUTION, WISHRND, RANDSPD

%   $ Hyunwoo J. Kim $  $ 2015/09/15 13:21:30 (CDT) $

    mu = rand(k, dim);
    Sigma = zeros(dim,dim,k);
    for i = 1:k
        Sigma(:,:,i) = wishrnd(eye(dim),dim+2);
    end
    p = rand(k,1);
    p = p./sum(p);
    gmm = gmdistribution(mu, Sigma, p);
end