function kgmm0 = initkgmm_diag(gmm, k,varargin)
%INITKGMM_DIAG initialize kgmm by heuristic.
%
%   INITKGMM_DIAG generates samples from GMM and fit kmeans with K components.
%   The centroids and variance of each cluster are used for initilization
%   of kGMM. Sigma is diagnonal matrix. 1 x NDimensions x K.
%
%   KGMM0 = INITKGMM_DIAG(GMM, K)
%   KGMM0 = INITKGMM_DIAG(GMM, K, NSAMPLES)
%
%   See Also: INITKGMM, RANDOM, KMEANS, GMDISTRIBUTION, OBJ2STRUCTGMM

%   $ Hyunwoo J. Kim $  $ 2014/11/07 01:30:18 (CST) $
 
    gm = gmdistribution(gmm.mu,gmm.Sigma, gmm.PComponents);
    NDimensions = gmm.NDimensions;

    % We will find a good initialization for kGMM.
    nsamples = 1000;
    if length(varargin) > 0 
        nsamples = varargin{1};
    end
    
    samples = random(gm,nsamples);
    [idx, C]  = kmeans(samples,k,'MaxIter',1000);

    mu0 = C;
    Sigma0 = zeros(1, NDimensions, k);
    PComponents0 = ones(1,k)*1/k; % Even weights.

    for i=1:k
        Sigma0(1,:,i) = var(samples(idx==i,:));
    end

    gm0 = gmdistribution(mu0,Sigma0, PComponents0);
    kgmm0 = obj2structGMM(gm0);
end