function model = gausskde(samples, bandwidth, covtype)
%GAUSSKDE performs kernel density estimation with Guassian kerenl.
%
%   Samples are row vectors.
%
%
%   See Also: GMDISTRIBUTION, KSDENSITY 

%   $ Hyunwoo J. Kim $  $ 2015/03/22 18:29:48 (CDT) $
    
    nsamples = size(samples,1);
    ndim = size(samples,2);
    p = ones(nsamples,1)/nsamples;
    if nargin==2 || strcmpi(covtype,'diag')
        model = gmdistribution(samples, ones(1, ndim, nsamples)*bandwidth,p);
    elseif strcmpi(covtype,'full')
        model = gmdistribution(samples, repmat(eye(ndim)*bandwidth,[1,1,nsamples]),p);
    end
end
