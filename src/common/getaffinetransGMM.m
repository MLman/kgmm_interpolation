function T = getaffinetransGMM(gmm1, data1, gmm2, data2, varargin)
%GETAFFINETRANSGMM returns a mixture of affine transformation. The weight
%of affine transform is calculated by responsibility in GMM1. Data should
%be in GMM1 space.
%
%
%
%   See Also: GETAFFINETRANSGAUSS, AFFINETRANSGMM, MOAtransform

%   $ Hyunwoo J. Kim $  $ 2015/04/01 22:35:52 (CDT) $

    assert(gmm1.NComponents == gmm2.NComponents,'The Number of components is NOT same.');
    assert(isobject(gmm1) && isobject(gmm2),'Two input gmms MUST be objects');
    
    T.NComponents = gmm1.NComponents;
    K = T.NComponents;
    
    if isstruct(data1) 
        X1 = data1.x;
    else
        X1 = data1;
    end
    if isstruct(data2) 
        X2 = data2.x;
    else
        X2 = data2;
    end
    if isempty(varargin)
        mode = [];
    else
        mode = varargin{1};
    end
    T.As = zeros(size(X1,2),size(X1,2),K);
    T.bs = zeros(size(X1,2),K);
    
    P1 = posterior(gmm1,X1);
    P2 = posterior(gmm2,X2);    
    
    d1.x = X1;
    d2.x = X2;
    rX1 = cluster(gmm1, X1);
    rX2 = cluster(gmm2, X2);
	for i = 1:K
        g1.mu = gmm1.mu(i,:);
        g1.Sigma = gmm1.Sigma(:,:,i);
        g2.mu = gmm2.mu(i,:);
        g2.Sigma = gmm2.Sigma(:,:,i);

       d1.w = P1(:,i);
       d2.w = P2(:,i);
%        d1.w = zeros(size(d1.x,1),1);
%        d1.w(rX1==i) = 1;
%        d2.w = zeros(size(d2.x,1),1);
%        d2.w(rX2==i) = 1;

%        d1.w = ones(size(d1.x,1),1);
%        d2.w = ones(size(d2.x,1),1);


        [A, b] = getaffinetransGauss(g1, d1, g2, d2, mode);
        
        T.As(:,:,i) = A;
        T.bs(:,i) = b;
    end
    T.gmm = gmm1;
    
end