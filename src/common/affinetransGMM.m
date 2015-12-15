function gmmy = affinetransGMM(gmmx, A, varargin)
%AFFINETRANSGMM performs affine transformation of GMM.
%
%   gmmy = affinetransGMM(gmmx, A) % By default, b is zero vector.
%   gmmy = affinetransGMM(gmmx, A, b)
%   Let x ~ N(m, Sigma). y = Ax+b
%   Then y ~ N(b+Am, ASigmaA')
%
%   After affine transformation of GMM, it will get full covariance
%   matrices.
%
%   See Also: GETAFFINETRANSGMM, GETAFFINETRANSGAUSS

%   $ Hyunwoo J. Kim $  $ 2015/04/13 16:36:03 (CDT) $
%   $ Hyunwoo J. Kim $  $ 2015/03/24 15:22:59 (CDT) $

    gmmy = gmmx;
    if nargin < 3
        gmmy.mu = (A*gmmx.mu')'; %b = zeros(size(A,1),1); by default zero.
    else
        b = varargin{1};
        gmmy.mu = (A*gmmx.mu'+repmat(b,1,gmmx.NComponents))';
    end
    
    
    
    % Full covariance matrix.
    if size(gmmx.Sigma,1)==gmmx.NDimensions
        for i =1:gmmx.NComponents
            gmmy.Sigma(:,:,i) = A*gmmy.Sigma(:,:,i)*A';
        end
    else % Diagonal covariance matrix
        assert(size(gmmx.Sigma,1));
        gmmy.Sigma = zeros(gmmx.NDimensions, gmmx.NDimensions, gmmx.NComponents);
        for i =1:gmmx.NComponents
            gmmy.Sigma(1,:,i) = A*diag(gmmy.Sigma(1,:,i))*A';
        end
    end
     
end