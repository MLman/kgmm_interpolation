function v = l2normGMM(GMM)
%L2NORMGMM calculates L2 norm of Gaussian Mixture Model.
%
%   See Also: INNDERPRODUCTGMM, L2DISTGMMS, L2DISTGMM
% TODO : integration with mex files for speed up.

%   $ Hyunwoo J. Kim $  $ 2015/09/15 13:40:35 (CDT) $
%   $ Hyunwoo J. Kim $  $ 2015/03/05 23:39:35 (CST) $

%    v = sqrt(ComputeInnderProductGMM(GMM,GMM));
    v = sqrt(innerprodGMM(GMM,GMM));
end
