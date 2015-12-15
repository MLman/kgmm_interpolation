function ngmm = l2normalizeGMM(gmm)
%L2NORMALIZEGMM normalize a GMM and devides the PComponents by the norm of
%the GMM.
%
%
%   See Also: L2NORMGMM

%   $ Hyunwoo J. Kim $  $ 2015/12/12 14:17:28 (CLT) $
    ngmm = gmm;
    ngmm.PComponents = gmm.PComponents./l2normGMM(gmm);
end