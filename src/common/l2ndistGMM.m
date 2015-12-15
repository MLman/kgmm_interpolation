function d = l2ndistGMM(gmm1, gmm2)
%L2NDISTGMM calculates the distance between unnormalized GMMs.
%
%
%
%   See Also: L2NORMGMM, INNERPRODGMM, L2DISTGMM

%   TODO gmm1 is full and gmm2 diag or vice versa.
%   $ Hyunwoo J. Kim $  $ 2015/03/25 02:15:48 (CDT) $

    ngmm1 = l2normalizeGMM(gmm1);
    ngmm2 = l2normalizeGMM(gmm2);
    d = l2distGMM(ngmm1, ngmm2);
end