function d = ngeodistGMM(gmm1, gmm2)
%NGEODISTGMM computes geodesic distance between normalized GMMs.
%
%
%   See Also: NL2DISTGMM, L2NORMGMM, INNERPRODGMM, L2DISTGMM

%   $ Hyunwoo J. Kim $  $ 2015/09/15 13:53:26 (CDT) $
    p1 = l2normGMM(gmm1);
    gmm1.PComponents = gmm1.PComponents/p1;
    p2 = l2normGMM(gmm2);
    gmm2.PComponents = gmm2.PComponents/p2;
    d = acos(innerprodGMM(gmm1, gmm2));
end