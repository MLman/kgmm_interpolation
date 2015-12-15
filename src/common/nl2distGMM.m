function d = nl2distGMM(gmm1, gmm2)
%NL2DISTGMM calculates the distance between normalized GMMs.
%
%
%   See Also: L2NORMGMM, INNERPRODGMM, L2DISTGMM

%   $ Hyunwoo J. Kim $  $ 2014/11/05 11:53:09 (CST) $
    p1 = l2normGMM(gmm1);
    gmm1.PComponents = gmm1.PComponents/p1;
    p2 = l2normGMM(gmm2);
    gmm2.PComponents = gmm2.PComponents/p2;
    d = l2distGMM(gmm1, gmm2);
end