function val = innerprodGaussian(mu1, mu2, Sigma1, Sigma2)
%INNERPRODGAUSSIAN calculates innerproduct of two Gaussian distribution.
%
%   See Also: L2NORMGMM, INNERPRODGMM, GMDISTRIBUTION, MVNPDF

%   $ Hyunwoo J. Kim $  $ 2014/10/26 19:06:57 (CDT) $

    val = mvnpdf(mu1',mu2',Sigma1+Sigma2);
end