function d = l2distGaussian(mu1, Sigma1, mu2, Sigma2)
%L2DISTGMM returns L2 distance between two Gaussian distributions.
%
%
%   See Also: INNERPRODGAUSSIAN, L2NORMGMM, INNERPRODGMM, GMDISTRIBUTION, MVNPDF

%   $ Hyunwoo J. Kim $  $ 2015/04/17 13:43:14 (CDT) $
    if norm(mu1-mu2) + norm(Sigma1-Sigma2) < 1e-8
        d = 0 ;
        return 
    end
    d = sqrt(innerprodGaussian(mu1, mu1, Sigma1, Sigma1)...
        -2*innerprodGaussian(mu1, mu2, Sigma1, Sigma2) ...
        + innerprodGaussian(mu2, mu2, Sigma2, Sigma2));
	try
        assert(isreal(d))
    catch
        d
        disp(sprintf('Numerial problem in %s', mfilename));
        d =real(d);
    end
end