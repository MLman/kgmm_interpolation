function val = innerprodGMM(gmm1, gmm2)
%INNERPRODGMM calculates the innerproduct between two GMMs.
%   
%   INNERPRODGMM allows to calculates two Gaussian mixture models with
%   different numbers of components.
%   Assert is removed.
%
%   See Also: L2NORMGMM, GMDISTRIBUTION, MVNPDF, INNERPRODGAUSSIAN

%   $ Version 2$
%   $ Hyunwoo J. Kim $  $ 2014/11/12 19:57:57 (CST) $

    
    k1 = gmm1.NComponents;
    k2 = gmm2.NComponents;
    mus1 = gmm1.mu';
    Sigmas1 = gmm1.Sigma;
    pis1 = gmm1.PComponents; % Weights of components
    mus2 = gmm2.mu';
    Sigmas2 = gmm2.Sigma;
    pis2 = gmm2.PComponents; % Weights of components
    val=0;
    for j = 1:k1
        for jj = 1:k2
%            fprintf('\n Inner product between component %d and %d = %f', j, jj, innerprodGaussian(mus1(:,j), mus2(:,jj), Sigmas1(:,:,j), Sigmas2(:,:,jj)));
            val = val + pis1(j)*pis2(jj)*innerprodGaussian(mus1(:,j), mus2(:,jj), Sigmas1(:,:,j), Sigmas2(:,:,jj));
            %assert(innerprodGaussian(mus1(:,j), mus2(:,jj), Sigmas1(:,:,j), Sigmas2(:,:,jj)) ...
            % == innerprodGaussian(mus2(:,jj), mus1(:,j), Sigmas2(:,:,jj),Sigmas1(:,:,j)));
        end
    end
    if val < 1e-15
        val = 0;
    end
end