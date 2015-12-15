function d = l2distGMM(gmm1, gmm2)
%L2DISTGMM calculates the distance between unnormalized GMMs.
%
%
%
%   See Also: L2NORMGMM, INNERPRODGMM, L2NDISTGMM

%   TODO gmm1 is full and gmm2 diag or vice versa.
%   TODO ComputeInnderProductGMM is not defined.
%   $ Hyunwoo J. Kim $  $ 2015/03/25 02:15:48 (CDT) $

    % Diagonal case
    if size(gmm1.Sigma,1) == 1 && size(gmm2.Sigma,1) == 1 && ismac ==0
        d2 = ComputeInnderProductGMM(gmm1,gmm1) +  ComputeInnderProductGMM(gmm2,gmm2) ...
        - 2*ComputeInnderProductGMM(gmm1, gmm2);
    
    % Numerical problem. 
        if d2 < 1e-20 || d2 < 0
            d = 0;
            return 
        end
        
    % Full Covariance case    
    elseif (size(gmm1.Sigma,1) > 1 && size(gmm1.Sigma,1) == size(gmm2.Sigma,1)) || ismac || ispc
        d2 = innerprodGMM(gmm1, gmm1)+ innerprodGMM(gmm2, gmm2) - 2*innerprodGMM(gmm1, gmm2);
    end
    d = sqrt(d2);
end