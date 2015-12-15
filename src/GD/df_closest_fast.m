function [dpi, dm, dS] = df_closest_fast(gmm1, gmm2)
%DF_CLOSEST_FAST returns the first derivative w.r.t pi1, m1, and S1.
%
%
%   See Also: DFDPI_CLOSEST, DFDM_CLOSEST, DFDS_CLOSEST


%   $ Hyunwoo J. Kim $  $ 2015/12/13 01:23:52 (CLT) $ 
%                       $   Comments are removed.   $
%   $ Hyunwoo J. Kim $  $ 2014/11/08 19:31:10 (CST) $

    K1 = gmm1.NComponents;
    K2 = gmm2.NComponents;
    NDimensions = gmm1.NDimensions;
    
    % Access M_{1,2}^{j,j'} = M11(j,:,j') %Row vector is returned.
    M11  = repmat(gmm1.mu,[1 1 K1]) - repmat(permute(gmm1.mu,[3 2 1]),[K1 1 1]);
    M12  = repmat(gmm1.mu,[1 1 K2]) - repmat(permute(gmm2.mu,[3 2 1]),[K1 1 1]);
    
    S11 = repmat(permute(gmm1.Sigma,[3,2,1]),[1 1 K1]) +repmat(gmm1.Sigma,[K1,1,1]);
    S12 = repmat(permute(gmm1.Sigma,[3,2,1]),[1 1 K2]) +repmat(gmm2.Sigma,[K1,1,1]);
       
    U11  = M11./S11;
    U12  = M12./S12;
    
    C11 = zeros(K1,K1);
    C12 = zeros(K1,K2);
    myzeros = zeros(K1,NDimensions);
    for jj = 1:K1
        C11(:,jj) = mvnpdf(M11(:,:,jj), myzeros, permute(S11(:,:,jj),[3 2 1]));
    end
    for jj = 1:K2
        C12(:,jj) = mvnpdf(M12(:,:,jj), myzeros, permute(S12(:,:,jj),[3 2 1]));
    end
        
    %%
    P1 = gmm1.PComponents;
    P1 = repmat(P1,[K1 1]);
    P2 = gmm2.PComponents;
    P2 = repmat(P2,[K1 1]);
    Z11 = P1.*C11;
    Z12 = P2.*C12;
   
    W11 = permute(repmat(Z11,[1,1,NDimensions]),[1, 3,2]).*U11;
    W12 = permute(repmat(Z12,[1,1,NDimensions]),[1, 3,2]).*U12;
        
    dpi = 2*(sum(Z11,2)-sum(Z12,2))';
    P1 = repmat(gmm1.PComponents',[1 NDimensions]);
    dm = -2*P1.*squeeze(sum(W11,3) - sum(W12,3));
    %%
    dS = -P1.*(sum(permute(repmat(Z11,[1 1 NDimensions]),[1 3 2]).*(1./S11-U11.^2),3)...
        -sum(permute(repmat(Z12,[1 1 NDimensions]),[1 3 2]).*(1./S12-U12.^2),3));
    dS = permute(dS,[3 2 1]);
end