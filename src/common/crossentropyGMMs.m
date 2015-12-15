function H = crossentropyGMMs(gmm1, gmm2)
%CROSSENTROPYGMMS
%
%
%   See Also:

%   $ Hyunwoo J. Kim $  $ 2015/10/13 14:21:56 (CDT) $

    %% Calculate responsibilies.
    [gammas, xentropy] =  getGamma_xentropy(gmm1, gmm2);
    try
        H = sum(gmm1.ComponentProportion(:).*sum(gammas.* xentropy,2));
        %H = H + sum(sum(repmat(gmm1.ComponentProportion',1,gmm2.NComponents).*log(gammas)));
        %H = H - sum(gmm2.ComponentProportion(:).*log(gmm2.ComponentProportion(:)));
    catch
        H = sum(gmm1.PComponents(:).*sum(gammas.* xentropy,2));
        %H = H + sum(sum(repmat(gmm1.PComponents',1,gmm2.NComponents).*log(gammas)));
        %H = H - sum(gmm2.PComponents(:).*log(gmm2.PComponents(:)));
    end
end

%% GMM2 is data. N-GMM.
function [gammas, xentropy] = getGamma_xentropy(gmmi, gmmj)
% gmm2 is data.

    Ki = gmmi.NComponents;
    Kj = gmmj.NComponents;
    
    gammas = zeros(Ki,Kj);
    xentropy = zeros(Ki,Kj);
    
    for i = 1:Ki
        for j = 1:Kj
            xentropy(i, j) = crossentropy(gmmi.mu(i,:)', gmmi.Sigma(:,:,i), gmmj.mu(j,:)', gmmj.Sigma(:,:,j));
            gammas(i,j) = exp(-xentropy(i, j));
        end
    end

    gammas = repmat(gmmj.PComponents,Ki,1).*gammas;
    gammas = gammas./repmat(sum(gammas,2),1,Kj);
end