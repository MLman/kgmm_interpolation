function gmmbar = l2meanGMMs(gmms, varargin)
%L2MEANGMMS calculates the l2mean of GMMS.
%
%   L2MEANGMMS calculates the l2mean of GMMS with k' = \sum_i^N k_i
%   components, where N is the number of GMMs and k_i is the number of
%   components for each GMM_i.
%
%   GMMS is a array of cells. 
%   Diagonal version needs to be implemented.
%
%      Sigma        An array or a matrix containing the component covariance
%                   matrices.  Sigma is one of the following
%
%                      * A D-by-D-by-K array if there are no restrictions on
%                        the form of covariance.  In this case, Sigma(:,:,J)
%                        is the covariance matrix of component J.
%                      * A 1-by-D-by-K array if the covariance matrices are
%                        restricted to be diagonal, but not restricted to be
%                        same across components.  In this case Sigma(:,:,J)
%                        contains the diagonal elements of the covariance
%                        matrix of component J.
%                      * A D-by-D matrix if the covariance matrices are
%                        restricted to be the same across components, but not
%                        restricted to be diagonal.  In this case, Sigma is
%                        the common covariance matrix.
%                      * A 1-by-D vector if the covariance matrices are
%                        restricted to be diagonal and to be the same across
%                        components.  In this case, Sigma contains the
%                        diagonal elements of the common covariance matrix.
%
% 
%   See Also: GMDISTRIBUTION

%   $ Hyunwoo J. Kim $  $ 2014/11/06 15:55:46 (CST) $
    if ~isempty(varargin)
        W = varargin{1};
        W = W/sum(W);
    end
    
    kk = 0;
    
    N = length(gmms);
    NDimensions = size(gmms{1}.mu,2);
    
    for i =1:N
        kk = kk + gmms{i}.NComponents;
    end
    
    % Initialization.
    gmmbar = gmms{1};%gmmbar.NDimensions = NDimensions; 
    gmmbar.NComponents = kk;
    gmmbar.PComponents = zeros(1, kk);
    gmmbar.mu = zeros(kk, NDimensions);
    
    %gmmbar.Sigma = zeros(NDimensions, NDimensions, kk);
    gmmbar.Sigma = zeros(size(gmms{1}.Sigma));
    % If CovType: 'diagonal'
    % gmmbar.Sigma = zeros(1, NDimensions, kk);
    
    nCompSt = 1;
    for i =1:N
        ki = gmms{i}.NComponents;
        nCompEnd = nCompSt + ki-1;
        if  ~isempty(varargin)
            gmmbar.PComponents(1, nCompSt:nCompEnd) = gmms{i}.PComponents*W(i);
        else
            gmmbar.PComponents(1, nCompSt:nCompEnd) = gmms{i}.PComponents/N;
        end
        gmmbar.mu(nCompSt:nCompEnd, :) = gmms{i}.mu;
        gmmbar.Sigma(:,:, nCompSt:nCompEnd) = gmms{i}.Sigma;
        
        % If CovType: 'diagonal'
        % gmmbar.Sigma(1,:,nCompsSt:nCompEnd) = gmms{i}.Sigma;
        
        nCompSt = nCompEnd + 1;
    end
    %gmmbar.PComponents = gmmbar.PComponents/N;
    
end
