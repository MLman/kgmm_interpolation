function [gmm1, stats] = em_gmm_closest_full(gmm2, K, varargin)
%EM_GMM_CLOSEST_FULL performs EMstyle GMM fitting with K components to
%GMM2.
%
%   EM_GMM_CLOSEST_FULL(GMM2, K)
%   EM_GMM_CLOSEST_FULL(GMM2, K, OPTION)
%
%   OPTION.maxiter = 100 default value.
%   OPTION.tol = 1e-6.
%   OPTION.debug = 1 returns detailed stats.
%   OPTION.debug = 2 show the intermediate fittings.
%   OPTION.getGamma = 'getGamma_xentropy', (default) responsibilities is evaluated by
%   cross entproy.
%   OPTION.getGamma = 'getGamma_inprod', responsibilities is evaluated by
%   innerproduct.
%   See Also : GD_GMMS_CLOSEST_FAST

%   $ Hyunwoo J. Kim $  $ 2015/04/04 18:25:56 (CDT) $
    
    %% Random initialization
    % pick k Gaussians from GMM2

    
    idx = randsample(1:gmm2.NComponents,K);
    %% Debugging
    


    % Full covariance matrix.
    assert(covtype(gmm2.Sigma) == 1, 'Covariance matrix MUST a full matrix.');
    if nargin <= 3
        gmm1 = gmdistribution(gmm2.mu(idx,:), gmm2.Sigma(:,:,idx), ones(1,K)/K);
        gmm1 = obj2structGMM(gmm1);
    else
        gmm1 = varargin{2};
        if ~isstruct(gmm1);
            gmm1 = obj2structGMM(gmm1);
        end
    end
    
    
    N = gmm2.NComponents;
    dims = size(gmm2.mu,2);
    
    if nargin >= 3
        option = varargin{1};
    else 
        option = [];
    end
    if isfield(option,'maxiter')
        maxiter = option.maxiter;
    else 
        maxiter = 300; 
    end

    if isfield(option, 'debug') && (option.debug >= 1)
        debug = option.debug;
    else 
        debug = 0;
    end
    
    if isfield(option, 'tol')
        tol = option.tol;
    else
        tol = 1e-6;
    end
    
	if isfield(option, 'dflag')
        dflag = option.dflag;
	else
        dflag = ones(1,3);
	end
    
    if isfield(option,'getGamma') && strcmpi(option.getGamma,'getGamma_inprod')
        getGamma = @getGamma_inprod;
    else
        getGamma = @getGamma_xentropy; 
    end
    
    if isfield(option, 'monotonicity')
        monotonicity = option.monotonicity;
    else
        monotonicity = false;
    end
    
    l2history = zeros(maxiter,1);
    stats =[];
    stats.terminate = 'Not converged.';
    diffnorm = 10000; % Inf
    
    %% debug = 1;
    
    if length(idx) > 1
        %gmm2.Sigma(:,:,idx);
        debug = 1;
    end
    for iter = 1:maxiter
        %iter
        if debug >= 1
            l2history(iter) = l2distGMM(gmm1,gmm2);
            fprintf('L2 distance %f, parameter change %e tol %e \n', l2history(iter), diffnorm, tol);
            %gmm1.Sigma;
            %gmm1.PComponents;
        end

        %%% E-step %%%
        %% Calculate responsibilies.
        gammas =  getGamma(gmm1, gmm2);
        
        gmm1new = gmm1;
        %%% M-step %%%
        %% Calculate Mean and Variance of each components.
        % Mean update (Weighted mean of means of all Gaussians)
        % mubar = sum mu*gamma / sum gamma
        
        W2 = gmm2.PComponents(:); % Weight of gmm2. Data Gaussian mixture models.
        if dflag(2) == 1
            % mu update 
            for j = 1:K
                W2J = W2.*gammas(:,j);
                gmm1new.mu(j,:) = sum(repmat(W2J,1,dims).*gmm2.mu,1)/sum(W2J);
            end
        end
    
        % Sigma udpate
        if dflag(3) == 1
            for j = 1:K
                % sum alpha_i C_i + sum alpha_i (mu_i-mubar)(mu_i-mubar)'
                Sigma = zeros(dims);
                W2J = W2.*gammas(:,j);
                for i = 1:N
                    Sigma = Sigma + gmm2.Sigma(:,:,i)*W2J(i);
                end
                mSigma = Sigma/sum(W2J);
                M = gmm2.mu - repmat(gmm1new.mu(j,:), N, 1);
                gmm1new.Sigma(:,:,j) = mSigma + M'*diag(W2J)*M/sum(W2J);
            end
        end
        
        % pi update (PComponents)
        if dflag(1) == 1
            pi2gammas = repmat(W2,1,K).*gammas;
            w = sum(pi2gammas, 1);
            gmm1new.PComponents = w/sum(w);
        end
        

        % Stopping condition
        diffmu = gmm1new.mu - gmm1.mu;
        diffS = gmm1new.Sigma - gmm1.Sigma;
        diffpi = gmm1new.PComponents - gmm1.PComponents;
        gmm1 = gmm1new;
        diffnorm = norm(diffmu)+norm(diffS(:))+norm(diffpi);
        if diffnorm < tol 
            stats.terminate = 'Converged.';
            break;
        end
        
        if monotonicity && iter >1 && (l2history(iter-1) < l2history(iter))
            stats.terminate = 'Coverged.';
            break;
        end
        
        if mod(iter,10) == 1 && debug >= 2
            gmmhat = struct2objGMM(gmm1);
            figure
            hold on
            ezcontourf(@(x,y)log(pdf(gmmhat,[x y])),[-0.3 1],[-0.3 1]);
            plot(gmm2.mu(:,1), gmm2.mu(:,2), 'b.');
            title(sprintf('K means style GMM fitting iter: %d.',iter))
            hold off
        end        
    end
    
    stats.iter = iter;
    if debug
        stats.l2history = l2history(1:iter);
    end
end

function gammas = getGamma_inprod(gmm1, gmm2)
    K = gmm1.NComponents;
    N = gmm2.NComponents;
    gammas = zeros(N,K);
    for i = 1:N
        for j = 1:K 
            gammas(i,j) = innerprodGaussian(gmm2.mu(i,:), gmm1.mu(j,:), gmm2.Sigma(:,:,i), gmm1.Sigma(:,:,j));
        end
    end
    gammas = repmat(gmm1.PComponents,N,1).*gammas;
    gammas = gammas./repmat(sum(gammas,2),1,K);
end

%% GMM1 is variable. K-GMM.
%% GMM2 is data. N-GMM.
function gammas = getGamma_xentropy(gmm1, gmm2)
    K = gmm1.NComponents;
    N = gmm2.NComponents;
    gammas = zeros(N,K);

    for i = 1:N
        for j = 1:K 
            gammas(i,j) = exp(-crossentropy(gmm2.mu(i,:)', gmm2.Sigma(:,:,i), gmm1.mu(j,:)', gmm1.Sigma(:,:,j)));
        end
    end

    gammas = repmat(gmm1.PComponents,N,1).*gammas;
    gammas = gammas./repmat(sum(gammas,2),1,K);
end

