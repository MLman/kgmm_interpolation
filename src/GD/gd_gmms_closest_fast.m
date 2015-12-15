function [kgmm, fval, gnorm, fval_history, status, gmm1init] = gd_gmms_closest_fast(gmm2, k, maxiter,varargin)
%GD_GMMS_CLOSEST performs the gradient descent method to find the closest KGMM to GMM2.
%
%   [KGMM, FVAL, GNORM, FVAL_HISTORY] = GD_GMMS_CLOSEST(GMM2, K, MAXITER,VARARGIN)
%
%   See Also: L2DISTGMMS, L2MEANGMM, DFDPI_CLOSEST, DFDM_CLOSEST, DFDS_CLOSEST

%   $ Hyunwoo J. Kim $  $ 2014/11/07 00:02:24 (CST) $


% Initialization is very important.

    verbose = false;
    maxalph = 10000;
    minalph = 1e-15;
    nsamples = 100000;
    MAXGNORM = 1e10;
    % Initialization by the first element.
    % There is much better initialization.
    % For example, convex combination of gmm0, gmmN and project to the feasible
    % set. This requires inner product between hetregeneous GMMs with different
    % number of components.

    if length(varargin) >= 1
        flags = varargin{1}; % flags [dfdpi dfdm dfdS]
    else
        flags = ones(1,3);
    end
    
    %init 3
    if length(varargin) >= 2
        gmm1 = varargin{2}; % Initialization outside
    else
        % Random initialization.
        gmm1 = initkgmm_diag(gmm2, k, nsamples);
    end
    
   
    gmm1init = gmm1; % For debuggin
    fval_history = zeros(maxiter,1);
    alph = 1;
    isdone =0;
    updated = true;
    for iter = 1:maxiter
         if mod(iter,10) == 1 && alph < 1
             alph = 1;
         end
        
        fval = l2distGMM(gmm1,gmm2);
        fval_history(iter) = fval;

        % Calculate gradient
        [grad_pi, grad_m, grad_s] = df_closest_fast(gmm1, gmm2);
        % Optimality condition
        gnorm2 = norm([flags(1)*grad_pi(:);flags(2)*grad_m(:);flags(3)*grad_s(:);]);
        gnorm1 = norm([flags(2)*grad_m(:);flags(3)*grad_s(:);]);

        if verbose
            fprintf('[%d][%s] fval  %e, gnorm1 : %e, gnorm2 %e, mu %f, Sigma %f, log10(alpha) %f\n', iter,question(updated,'up','no'), fval_history(iter), gnorm1, gnorm2, gmm1.mu(1), gmm1.Sigma(1), log10(alph));
        end
        if gnorm1*alph < 1.0000e-16 || ~updated || gnorm1  < 1.0000e-10
            isdone =1;
            break;
        end
        if gnorm2 > MAXGNORM
            grad_pi = grad_pi/gnorm2;
            grad_m = grad_m/gnorm2;
            grad_s = grad_s/gnorm2;
        end
        
        gmm1_bak = gmm1;
        while ~isdone
            updated = false;
            % dfdpi
            if flags(1)
                gmm1.PComponents = gmm1.PComponents - grad_pi*alph;
                % Projection, Nonegativity.
                gmm1.PComponents = max(gmm1.PComponents, 0); 
                % Projection, Sum is 1.
                %gmm1.PComponents = gmm1.PComponents/sum(gmm1.PComponents); 
            end

            % dfdm
            if flags(2)
                gmm1.mu = gmm1.mu - grad_m*alph;
            end
            % dfds
            if flags(3)
                gmm1.Sigma = gmm1.Sigma - grad_s*alph;  
                gmm1.Sigma = max(gmm1.Sigma, eps); 
            end
      
            fval_new = l2distGMM(gmm1, gmm2);
            if isnan(fval_new) && alph >= minalph
                gmm1 = gmm1_bak;
                alph = alph / 2;
                continue
            end
            if alph < minalph 
                gmm1 = gmm1_bak;
                isdone = 1;
                break
            end
            if l2distGMM(gmm1, gmm2) > fval
                alph = alph / 2;
                gmm1 = gmm1_bak;
            else
                alph = min(alph * 2, maxalph);
                updated = true;
                break
            end
        end
    end
    
    if isdone ==1;
        fval_history(1:iter);
        status = 'Converged.';
    else
        status = 'Maxiter';
    end
    if verbose
        disp('Finished.')
    end
    gnorm = [gnorm1, gnorm2];
    % Projection, Sum is 1.
    gmm1.PComponents = gmm1.PComponents/sum(gmm1.PComponents); 
    fval = l2distGMM(gmm1,gmm2);
    fval_history = [fval_history;fval];
    kgmm = gmm1;
end

function val = question(statement, val1, val2)
    if statement
        val = val1;
    else
        val = val2;
    end
end