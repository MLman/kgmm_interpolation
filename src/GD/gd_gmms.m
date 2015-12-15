function [gmms, fval, gnorm, fval_history] = gd_gmms(gmm0, gmmN, npts, maxiter,varargin)
%GD_GMMS performs the gradient descent method.
%
%   [GMMS, FVAL, GNORM, FVAL_HISTORY] = GD_GMMS(GMM0, GMMN, NPTS, MAXITER,VARARGIN)
%   GD_GMMS(GMM0, GMMN, NPTS, MAXITER,C1,FLAGS, GMMSINIT)
%
%   See Also: L2DISTGMMS, DFDPI, DFDM, DFDS

%   $ Hyunwoo J. Kim $  $ 2014/11/02 00:02:02 (CDT) $

% Initialization
global verbose maxalph minalph C1

    N = npts+2;
    gmms = cell(1,N);
    gmms{1} = gmm0;
    gmms{N} = gmmN;

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
    ndim = gmm0.NDimensions;
    K = gmm0.NComponents;

    if length(varargin) >= 2
        gmms = varargin{2}; % Initialization outside
    else
        for i = 2:N-1
            gmms{i} = gmmN;
            gmms{i}.PComponents = ones(1,K)/K;
            gmms{i}.mu = C1*rand(size(gmms{i}.mu));
            gmms{i}.Sigma = repmat(eye(ndim),[1 1 K]);
        end
    end

    fval_history = zeros(maxiter,1);
    alph = 1;
    isdone =0;
    %minalph = 1e-20;
    %maxalph = 1000;
    updated = false;
    for iter = 1:maxiter
        if mod(iter,10) == 1 && alph < 1
            alph = 1;
        end

        fval = l2distGMMs(gmms)/npts;
        fval_history(iter) = fval;

        % Calculate gradient
        grad_pi = dfdpi(gmms)/npts;
        grad_m = dfdm(gmms)/npts;
        grad_s = dfds(gmms)/npts;
        % Optimality condition
        gnorm2 = norm([flags(1)*grad_pi(:);flags(2)*grad_m(:);flags(3)*grad_s(:);]);
        gnorm1 = norm([flags(2)*grad_m(:);flags(3)*grad_s(:);]);

        if verbose
            %fprintf('[%d][%s] fval  %f, gnorm1 : %f, gnorm2 %f, mu %f, Sigma %f, alpha %f\n', iter,question(updated,'up','no'), fval_history(iter), gnorm1, gnorm2, gmms{2}.mu, gmms{2}.Sigma, alph);
            fprintf('[%d][%s] fval  %f, gnorm1 : %f, gnorm2 %f, mu %f, Sigma %f, alpha %f\n', iter,question(updated,'up','no'), fval_history(iter), gnorm1, gnorm2, gmms{2}.mu(1), gmms{2}.Sigma(1), alph);
        end
        if gnorm1 < 1e-15
            isdone =1;
            break;
        end
        
        gmms_bak = gmms;
%          if iter > 20
%              disp('Check this out.');
%          end
        while true
            updated = false;
            for i = 2:N-1
                % dfdpi
                if flags(1)
                    gmms{i}.PComponents = gmms{i}.PComponents - grad_pi(i-1,:)*alph;
                    % Projection, Nonegativity.
                    gmms{i}.PComponents = max(gmms{i}.PComponents, 0); 
                    % Projection, Sum is 1.
                    %gmms{i}.PComponents = gmms{i}.PComponents/sum(gmms{i}.PComponents); 
                end

                % dfdm
                if flags(2)
                    gmms{i}.mu = gmms{i}.mu - grad_m(:,:,i-1)*alph;
                end
                % dfds
                if flags(3)
                    gmms{i}.Sigma = gmms{i}.Sigma - (grad_s(:,:,:,i-1))*alph;    
                    Sigmanorm = norm(gmms{i}.Sigma(:));
%                     if Sigmanorm > C2
%                         gmms{i}.Sigma = C2*gmms{i}.Sigma./Sigmanorm;    
%                     end
                end
            end
            fval_new = l2distGMMs(gmms)/npts;
            if isnan(fval_new) && alph >= minalph
                gmms = gmms_bak;
                alph = alph / 2;
                continue
            end
            if alph < minalph 
                gmms = gmms_bak;
                isdone = 1;
                break
            end
            if l2distGMMs(gmms)/npts > fval
                alph = alph / 2;
                gmms = gmms_bak;
            else

                alph = min(alph * 2, maxalph);
                updated = true;
                break
            end

        end
        if isdone ==1;
            fval_history(1:iter);
            break
        end
    end
    disp('Finished.')
    gnorm = [gnorm1, gnorm2];
    
    for i = 2:N-1
        % Projection, Nonegativity.
        gmms{i}.PComponents = max(gmms{i}.PComponents, 0); 
        gmms{i}.PComponents = gmms{i}.PComponents/sum(gmms{i}.PComponents); 
    end
    fval_history = [fval_history;l2distGMMs(gmms)/npts];
                
end
function val = question(statement, val1, val2)
    if statement
        val = val1;
    else
        val = val2;
    end
end

