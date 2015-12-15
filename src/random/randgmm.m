function gmm = randgmm(d, k, varagin)
%RANDGMM generate a random GMM.
%
%
%   See Also: GMDISTRIBUTION

%   $ Hyunwoo J. Kim $  $ 2015/12/12 13:16:56 (CLT) $

    mu = randn(k,d);
    P = rand(1,k);
    P = P./sum(P);

    Sigma = eye(d);
    df = ceil(1.5*d);
    
    if nargin < 3
        covtype = 'fullcov'; % 
    else 
        covtype = varargin{1};
    end

        
    if strcmpi(covtype,'fullcov') == 1
        S = zeros(d,d,k);
        for i = 1:k
            S(:,:,i) = wishrnd(Sigma,df)/df;
        end
    end

    if strcmpi(covtype,'diagcov') == 1
        S = zeros(1,d,k);  % Double check
        for i = 1:k
            S(:,:,i) = abs(randn(1,k)); 
        end
    end

    gmm = obj2structGMM(gmdistribution(mu, S, P));
end
