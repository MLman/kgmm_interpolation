function p = proj_M_spd(X,varargin)
%PROJ_M_SPD converts a matrix to a SPD matrix.
%
%
%   See Also:

%   $ Hyunwoo J. Kim $  $ 2015/01/28 19:49:31 (CST) $

    if ~myisreal(X)
        p = NaN;
        return 
    end
    if nargin == 2
        c = varargin{1};
    else 
        c = 1e-10;
    end
    
    dim = size(X,1);
    % Make a matrix symmetric positive definite.
    if norm(X-X') > eps
        X = (X+X')/2;
    end
    [V, D ] = eig(X);
    D = diag(D);
    p = zeros(size(X));
    for i =1:length(D)
        if D(i) > 0
            p = p + D(i)*V(:,i)*V(:,i)';
        end
    end
    
    % Now X is spd
    % Make psd matrix
    % NOTE I don't understand this weird code.
    
    a = 1e-3; 
    pnew = p;

    % NOTE isspd tolerance 1e-10
    %while ~isspd(pnew, 1e-10)
    while ~isspd(pnew, c)
        pnew = p + a*eye(dim);
        a = 2*a;
    end
    p = pnew;

end