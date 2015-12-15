function type = covtype(mx)
%COVTYPE checks the type of Covariance matrix SIGMA.
%
%   COVTYPE(Sigma) : 1 full
%                  : 2 diag
%                  : 3 non-square matrix
%                  : -1 unknown
%   See Also: DIAG3D
    
%   $ Hyunwoo J. Kim $  $ 2015/03/26 14:59:05 (CDT) $
    dims = size(mx);
    if dims(1) >= 1 && dims(1) == dims(2)
        type = 1;
    elseif dims(1) == 1 && dims(1) < dims(2)
        type = 2;
    elseif dims(1) ~= dims(2)
        type = 3; 
    else
        type = -1;
    end
end