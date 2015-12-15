function pd = getPrincipalDirection(data)
%GETPRINCIPALDIRECTION returns the principal direction of data. Data is a
%set of row vectors.
%
%
%   See Also: EIG, EIGS, ROTATION2D

%   $ Hyunwoo J. Kim $  $ 2015/04/01 12:38:16 (CDT) $

    data = data-repmat(mean(data,1), size(data,1),1);
    [U, D] = eig(data'*data);
    [~, idx] = max(diag(D));
    pd = U(:,idx);
end