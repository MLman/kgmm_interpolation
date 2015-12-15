function TX = MOAtransform(T, X)
%MOATRANSFORM transforms X. T is structure for the mixture of affine
%transformations. 
%
%
%   WARNING: X MUST be in the same space as T.gmm
%
%
%   See Also: GETAFFINETRANSGMM, GETAFFINETRANSGAUSS, AFFINETRANSGMM

%   $ Hyunwoo J. Kim $  $ 2015/04/02 02:08:12 (CDT) $

    assert(isobject(T.gmm),'GMM MUST be an object of gmdistribution.');
    TX = size(X);
    W = posterior(T.gmm,X);
    for i = 1:size(X,1)
        A = zeros(size(X,2));
        b = zeros(size(X,2),1);
        for j = 1:T.NComponents
           A = A + T.As(:,:,j)*W(i,j);
           b = b + T.bs(:,j)*W(i,j);
        end
        TX(i,:) = (A*X(i,:)' + b)';

%        [val, idx] = max(W(i,:));
%        TX(i,:) = (T.As(:,:,idx)*X(i,:)' + T.bs(:,idx))';
    end
end