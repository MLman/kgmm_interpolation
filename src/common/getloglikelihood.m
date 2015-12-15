function ll = getloglikelihood(D, gmm)
%GETLOGLIKELIHOOD calculates loglikelihood of data D.
%   D is a set of row vectors.
%   gmm is an distribution object or my structure.
%   See Also: PDF

%   $ Hyunwoo J. Kim $  $ 2015/02/24 17:00:00 (CST) $

    if ~isobject(gmm)
        gmm = struct2objGMM(gmm);
    end
    ll = 0;
    for i = 1:size(D,1)
        ll = ll + log(pdf(gmm, D(i,:)));
    end
end