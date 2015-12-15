function maxK = maxKgmms(gmms)
%MAXKGMMS returns the maximum components of gmms.
%
%
%   See Also: EAP_AFFINE_DEMO

%   $ Hyunwoo J. Kim $  $ 2015/04/13 17:38:33 (CDT) $
    Ks = zeros(1,length(gmms));
    for i = 1:length(gmms)
        if isempty(gmms{i})
            continue
        end
        Ks(i) = gmms{i}.NComponents;
    end
    maxK = max(Ks);
end