function [l, d, f] = l2distGMMs(gmms)
%L2DISTGMMS calculates the length of the path along GMMS.
%   
%   L2DISTGMMS calculates the sum of l2dist between gmm_i+1, gmm_i.
%   l is the length of the path. d is the l2dist between gmm_i+1, gmm_i.
%   If length(GMMS) == n, then d is a n-1 length vector.
%   f is a objective function value. The sum of the squared l2 distance,
%   i.e., f = sum(d.^2).
%
%   [l, d, f] = l2distGMMs(GMMS) 
%   [l, ~, ~] = l2distGMMs(GMMS) % Path length
%   [~, ~, f] = l2distGMMs(GMMS) % Function evaluation.
%
%   See Also: L2DISTGMM, L2NORMGMM, INNERPRODGMM

%   $ Hyunwoo J. Kim $  $ 2014/10/27 20:22:39 (CDT) $

    d = zeros(1,length(gmms)-1);
    for i = 1:length(gmms)-1
        d(i) = l2distGMM(gmms{i+1}, gmms{i});
    end
    l = sum(d);
    f = sum(d.^2);
end