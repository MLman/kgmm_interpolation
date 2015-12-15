function gmmhat = interpolationGMMbyEM(gmms, w, K, varargin)
%INTERPOLATIONGMMBYEM
%
%   gmmhat = interpolationGMMbyEM(gmms, w, K)
%   gmmhat = interpolationGMMbyEM(gmms, w, K, scaling)
%
%   See Also: L2MEANGMMS

%   $ Hyunwoo J. Kim $  $ 2015/04/13 17:10:34 (CDT) $
    if isempty(gmms)
        gmmhat = [];
        return
    end
    gmml2bar = l2meanGMMs(gmms, w);

    if nargin >= 4 
       scaling = varargin{1};
       gmml2bar.Sigma = gmml2bar.Sigma * scaling;
    end
    option.tol = 1e-18;
    gmmhat = em_gmm_closest_full(gmml2bar, K, option);
   
    if nargin >= 4
       gmmhat.Sigma = gmmhat.Sigma/scaling;
    end
end
 